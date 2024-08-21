
#' get_taus_karyo Function
#'
#' This function allows you to obtain simulated data for copy number events (tau and karyo input for get_simulation or simulate_and_fit) inputting the number of events and the proportions of taus and karyos.
#' @param number_events number of events to be simulated
#' @param vector_taus unique tau to be present in the simulation
#' @param vector_karyo unique karyotype to be present in the simulation
#' @param weignths_taus proportion of tau to generate for each unique value
#' @param weights_karyo proportion of karyo to generate for each unique karyo
#' @keywords simulation
#' @export
#' @examples
#' get_taus_karyo()

get_taus_karyo = function (number_events,
                                vector_tau, vector_karyo,
                                weignths_tau, weights_karyo ){

  taus <- sample( vector_tau[1:length(vector_tau)], number_events, replace=TRUE, prob=weignths_tau )
  karyo <- sample( vector_karyo[1:length(vector_karyo)], number_events, replace=TRUE, prob=weights_karyo )

  return(list(taus=taus,karyo=karyo))

}


#' get_simulation Function
#'
#' This function allows you to obtain simulated data for copy number events.
#' @param taus tau
#' @param karyo karyotype
#' @param purity sample purity
#' @keywords simulation
#' @export
#' @examples
#' get_simulation()

get_simulation = function (taus, karyo, purity = 1){
  #taus = a vector of doubles between [0,1]
  #karyo = a vector of the same length of taus of strings of the type "2:1", "2:0", "2:2"

  S = length(karyo)
  names <- paste("segment", 1:S,sep = " ")
  all_sim <- dplyr::tibble()


  for (j in 1:length(taus)) {
    tau <- taus[j]
    print(tau)
    sim <- simulate_mutations(karyo[j], 20, tau = tau, l = 1e7, mu = 1e-4, 1e-2, segment_id = j)
    sim <- add_DP_and_NV(karyo[j], sim, coverage = 100, purity = purity)

    sim$tau = tau
    sim$j = names[j]
    all_sim <- dplyr::bind_rows(all_sim, sim)
  }

  peaks <- matrix(0, nrow = length(karyo), ncol = 2)
  for (i in 1:length(karyo)) {
    peaks[i,] <- get_clonal_peaks(karyo[i], purity)
  }


  return(all_sim)
}





#' simulate_mutations Function
#'
#' This function allows you to obtain the posterior for the single segment inference.
#' @param karyotype karyotype
#' @param time_interval time
#' @param tau tau generated
#' @param l l
#' @param mu mu
#' @param w w
#' @param segment_id segment id
#' @keywords mutations
#' @export
#' @examples
#' simulate_mutations()


simulate_mutations = function(karyotype, time_interval, tau, l, mu, w, segment_id = "segment_id") {
  available_karyotypes = c('2:0', '2:1', '2:2')
  # Check input
  if (!(time_interval > 0)) stop("time_interval must be a positive real number!")
  if (!(dplyr::between(tau, 0, 1))) stop("tau must be in [0,1]!")
  if (!(l > 0)) stop("l, which is the length of the segment considered, must be a positive real number!")
  if (!(mu > 0)) stop("mu, the mutation rate, must be a positive real number!")
  if (!(w > 0)) stop("w, the cell division rate, must be a positive real number!")
  if (!(karyotype %in% available_karyotypes)) stop("input karyotype is not available!")

  # Create time intervals
  t1 = time_interval * tau
  t2 = time_interval - t1

  # Extract nA and nB
  multiplicities = strsplit(karyotype, ":") %>% unlist() %>% as.numeric()
  nA = multiplicities[1]
  nB = multiplicities[2]

  # Init vector of used mutations
  used = c()
  all_mutations_df <- data.frame(mutation = c(), allele = c(), type = c())

  generate_mutations <- function(allele, l, mu, w, dt, type) {

    rate <- l * mu * w * dt
    n_mutations <- rpois(1, lambda = rate)

    if (n_mutations >= 1) {
      mutations <- lapply(c(1:n_mutations), function(i) {
        repeat{
          id = sample(LETTERS, 8, replace = TRUE) %>% paste(collapse = "")
          if(!id %in% used ) break
        }
        used <<- c(used, id) # <<- changes global value
        id
      })

      mutations_df <- data.frame(
        mutation = mutations %>% unlist(),
        allele = allele,
        type = type
      )
      return(mutations_df)
    } else {
      return(data.frame(mutation = c(), allele = c(), type = c()))
    }
  }

  # Simulate mutations on A
  if (nA == 2) {
    m1 <- generate_mutations("A1A2", l, mu, w, t1, "shared")
    m2 <- generate_mutations("A1", l, mu, w, t2, "private")
    m3 <- generate_mutations("A2", l, mu, w, t2, "private")

    all_mutations_df <- bind_rows(all_mutations_df, m1, m2, m3)
  }

  # Simulate mutations on B
  if (nB == 2) {
    m1 <- generate_mutations("B1B2", l, mu, w, t1, "shared")
    m2 <- generate_mutations("B1", l, mu, w, t2, "private")
    m3 <- generate_mutations("B2", l, mu, w, t2, "private")

    all_mutations_df <- bind_rows(all_mutations_df, m1, m2, m3)
  } else if (nB == 1) {
    m1 <- generate_mutations("B1", l, mu, w, time_interval, "private")
    all_mutations_df <- bind_rows(all_mutations_df, m1)
  }

  if (!(nrow(all_mutations_df) > 0)) stop("no mutations were simulated with the given parameters!")
  all_mutations_df <- all_mutations_df %>% mutate(karyotype = karyotype, segment_id = segment_id)
  all_mutations_df
}



#' add_DP_and_NV Function
#'
#' This function allows you to simulate the NV and DP for the CN events.
#' @param karyotype karyotype
#' @param simulated_mutations mutations
#' @param coverage coverage
#' @param purity purity
#' @keywords NV
#' @export
#' @examples
#' simulate_mutations()


add_DP_and_NV = function(karyotype, simulated_mutations, coverage, purity) {
  multiplicities = strsplit(karyotype, ":") %>% unlist() %>% as.numeric()
  nA = multiplicities[1]
  nB = multiplicities[2]

  peak <- function(m, n, purity) {
    num <- m * purity
    den <- (1 - purity) * 2 + n * purity
    return(num/den)
  }
  # Simulate DP and NV values for shared mutation
  p <- peak(nA, nA + nB, purity)
  shared_mutations <- simulated_mutations %>% filter(type == "shared")

  n <- nrow(shared_mutations)
  shared_mutations$DP <- rpois(n, coverage)
  shared_mutations$NV <- rbinom(n, size = shared_mutations$DP, prob = p)

  # Simulate DP and NV values for private mutation
  p <- peak(1, nA + nB, purity)
  private_mutations <- simulated_mutations %>% filter(type == "private")

  n <- nrow(private_mutations)
  private_mutations$DP <- rpois(n, coverage)
  private_mutations$NV <- rbinom(n, size = private_mutations$DP, prob = p)

  return(bind_rows(private_mutations, shared_mutations))
}
