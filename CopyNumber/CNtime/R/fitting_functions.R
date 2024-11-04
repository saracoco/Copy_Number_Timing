library(rstan)
library(tidyr)
library(ggthemes)
library(matrixStats)
library(dplyr)

#' fit_model_selection_best_K Function
#'
#' being able to perform model selection on already simulated data or real data and add initialization to the inference 
#' @param data data
#' @param karyo karyotype
#' @param purity sample purity
#' @param max_attempts max number of repeated inference for ADVI
#' @param INIT initialization values list (obtained by get_init)
#' @keywords fit
#' @export
#' @examples
#' fit_model_selection()


fit_model_selection_best_K = function(data, karyo, purity, max_attempts=4, INIT=TRUE, tollerance = 0.01, compute_external_metric = FALSE){

  if (INIT==TRUE){
    #tau_single_inference <- fit_single_segments(data, purity=purity) # uncomment to run initialization with single segment model
  }

#MODEL SELECTION
model_selection_tibble <- dplyr::tibble()

#  if (length(karyo) <= 15){
#     k_max = length(karyo)
#  } else { 
#     k_max = (length(karyo))/2  
#  }
 if (length(karyo) <= 2){
   k_max = (length(karyo)) 
 } else if (length(karyo) <= 7){
    k_max = (length(karyo)-1) 
 } else if (length(karyo) <= 15) { 
    k_max = ((floor(length(karyo)/2))-1)
 } else{
    k_max = ceiling(sqrt(length(karyo))) 
 }
  input_data <- prepare_input_data(data, karyo, 1, purity=purity)
  accepted_mutations = readRDS("results/accepted_mutations.rds") #potrei prenderli direttamente da input_data quando non farò model validation      
  S <- length(unique(accepted_mutations$segment_id))

  entropy_per_segment_matrix = matrix(0, k_max, S) # k_max rows and S columns
  entropy_per_segment_matrix_norm = matrix(0, k_max, S)
  
 for (K in 1:k_max) {
    input_data <- prepare_input_data(data, karyo, K, purity=purity)

    if (INIT==TRUE){
      # inits_chain <- get_init(tau_single_inference, K)     # uncomment to run initialization with single segment model
      inits_chain <- get_init_simpler(accepted_mutations, K, purity = purity)
      res <- fit_variational(input_data, max_attempts=max_attempts, initialization = inits_chain, INIT = TRUE, tollerance = tollerance)
    } else {
      res <- fit_variational(input_data, max_attempts=max_attempts, INIT = FALSE, tollerance = tollerance)
    }


    # Plot ELBO values over iterations
    output_files <- res$latent_dynamics_files()
    print(paste0("output_files ", output_files,"\n"))
    elbo_data <- read.csv(output_files, header = FALSE, comment.char = "#")
    colnames(elbo_data) <- c("iter", "time_in_seconds", "ELBO")
    iterations <- elbo_data$iter  # iteration column
    elbo_values <- elbo_data$ELBO  # ELBO column
    elbo_df <- data.frame(iteration = iterations, elbo = elbo_values)
    # saveRDS(elbo_df, paste0("./elbo_df_",K,".rds"))    
    p <- ggplot(elbo_df, aes(x = iteration, y = elbo)) +
    geom_line(color = "blue") +
    labs(title = "ELBO Values over Iterations",
        x = "Iteration",
        y = "ELBO")
    saveRDS(p, paste0("./elbo_vs_iterations_",K,".rds"))  


    if (compute_external_metric){
      compute_external_metric(accepted_mutations, res, K)   
    }

    ############################################
    # SELECT THE BEST MODEL 
    draws = 1000 

    stanfit <- rstan::read_stan_csv(res$output_files())
    S <- length(unique(accepted_mutations$segment_id))
    total_number_params <- K+((K-1)*S)+2                    # tau = K, w = K*S, phi, kappa (dirichlet reparametrization)
    N <- nrow(accepted_mutations)

    # Check log likelihood values 
    lp__ <- res$draws("lp__")
    lp_approx__ <- res$draws("lp_approx__")
    log_lik_matrix <- extract_log_lik(stanfit, parameter_name = "log_lik", merge_chains = TRUE) # merge chains potrebbe non servire
    
    log_lik_total_per_sample <- rowSums(log_lik_matrix)
    L <- median(log_lik_total_per_sample)

    BIC <- ((total_number_params * log(N)) - 2 * L) # %>% unname()
    AIC <- 2 * total_number_params - 2 * L

    loo_result<-loo(log_lik_matrix)
    loo_value <- loo_result$estimates[3, "Estimate"]  # LOO-CV estimate

    # Calculate ICL --> move outside the main function
    # Generate names for w[sim_params_num, K] format
    names_weights <- outer(1:input_data$S, 1:K, 
                FUN = function(i, j) paste0("w[", i, ",", j, "]"))
    names_weights <- as.vector(names_weights)
    w_ICL <- as.matrix(stanfit, pars = names_weights)
    dim(w_ICL) = c(draws,S*K)
    w_ICL <- apply(w_ICL, 2, median) # median check over draws
    dim(w_ICL) = c(S,K) # check by row
    w_ICL = t(w_ICL)

    num_mutations_all <- c()
    for (i in seq_along(unique(accepted_mutations$segment_id))) {
        segment <- unique(accepted_mutations$segment_id)[i]
        num_mutations_single <- nrow(accepted_mutations %>% filter(segment_id == segment))
        num_mutations_all <- c(num_mutations_all, num_mutations_single)
    }
    mut_per_seg = num_mutations_all

    res_entropy = 0
    post = w_ICL
    for (k in 1:K ){
        post_k = post[k]
        log_post_k = log(post_k + 0.000001)
        post_k_entr = post_k * log_post_k * mut_per_seg
        post_k_entr = sum(post_k_entr)
        post_k_entr = -1 * (post_k_entr)            
        res_entropy = res_entropy + post_k_entr
    }
    entropy = res_entropy
    ICL = BIC + entropy


    model_selection_tibble <- dplyr::bind_rows(model_selection_tibble, dplyr::tibble(K = K, BIC = BIC, AIC = AIC, LOO = loo_value, Log_lik = L, ICL = ICL))


    # ICL PER SEGMENT to see its behaviour with increasing K
    post_Segments = t(w_ICL)
    entropy_per_segment = c()
    entropy_per_segment_norm = c()
    for (s in 1:S ){
      post_s = post_Segments[s]
      log_post_s = log(post_s + 0.000001)
      post_s_entr = post_s * log_post_s * mut_per_seg
      post_s_entr = sum(post_s_entr)
      post_s_entr = -1 * (post_s_entr)            
      entropy_per_segment = c(entropy_per_segment, post_s_entr)

      post_s_entr_norm = post_s * log_post_s
      post_s_entr_norm = sum(post_s_entr_norm)
      post_s_entr_norm = -1 * (post_s_entr_norm)            
      entropy_per_segment_norm = c(entropy_per_segment_norm, post_s_entr_norm)

    }

    print("entropy per segment: ")
    print(entropy_per_segment)

    print("entropy per segment normalized: ")
    print(entropy_per_segment_norm)

    entropy_per_segment_matrix_norm[K,] = entropy_per_segment_norm
    entropy_per_segment_matrix[K,] = entropy_per_segment




    # CHECK PRIORS (function outside the main function)
    draws_matrix <- res$draws(format = "matrix")
    color_scheme_set("teal")
    #manage to get the K from the model fit directly rather than as input
    names_tau <- paste("tau_prior[", 1:K, "]", sep = "")
    areas_tau <- mcmc_areas(
        draws_matrix,
        pars = names_tau,
        prob = 0.8, # 80% intervals
        #prob_outer = 0.95, # 99%
        point_est = "median"
    )+
        labs(
        title = " Prior distributions",
        subtitle = "with median and 80% and 95% intervals"
        )+
        xlim(0, 1) + # + scale_x_continuous(breaks = c(1:5), labels = c("A", "B", "C", "D", "E"))
        theme(plot.background = element_rect(fill = "white"),text = element_text(size = 25))

    dim_1 = (8 + (input_data$S/2))
    ggsave(paste0("./plots/priors_tau_",K,".png"),  width = dim_1, height = dim_1, limitsize = FALSE, device = png, plot=areas_tau)


    intervals_weigths_per_tau <- list()
    for (k in 1:K){
            names_weights <- paste("w_prior[",1:input_data$S,",", k, "]", sep = "") 
            intervals_weigths_per_tau[[k]] <- mcmc_areas_ridges(draws_matrix, pars = names_weights, point_est = "median", prob = 0.8)+
                        labs(
                        title =  str_wrap( paste0("Prior distributions of the weigths for tau ",k), width = 30 + K + sqrt(input_data$S)),
                        subtitle = "with median and 80% and 95% intervals"
                        ) +
                        theme(plot.background = element_rect(fill = "white"),text = element_text(size = 15))

    }
    areas_w <- gridExtra::grid.arrange(grobs = intervals_weigths_per_tau, ncol=K) #add global title

    dim_2 = (15 + (input_data$S/2))
    ggsave(paste0("./plots/priors_w_",K,".png"),  width = dim_2, height = dim_1, limitsize = FALSE, device = png, plot=areas_w)



    # move these somewhere else
    saveRDS(res, paste0("results/res",K,".rds"))
    saveRDS(input_data, paste0("results/input_data_",K,".rds"))

    # plotting
    p <- plotting_fit(res,input_data, data,K)
    ggsave(paste0("./plots/plot_inference_",K,".png"), width = (12 + (input_data$S/2)), height = (16 + (input_data$S/2)), limitsize = FALSE, device = png, plot=p)

    plot_partition = plotting_cluster_partition(res, K, VALIDATION=TRUE)
    ggsave(paste0("./plots/inferred_partition_",K,".png"), width = 8, height = 5, limitsize = FALSE, device = png, plot=plot_partition)

    }
  
   # plot the entropy behaviour 
   if (k_max>2){
    print(paste0("k_max: ", k_max))
    plot_entropy <- plotting_entropy(entropy_per_segment_matrix, entropy_per_segment_matrix_norm, K)
    ggsave(paste0("./plot_entropy.png"),width = (6 + (K*2)), height = (5), plot = plot_entropy)
   }



   model_selection_tibble_temp <- model_selection_tibble[1:2, bycol= TRUE]
   best_K_temp <- model_selection_tibble_temp %>% dplyr::filter(BIC == min(BIC)) %>% pull(K)   

     if (best_K_temp!=1){
          if (k_max==2){
          best_K <- 2
          }else{
            while(mean(entropy_per_segment_matrix_norm[best_K_temp+1,]) - mean(entropy_per_segment_matrix_norm[best_K_temp,]) < 0 & best_K_temp < k_max ){
               best_K_temp = best_K_temp + 1
               if ( best_K_temp == k_max ){
                break
               }}}
      } else {
      best_K <- 1
      }
    best_K <- best_K_temp

   if(best_K==k_max){
    cli::cli_alert_info("The algorithm should be run with more Components ")
   }
   

   input_data <- prepare_input_data(data, karyo, best_K, purity=purity, VALIDATION == FALSE)
   
   if (INIT==TRUE){
    accepted_mutations = readRDS("results/accepted_mutations.rds")
    inits_chain <- get_init_simpler(accepted_mutations, best_K, purity = purity)
    #  inits_chain <- get_init(tau_single_inference, best_K)
     print(paste0("These are the values used for initializing the model ",inits_chain))
     res <- fit_variational(input_data, max_attempts=max_attempts, initialization = inits_chain, INIT=TRUE, tollerance = tollerance)
   } else {
     res <- fit_variational(input_data, max_attempts=max_attempts, INIT=FALSE, tollerance = tollerance)
   }
   
    if (compute_external_metric){
    compute_external_metric(accepted_mutations, res, best_K, best_K = TRUE)   
    }

   p_best_K <- plotting_fit(res,input_data, data, best_K)
   ggsave(paste0("./plots/plot_inference_",best_K,"_best_K.png"), width = (12 + (input_data$S/2)), height = (16 + (input_data$S/2)), limitsize = FALSE, device = png, plot=p_best_K)    

   p_elbo_iter <- plotting_elbo(k_max)
   ggsave(paste0("./elbo_vs_iterations_.png"),width = (20 + (input_data$S/2)), height = (10 + (input_data$S/2)), plot = p_elbo_iter)

  return(list(data = data, model_selection_tibble = model_selection_tibble, res_best_K=res, best_K=best_K, input_data=input_data, accepted_mutations=accepted_mutations
))

}












#' fit_variational Function
#'
#' This function performs the inference for CN timing using the ADVI algorithm.
#' @param input_data data
#' @param max_attempts max number of repeated inference for ADVI
#' @param init list of initializing values
#' @param INIT boolean value if to initialize the inference with specific values
#' @param initial_iter description
#' @param grad_samples description
#' @param elbo_samples description
#' @keywords fit
#' @export
#' @examples
#' fit_variational()

##### fit without filtering step ############################################################################################################################
library(cmdstanr)

#fit variational taking the best run 
fit_variational <- function(input_data, max_attempts = 5, initialization = NULL, INIT = TRUE, initial_iter = 1000, grad_samples = 10, elbo_samples = 100, tollerance = 0.01) {
  # Load the Stan model
  model <- cmdstanr::cmdstan_model("../../CopyNumber/models/timing_mixed_simple.stan")
  best_elbo <- -Inf  # Start with the worst possible ELBO
  best_fit <- NULL  # To store the best model fit
  total_attempts <- 0
  
  for (attempt in 1:max_attempts) {
    message("Attempt ", attempt, " of ", max_attempts)
    
    # retries <- 0  # Initialize retries for each attempt
    fit_successful <- FALSE  # Reset success flag for each attempt
    iter <- initial_iter  # Reset iterations for each new attempt
    
    retries <- 0

    while (!fit_successful && retries < 3) {  # Set a max number of retries for each attempt

      # message("Running inference, retry ", retries + 1)

      # Increment total_attempts at the start of each retry
      total_attempts <- total_attempts + 1
      retries <- retries + 1

      result <- tryCatch({
        # Attempt variational inference
        # jacobian: the default is FALSE, meaning optimization yields the (regularized) maximum likelihood estimate. 
        # Setting it to TRUE yields the maximum a posteriori estimate.
        if (INIT == TRUE) {
          res <- model$variational(
            data = input_data, 
            init = list(initialization),  # Use the provided initialization
            iter = iter, 
            grad_samples = grad_samples, 
            elbo_samples = elbo_samples,
            save_latent_dynamics = TRUE,
            draws = 1000,
            # output_dir = "./",
            eval_elbo = 1,
            tol_rel_obj = tollerance 
          )
          print(res$init())


          # Shuffle and moderately perturb initialization values
          ###############################################################################################################################################################
          shuffled_taus <- sample(initialization$tau)
          # print(paste0("init_tau: ", initialization$tau ))
          perturbed_taus <- shuffled_taus + rnorm(length(shuffled_taus), 0, 0.1) # Adjust the perturbation scale as needed
          perturbed_taus[perturbed_taus >= 0.88] <- 0.88  # Check for elements greater than 1 and replace them with 1 otherwise fit fails
          perturbed_taus[perturbed_taus <= 0] <- 0.07  

          initialization$tau = perturbed_taus
          # phi does not change

          # Apply the perturbation to w as wells
          print (paste0("w_init after inference: ",initialization$w," ",ncol(initialization$w) == 1))
          #######################################################################################################################################################################
          

        } else {
          res <- model$variational(
            data = input_data, 
            iter = iter, 
            grad_samples = grad_samples, 
            elbo_samples = elbo_samples,
            save_latent_dynamics = TRUE,
            draws = 1000,
            # output_dir = "./",
            eval_elbo = 1,
            tol_rel_obj = tollerance 
          )
        }

        output_files <- res$latent_dynamics_files()
        elbo_data <- read.csv(output_files, header = FALSE, comment.char = "#")
        colnames(elbo_data) <- c("iter", "time_in_seconds", "ELBO")
        elbo_values <- elbo_data$ELBO  # The ELBO column
        elbo <- elbo_values[length(elbo_data)]
        # Update the best model if this ELBO is better
        if (!is.na(elbo) && elbo > best_elbo) {
          best_elbo <- elbo
          best_fit <- res
        }

        message("ELBO for this run: ", elbo)

        # Check for invalid log evaluations
        output_files <- res$output_files()
        fit_successful <- TRUE  # Mark this fit as successful
        
      }, error = function(e) {
        message("An error occurred during inference: ", e$message)
        # retries <- retries + 1  # Increment retry count
        fit_successful <- FALSE  # Mark fit as unsuccessful
        NULL  # Ensure NULL is returned so loop can continue
      })
      
      # No return from `tryCatch` - just continue the loop if not successful
    }

    if (retries == 3) {
      message("Max retries reached for attempt ", attempt)
    }
  }
  
  if (!is.null(best_fit)) {
    # message("Best ELBO after ", total_attempts, " attempts: ", best_elbo)
    return(best_fit)  # Return the result with the best ELBO
  } else {
    message("Inference could not be completed successfully.")
    return(NULL)  # Return NULL if no fit was successful
  }
}


#' prepare_input_data Function
#'
#' This function obtains the list of input data to be used for CN timing.
#' @param data  simulated data
#' @param karyo karyotype posso togliere lo ricavo da data XXXXXXX
#' @param K number of components
#' @param purity sample purity
#' @keywords input
#' @export
#' @examples
#' prepare_input_data()

prepare_input_data = function(data, karyo, K, purity, VALIDATION = TRUE){
  
  data <- data[order(data$segment_id), ]

  karyotype <- data %>%
  group_by(segment_id) %>%
  summarise(karyotype = unique(karyotype)[1]) %>%
  pull(karyotype)

  seg <- data %>%
  group_by(segment_id) %>%
  summarise(karyotype = unique(karyotype)[1]) %>%
  pull(segment_id)

  peaks <- matrix(0, nrow = length(karyotype), ncol = 2)
  for (i in 1:length(karyotype)) {
    peaks[i,] <- get_clonal_peaks(karyotype[i], purity=purity)
  }

  cat("these are the peaks used to filter out mutations \n")
  print(peaks)
  cat("These are the Karyotypes extracted from the data, used in peaks\n")
  print(karyotype)
  print(seg)
  cat("These are the Karyotypes directly from the simulation, used for peaks in the model\n ")
  print(karyo)

  alpha = 0.05
  min_mutations_number = 2
  accepted_mutations <- data.frame()
  
  if (nrow(data) > 0) {
    # Check if mutation is inside CI
    probs <- c(alpha/2, 1 - alpha/2)
    
    DP <- data$DP
    NV <- data$NV

    print(paste0("DP ", DP[1:4]))
    

    # require(reshape2)
    # mat$id <- rownames(mat) 
    # melt(mat)
    
    accepted_idx <- lapply(1:length(DP), function(i) {
      for (j in unique(data$segment_id)) {
          if (data$segment_id[i]==j){
          quantiles_1 <- qbinom(probs, DP[i], peaks[j,1])
          quantiles_2 <- qbinom(probs, DP[i], peaks[j,2])
        if (((NV[i] >= quantiles_1[1]) && (NV[i] <= quantiles_1[2])) | ((NV[i] >= quantiles_2[1]) && (NV[i] <= quantiles_2[2]))) {
          return(i)
        }
        }
      }
    }) %>% unlist()
    
    print(paste0(data$segment_name_real[accepted_idx]))
    # Get only good mutations
      accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx], segment_id=data$segment_name[accepted_idx], karyotype=data$karyotype[accepted_idx], tau=data$tau[accepted_idx] , segment_index=data$segment_id[accepted_idx])
    
  }
  
  counts <- table(accepted_mutations$segment_id)
  minimum_number_per_segment <- all(counts >= min_mutations_number)
  
  # if (minimum_number_per_segment == TRUE) {
  
  input_data <- list(
    S = length(unique(data$segment_name[accepted_idx])),
    K = K,
    N = nrow(accepted_mutations),
    karyotype = lapply(karyo, karyo_to_int) %>% unlist(),
    seg_assignment = data$segment_id[accepted_idx],
    peaks = peaks, #peaks_inference(karyo,purity),
    NV = NV[accepted_idx],
    DP = DP[accepted_idx]
  )
  
  # cat("These are the peaks which are passed in the inference task, more reliable let's say\n")
  # print(input_data$peaks)
  # print(karyotype)

  print("saving_RDS")
  saveRDS(accepted_mutations, paste0("results/accepted_mutations.rds"))

  return(input_data)
# } else {         cli::cli_alert_info("Segment with index {.val {segment_idx}} has less that ")}
}




#' get_init Function
#'
#' Perform cmeans on the "single inference" result and retrieve the centroids and u matrix to be used as initialization parameteres
#' @param K number of clusters
#' @param data dataframe of points to be clustered
#' @param phi parameters of reparametrization initialized uninformatively
#' @param kappa parameters of reparametrization
#' @keywords cluster
#' @export
#' @examples
#' get_init()

get_init = function(data, K, phi=c(), kappa=5){

  myReps <- function(x, y, n) rep(x, (x %in% y) * (n-1) +1)
  phi = myReps(1/K, 1/K, K)

  res.fcm <- fcm(data, centers=K)

  init_taus <- c(res.fcm$v)
  init_taus[init_taus >= 1] <- 0.99  # Check for elements greater than 1 and replace them with 1 otherwise fit fails
  init_taus[init_taus == 0] <- 0.01   
  init_w <- as.matrix(res.fcm$u)
  epsilon <- 1e-4
  perturbed_probabilities <- init_w + epsilon

  # Rinormalizza i pesi lungo l'asse 1 (per riga)
  normalized_probabilities <- (apply(perturbed_probabilities, 1, function(x) x / sum(x)))
  # Verifica che i pesi siano stati rinormalizzati correttamente
  #print(rowSums(normalized_probabilities))
  # Stampa i pesi normalizzati
  init_w <- (normalized_probabilities)
  if (K==1){
    init_w = t(init_w)
  }

  inits_chain <- list(w = t(init_w), tau = init_taus, phi=phi, kappa=kappa)

  cat(paste0("w = ",inits_chain$w,"\n "))
  cat(paste0("tau = ", inits_chain$tau, "\n "))
  cat(paste0("phi = ", inits_chain$phi, "\n "))
  cat(paste0("kappa = ", inits_chain$kappa))


  return(inits_chain)
}






#' get_init_simpler
#'
#' Perform cmeans on the "single inference" result and retrieve the centroids and u matrix to be used as initialization parameteres
#' @param K number of clusters
#' @param data dataframe of points to be clustered
#' @param phi parameters of reparametrization initialized uninformatively
#' @param kappa parameters of reparametrization
#' @keywords cluster
#' @export
#' @examples
#' get_init()

#take as input the data ready to be used for inference
get_init_simpler = function(accepted_mutations, K, phi=c(), kappa=5, purity){

    accepted_mutations <- accepted_mutations[order(accepted_mutations$segment_index), ]

      karyotype <- accepted_mutations %>%
      group_by(segment_index) %>%
      summarise(karyotype = unique(karyotype)[1]) %>%
      pull(karyotype)

      seg <- accepted_mutations %>%
      group_by(segment_index) %>%
      summarise(karyotype = unique(karyotype)[1]) %>%
      pull(segment_index)

      peaks <- matrix(0, nrow = length(karyotype), ncol = 2)
      for (i in 1:length(karyotype)) {
        peaks[i,] <- get_clonal_peaks(karyotype[i], purity=purity)
      }


      alpha = 0.05
      min_mutations_number = 2
      
      if (nrow(accepted_mutations) > 0) {
        # Check if mutation is inside CI
        probs <- c(alpha/2, 1 - alpha/2)
        
        DP <- accepted_mutations$DP
        NV <- accepted_mutations$NV

        
        print(unique(accepted_mutations$segment_index))
        alpha_beta_all <- lapply(1:length(DP), function(i) {
          for (j in unique(accepted_mutations$segment_index)) {
              if (accepted_mutations$segment_index[i]==j){
              quantiles_1 <- qbinom(probs, DP[i], peaks[j,1])
              quantiles_2 <- qbinom(probs, DP[i], peaks[j,2])

            if ((NV[i] >= quantiles_1[1]) && (NV[i] <= quantiles_1[2])) {
                return("omega1")
            }else if ((NV[i] >= quantiles_2[1]) && (NV[i] <= quantiles_2[2])){
                return("omega2")

            }
            }
          }
        }) %>% unlist()
    

    print(paste0("alpha_beta_all ",length(alpha_beta_all)))
    df <- data.frame(alpha_beta_all = alpha_beta_all, segment_id = accepted_mutations$segment_index, karyotype = accepted_mutations$karyotype)
    # Group by segment_id and calculate the proportions of alpha and beta
    proportions <- df %>%
      group_by(segment_id, karyotype) %>%
      summarise(
        proportion_alpha = mean(alpha_beta_all == "omega1"),
        proportion_beta = mean(alpha_beta_all == "omega2"),
        .groups = 'drop'
      )

  

  myReps <- function(x, y, n) rep(x, (x %in% y) * (n-1) +1)
  phi = myReps(1/K, 1/K, K)


  # 1) calcolo prima i tau e poi clusterizzo (perché dipendono dal karyotypo)
  # Initialize the tau_posterior column
  proportions$tau_posterior <- NA  # Create an empty column in the dataframe

  # Loop through each row to calculate tau_posterior
  for (i in 1:nrow(proportions)) {
    if (proportions$karyotype[i] == '2:1') {
      proportions$tau_posterior[i] <- 3 * (proportions$proportion_beta[i]+0.0001) / 
                                      (2 * (proportions$proportion_beta[i]+0.0001) + (proportions$proportion_alpha[i]+0.0001))
    } else {
      proportions$tau_posterior[i] <- 2 * (proportions$proportion_beta[i]+0.0001) / 
                                      (2 * (proportions$proportion_beta[i]+0.0001) + (proportions$proportion_alpha[i]))
    }
  }

  # Check if tau_posterior contains any NA values
  if (any(is.na(proportions$tau_posterior))) {
    stop("tau_posterior contains NA values, check the loop calculations.")
  }

  # Apply fuzzy c-means clustering
  res.fcm <- fcm(as.matrix(proportions$tau_posterior), centers = K)

  print(paste0("dim(proportions$tau_posterior): ", length(proportions$tau_posterior) ))
  if (length(proportions$tau_posterior)==K){
      init_taus <- (proportions$tau_posterior)  
    }else{
      init_taus <- c(res.fcm$v)
    }
  cat(paste0("init_taus from clustering  ",init_taus))

  # init_taus <- c(res.fcm$v)
  init_taus[init_taus >= 0.88] <- 0.88  # Check for elements greater than 1 and replace them with 1 otherwise fit fails
  init_taus[init_taus == 0] <- 0.07   
  init_w <- as.matrix(res.fcm$u)
  epsilon <- 1e-4
  if (all(init_w > 0.5) ){
    perturbed_probabilities <- init_w - epsilon
  } else {
    perturbed_probabilities <- init_w + epsilon
  }
  

  # Rinormalizza i pesi lungo l'asse 1 (per riga)
  normalized_probabilities <- (apply(perturbed_probabilities, 1, function(x) x / sum(x)))
  # Verifica che i pesi siano stati rinormalizzati correttamente
  #print(rowSums(normalized_probabilities))
  # Stampa i pesi normalizzati
  init_w <- (normalized_probabilities)
  if (K==1){
    init_w = t(init_w)
  }

  inits_chain <- list(w = t(init_w), tau = init_taus, phi=phi, kappa=kappa)

  cat(paste0("w = ",inits_chain$w,"\n "))
  cat(paste0("tau = ", inits_chain$tau, "\n "))
  cat(paste0("phi = ", inits_chain$phi, "\n "))
  cat(paste0("kappa = ", inits_chain$kappa))

  return(inits_chain)
}

}







#' peaks_inference Function
#'
#' This function obtains the peaks to be used for CN timing model from a list of karyotypes (used in prepare_input_data) (can be sub with lapply)
#' @param karyo karyotype
#' @param purity peaks
#' @keywords peaks
#' @export
#' @examples
#' peaks_inference()

peaks_inference = function(karyo, purity){
  peaks <- matrix(0, nrow = length(karyo), ncol = 2)
  for (i in 1:length(karyo)) {
    peaks[i,] <- get_clonal_peaks(karyo[i], purity=purity)
  }
  return(peaks)
}








#' karyo_to_int Function
#'
#' This function allows you to change the string representation of the CN events into numeric values 0 -1.
#' @param k karyotype string
#' @keywords karyotype
#' @export
#' @examples
#' karyo_to_int()

karyo_to_int <- function(k) {
  if (k == "2:1") return(1)
  return(0)
}


#' get_clonal_peaks Function
#'
#' This function allows you to obtain the theoretical peaks that would be observed in a VAF spectrum.
#' @param k karyotype
#' @param purity peaks
#' @keywords peaks
#' @export
#' @examples
#' get_clonal_peaks()

get_clonal_peaks = function(k, purity) {
  multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
  major <- multiplicities[1]
  n_tot <- sum(multiplicities)
  # get only Major and 1
  multiplicities <- c(1, major)
  peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
  return(sort(peaks))
}



#' get_tau_posteriors Function
#'
#' This function allows you to obtain the posterior for the single segment inference.
#' @param fit fit results
#' @param k karyotype
#' @keywords tau
#' @export
#' @examples
#' get_tau_posteriors()

get_tau_posteriors = function(fit, k) {
  omega1 <- fit$draws("omega[1]", format = "matrix") %>% as.numeric()
  omega2 <- fit$draws("omega[2]", format = "matrix") %>% as.numeric()
  if (k == '2:1') {
    tau_posterior <- 3 * omega2 / (2*omega2 + omega1)
  } else {
    tau_posterior <- 2 * omega2 / (2*omega2 + omega1)
  }
  tau_posteriors <- dplyr::tibble(tau = tau_posterior)
  tau_posteriors
}




# mode <- function(v) {
#    uniqv <- unique(v)
#    uniqv[which.max(tabulate(match(v, uniqv)))]
# }


















#' fit_single_segments Function
#'
#' -NEED TO ADD FILTERING STEP - This function allows you to obtain simulated data and perform the single segments inference for CN timing, could be used to obtain initialization values for tau and w.
#' @param data data
#' @param alpha confidece level
#' @param purity sample purity
#' @keywords fit
#' @export
#' @examples
#' fit_single_segments()

fit_single_segments = function(data, alpha = .05, purity ){
  # FIRST OPTION: run single segment inference for initialization of tau and w
  model_single <- cmdstanr::cmdstan_model("../../CopyNumber/models/mixture_CNA_timing_binomial.stan")

  #prepare data for single segment inference : S = number of segments
  S <- length(unique(data$segment_id))
  data_plot <- dplyr::tibble()
  inference_results <- dplyr::tibble()
  summarized_results <- dplyr::tibble()

  plots = c()
  tau = c()
  inference_results_tot = c()

  # fit the model for each segment
  for (segment_index in 1:S){

    print(segment_index)

    # filtering step

    data_single = data %>% filter(segment_id==segment_index)
    k = unique(data$karyotype[data$segment_id==segment_index])
    peaks_single <- get_clonal_peaks(k, purity=purity)

    probs <- c(alpha/2, 1 - alpha/2)

    DP <- data_single$DP
    NV <- data_single$NV

    accepted_idx <- lapply(1:length(data_single$DP), function(i) {
      for (p in peaks_single) {
        
        quantiles <- qbinom(probs, data_single$DP[i], p)

        if ((data_single$NV[i] >= quantiles[1]) && (data_single$NV[i] <= quantiles[2])) {
          return(i)
        }
      }
    }) %>% unlist()

    # Get only good mutations
    accepted_mutations <- data.frame(DP = data_single$DP[accepted_idx], NV = data_single$NV[accepted_idx])



    input_data <- list(
      N = length(data_single$NV[accepted_idx]),
      NV = data_single$NV[accepted_idx], #data$NV[data$segment_id==segment_index],
      DP = data_single$DP[accepted_idx], #data$DP[data$segment_id==segment_index],
      peaks = peaks_single
    )







    fit <- model_single$sample(data=input_data, iter_warmup=2000, iter_sampling=2000, chains=4, parallel_chains=4)


    # Compute tau posteriors
    tau_posteriors <- get_tau_posteriors(fit, k)$tau %>% unname() %>% as.numeric()

    q1 <- alpha / 2
    q2 <- 1 - alpha / 2
    tau_low <- quantile(tau_posteriors, q1) %>% unname()
    tau_high <- quantile(tau_posteriors, q2) %>% unname()
    tau_mean <- mean(tau_posteriors)

    inference_results <- dplyr::bind_rows(inference_results, dplyr::tibble(tau = tau_posteriors, segment = unique(data$segment_id[data$segment_id==segment_index]), karyotype = k))
    summarized_results <- dplyr::bind_rows(summarized_results, dplyr::tibble(tau_low = tau_low, tau_mean = tau_mean, tau_high = tau_high, segment = unique(data$segment_id[data$segment_id==segment_index]), karyotype = k))

    inference_results_tot = c(inference_results_tot, inference_results)
    data_plot <- dplyr::bind_rows(data_plot, dplyr::tibble(tau = tau_posteriors, segment = unique(data$segment_id[data$segment_id==segment_index]), karyotype = k))

  }


  inference_single_segment <- data_plot %>%
    ggplot(mapping = aes(x=tau, fill=as.factor(segment))) +
    geom_histogram(alpha=.5, position = "identity", bins = 100) +
    labs(title = "Posterior distribution of tau obtained from single segment inference")

  #plot_data <- plot_data + facet_wrap(vars(karyotype, segment))

  #posterior distribution of tau obtained from the single segment inference
  inference_single_segment
  ggsave("./plots/inference_single_segments.png", width = 8, height = 8,  device = png)

  tau_single_inference <- summarized_results$tau_mean
  return(tau_single_inference)
}








compute_external_metric <- function(accepted_mutations, res, K, best_K = FALSE){

  # EXTERNAL METRIC 
    #obtain score of simulation accuracy EXTERNAL METRICS (with known ground truth)
    #MAE evaluated on the tau values associated to each segment compare original with predicted/arrigned
    #RI compare the original partition of segments with respect to the tau they are associated with and the partition obtained from the tau assigned through the model
    
    # accepted_mutations <- readRDS("results/accepted_mutations.rds") #Reload as I modify it for visualization purposes. Give as input do not call it from inside the function for the package
    
    names_tau <- paste("tau[", 1:K, "]", sep = "")
    tau_inferred <- res$draws(names_tau, format = "matrix")
    tau_inferred_map <- lapply(1:ncol(tau_inferred), function(i) {median(tau_inferred[,i])} ) %>% unlist() 

    all_differences = c()
    identity_matrix_RI = c()
    for (i in seq_along(unique(accepted_mutations$segment_id))) {
        segment <- unique(accepted_mutations$segment_id)[i] #potrei mettere direttamente i 
        seg_modified <- segment %>%  str_split(" ")
        segment_number <- seg_modified[[1]][2]
        tau_original <- unique(accepted_mutations %>% filter(segment_id == segment)%>% select(tau)) %>% pull(tau)
        
        names_weights <- paste("w[",segment_number,",", 1:K, "]", sep = "")  #regex_pars = c("w")
        weights_inferred <- res$draws(names_weights, format = "matrix")
        weights_inferred_median <- lapply(1:ncol(weights_inferred), function(i) {median(weights_inferred[,i])} ) %>% unlist() 

        tau_index_assigned <- which.max(weights_inferred_median)
        tau_inferred_assigned =  tau_inferred_map[which.max(weights_inferred_median)]

        difference <- abs(tau_inferred_assigned - tau_original) # MAE - meglio MSE?
        all_differences <- c(all_differences, difference)
        identity_matrix_RI = c(identity_matrix_RI, tau_index_assigned) #extract the vector of taus (unordered) to which the ordered segments are assigned
    }

    MAE <- mean(all_differences)
    
    if (best_K == TRUE){
      saveRDS(MAE, paste0("results/MAE_best_K.rds"))
    } else {
      saveRDS(MAE, paste0("results/MAE_",K,".rds"))
    }

    real_assignment <- accepted_mutations %>%
                                  group_by(segment_id) %>%
                                  summarize(tau = first(tau)) %>%
                                  arrange(match(segment_id, unique(accepted_mutations$segment_id))) %>%
                                  pull(tau) # here real_assignment is the true tau associated to the segment
    real_assignment <- rank(real_assignment, ties.method = "min") # here is the "class" to which it is associated --> the first event, the second ecc

    model_assignment <- identity_matrix_RI        #i,2,3 ecc 
    RI <- rand.index(real_assignment,model_assignment)
    #ARI <- adj.rand.index(real_assignment,model_assignment)  #NaN
    #saveRDS(ARI, paste0("results/ARI_",K,".rds"))
    if (best_K == TRUE){
          saveRDS(RI, paste0("results/RI_best_K.rds"))
    } else {
          saveRDS(RI, paste0("results/RI_",K,".rds"))
    }

    MAE = round(MAE, 3)
    print(paste0("MAE: ", MAE))
    RI = round(RI, 3)
    print(paste0("RI: ", RI))


}






