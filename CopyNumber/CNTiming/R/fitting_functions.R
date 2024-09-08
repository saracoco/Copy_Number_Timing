

#' fit_model_selection_best_K Function
#'
#' being able to perform model selection on already simulated data or real data and add initialization to the inference 
#' @param all_sim data
#' @param karyo karyotype
#' @param purity sample purity
#' @param max_attempts max number of repeated inference for ADVI
#' @param INIT initialization values list (obtained by get_init)
#' @keywords fit
#' @export
#' @examples
#' fit_model_selection()

fit_model_selection_best_K = function(all_sim, karyo, purity=0.99, max_attempts=10, INIT=TRUE, simulation_params = simulation_params){

  
  if (INIT==TRUE){
    tau_single_inference <- fit_single_segments(all_sim, purity=purity)
  }


  #MODEL SELECTION
  model_selection_tibble <- dplyr::tibble()

  
  #Regola empirica: In generale, si può considerare come massimo plausibile un numero di cluster pari a 
  # sqrt(n/2), dove n è il numero di punti dati.
  

  if (length(karyo) <= 15){
    k_max = length(karyo)/2
  } else { k_max = sqrt(length(karyo))
    
  }
  
  for (K in 1:k_max) {
    input_data <- prepare_input_data(all_sim, karyo, K, purity)

    if (INIT==TRUE){
      inits_chain <- get_init(tau_single_inference, K)
      res <- fit_variational(input_data, max_attempts=max_attempts, initialization = inits_chain, INIT = TRUE)
    } else {
      res <- fit_variational(input_data, max_attempts=max_attempts, INIT = FALSE)

    }


    accepted_mutations = readRDS("results/accepted_mutations.rds") #passo in input a fit_model_selection_best_k? 

    S <- length(unique(accepted_mutations$segment_id))
    stanfit <- rstan::read_stan_csv(res$output_files())
  
    total_number_params <- K+(K*S)+2 # tau = K, w = K*S, phi, kappa (dirichlet reparametrization)
    k_inferred <- total_number_params
    N <- nrow(accepted_mutations)
    #L <- res_tibble %>% dplyr::filter(param == "lp__") %>% pull(median)
    log_lik<- extract_log_lik(stanfit)
    L <- mean(log_lik)
    BIC <- ((k_inferred * log(N)) - 2 * L) # %>% unname()
    AIC <- 2 * k_inferred - 2 * L

    #waic_result <- waic(log_lik_matrix)
    #waic_value <- waic_result$estimates[1, "Estimate"]  # WAIC estimate

    #LOO
    loo_result<-loo(log_lik)
    loo_value <- loo_result$estimates[3, "Estimate"]  # LOO-CV estimate

    #PSIS
    # log_ratios <- -1*(extract_log_lik(stanfit))
    # r_eff <- relative_eff(exp(-log_ratios))
    # psis_result <- psis(log_ratios, r_eff = r_eff)


    model_selection_tibble <- dplyr::bind_rows(model_selection_tibble, dplyr::tibble(K = K, BIC = BIC, AIC = AIC, LOO = loo_value, Log_lik = L))
    

    #instead of saving here save all together as soon as you can
    saveRDS(res, paste0("results/res",K,".rds"))
    saveRDS(input_data, paste0("results/input_data_",K,".rds"))
    
    p <- plotting(res,input_data, all_sim ,K, simulation_params)
    ggsave(paste0("./plots/plot_inference_",K,".png"), width = (12 + (simulation_params$number_events/2)), height = (16 + (simulation_params$number_events/2)), limitsize = FALSE, device = png, plot=p)
    
  }
  
   best_K <- model_selection_tibble %>% dplyr::filter(LOO == min(LOO)) %>% pull(K)
   input_data <- prepare_input_data(all_sim, karyo, best_K, purity)
   
    
   if (INIT==TRUE){
     inits_chain <- get_init(tau_single_inference, best_K)
     res <- fit_variational(input_data, max_attempts=max_attempts, initialization = inits_chain, INIT=TRUE)
   } else {
     res <- fit_variational(input_data, max_attempts=max_attempts, INIT=FALSE)
   }
   

  
   p_best_K <- plotting(res,input_data, all_sim, best_K, simulation_params)
   ggsave(paste0("./plots/plot_inference_",best_K,"_best_K.png"), width = (12 + (simulation_params$number_events/2)), height = (16 + (simulation_params$number_events/2)), limitsize = FALSE, device = png, plot=p_best_K)
   
  return(list(all_sim = all_sim, model_selection_tibble = model_selection_tibble, res_best_K=res, best_K=best_K, input_data=input_data, accepted_mutations=accepted_mutations
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

fit_variational <- function(input_data, max_attempts = 5, initialization = NULL, INIT = TRUE, initial_iter = 10000, grad_samples = 1, elbo_samples = 100) {
  model <- cmdstanr::cmdstan_model("../../CopyNumber/models/timing_mixed_simple.stan")
  attempt <- 0
  fit_successful <- FALSE
  res <- NULL
  iter <- initial_iter

  repeat {
    attempt <- attempt + 1
    message("Attempt ", attempt, " with ", iter, " iterations, ", grad_samples, " grad_samples, and ", elbo_samples, " elbo_samples")

    result <- tryCatch({
      if (INIT == TRUE) {
        res <- model$variational(data = input_data, init = list(initialization), iter = iter, grad_samples = grad_samples, elbo_samples = elbo_samples)
      } else {
        res <- model$variational(data = input_data, iter = iter, grad_samples = grad_samples, elbo_samples = elbo_samples)
      }

      # Verifica se l'output contiene il messaggio di errore specifico
      output_text <- capture.output(res$output())
      output_text <- paste(output_text, collapse = "\n")

      if (grepl("The number of dropped evaluations has reached its maximum amount", output_text)) {
        stop("Fit failed: Maximum number of discarded evaluations reached.")
      }

      # Estrazione dei draws
      draws <- res$draws(format = "matrix")

      # Controlla se ci sono NAs nei draws
      if (any(is.na(draws))) {
        stop("Fit failed: draws contain NAs.")
      }

      fit_successful <- TRUE
      res
    }, error = function(e) {
      message("An error occurred during inference: ", e$message)
      return(NULL) # Ritorna NULL in caso di errore
    })

    if (fit_successful) {
      message("Fit succeeded in attempt ", attempt)
      break
    }

    if (attempt >= max_attempts) {
      message("Maximum number of attempts reached. Stopping.")
      break
    }

    # Incrementa il numero di iterazioni, grad_samples e elbo_samples per il prossimo tentativo
    #iter <- iter * 2
    grad_samples <- grad_samples * 2
    elbo_samples <- elbo_samples * 2
  }

  # Controlla se l'inferenza è riuscita
  if (fit_successful) {
    message("Inference completed successfully.")
    return(res) # Ritorna il risultato se il fit è riuscito
  } else {
    message("The inference could not be completed successfully.")
    return(NULL) # Ritorna NULL se il fit non è riuscito
  }
}





#' prepare_input_data Function
#'
#' This function obtains the list of input data to be used for CN timing.
#' @param all_sim  simulated data
#' @param karyo karyotype posso togliere lo ricavo da all_sim XXXXXXX
#' @param K number of components
#' @param purity sample purity
#' @keywords input
#' @export
#' @examples
#' prepare_input_data()

prepare_input_data = function(all_sim, karyo, K, purity){
  
  all_sim <- all_sim[order(all_sim$segment_id), ]

  karyotype <- all_sim %>%
  group_by(segment_id) %>%
  summarise(karyotype = first(karyotype)) %>%
  pull(karyotype)

  seg <- all_sim %>%
  group_by(segment_id) %>%
  summarise(karyotype = first(karyotype)) %>%
  pull(segment_id)

  peaks <- matrix(0, nrow = length(karyotype), ncol = 2)
  for (i in 1:length(karyotype)) {
    peaks[i,] <- get_clonal_peaks(karyotype[i], purity)
  }

  

  cat("these are the peaks used to filter out mutations \n")
  print(peaks)
  cat("These are the Karyotypes extracted from the all_sim data, used in peaks\n")
  print(karyotype)
  print(seg)

  cat("These are the Karyotypes directly from the simulation, used for peaks in the model\n ")
  print(karyo)


  #peaks = peaks_inference(karyotype,purity)
  alpha = 0.05
  min_mutations_number = 2
  accepted_mutations <- data.frame()
  
  if (nrow(all_sim) > 0) {
    # Check if mutation is inside CI
    probs <- c(alpha/2, 1 - alpha/2)
    
    DP <- all_sim$DP
    NV <- all_sim$NV
    
    
    accepted_idx <- lapply(1:length(DP), function(i) {
      for (j in unique(all_sim$segment_id)) {
          if (all_sim$segment_id[i]==j){
          quantiles_1 <- qbinom(probs, DP[i], peaks[j,1])
          quantiles_2 <- qbinom(probs, DP[i], peaks[j,2])
        if (((NV[i] >= quantiles_1[1]) && (NV[i] <= quantiles_1[2])) | ((NV[i] >= quantiles_2[1]) && (NV[i] <= quantiles_2[2]))) {
          return(i)
        }
        }
      }
    }) %>% unlist()
    
    # Get only good mutations
    accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx], segment_id=all_sim$segment_name[accepted_idx], karyotype=all_sim$karyotype[accepted_idx], tau=all_sim$tau[accepted_idx] )
    
  }
  
  counts <- table(accepted_mutations$segment_id)
  minimum_number_per_segment <- all(counts >= min_mutations_number)
  
  if (minimum_number_per_segment == TRUE) {
  
  input_data <- list(
    S = length(unique(all_sim$segment_name[accepted_idx])),
    K = K,
    N = nrow(accepted_mutations),
    karyotype = lapply(karyo, karyo_to_int) %>% unlist(),
    seg_assignment = all_sim$segment_id[accepted_idx],
    peaks = peaks_inference(karyo,purity),
    NV = NV[accepted_idx],
    DP = DP[accepted_idx]
  )
  

  cat("These are the peaks which are passed in the inference task, more reliable let's say\n")
  print(input_data$peaks)
  print(karyotype)


  saveRDS(accepted_mutations, paste0("results/accepted_mutations.rds"))

  return(input_data)
} 
}






# INITIALIZATION 


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
  init_w <- as.matrix(res.fcm$u)
  epsilon <- 1e-5
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

  return(inits_chain)
}







#' fit_single_segments Function
#'
#' -NEED TO ADD FILTERING STEP - This function allows you to obtain simulated data and perform the single segments inference for CN timing, could be used to obtain initialization values for tau and w.
#' @param all_sim data
#' @param alpha confidece level
#' @param purity sample purity
#' @keywords fit
#' @export
#' @examples
#' fit_single_segments()

fit_single_segments = function(all_sim, alpha = .05, purity = 1){
  # FIRST OPTION: run single segment inference for initialization of tau and w
  model_single <- cmdstanr::cmdstan_model("../../CopyNumber/models/mixture_CNA_timing_binomial.stan")

  #prepare data for single segment inference : S = number of segments
  S <- length(unique(all_sim$segment_id))
  data_plot <- dplyr::tibble()
  inference_results <- dplyr::tibble()
  summarized_results <- dplyr::tibble()

  plots = c()
  tau = c()
  inference_results_tot = c()

  # fit the model for each segment
  for (segment_index in 1:S){

    print(segment_index)

    k = unique(all_sim$karyotype[all_sim$segment_id==segment_index])
    peaks_single <- get_clonal_peaks(k, purity)




    # filtering step







    input_data <- list(
      N = length(all_sim$segment_id[all_sim$segment_id==segment_index]),
      NV = all_sim$NV[all_sim$segment_id==segment_index],
      DP = all_sim$DP[all_sim$segment_id==segment_index],
      peaks = peaks_single
    )







    fit <- model_single$sample(data=input_data, iter_warmup=2000, iter_sampling=2000, chains=8, parallel_chains=8)


    # Compute tau posteriors
    tau_posteriors <- get_tau_posteriors(fit, k)$tau %>% unname() %>% as.numeric()

    q1 <- alpha / 2
    q2 <- 1 - alpha / 2
    tau_low <- quantile(tau_posteriors, q1) %>% unname()
    tau_high <- quantile(tau_posteriors, q2) %>% unname()
    tau_mean <- mean(tau_posteriors)

    inference_results <- dplyr::bind_rows(inference_results, dplyr::tibble(tau = tau_posteriors, segment = unique(all_sim$segment_id[all_sim$segment_id==i]), karyotype = k))
    summarized_results <- dplyr::bind_rows(summarized_results, dplyr::tibble(tau_low = tau_low, tau_mean = tau_mean, tau_high = tau_high, segment = unique(all_sim$segment_id[all_sim$segment_id==i]), karyotype = k))

    inference_results_tot = c(inference_results_tot, inference_results)
    data_plot <- dplyr::bind_rows(data_plot, dplyr::tibble(tau = tau_posteriors, segment = unique(all_sim$segment_id[all_sim$segment_id==i]), karyotype = k))

  }


  inference_single_segment <- data_plot %>%
    ggplot(mapping = aes(x=tau, fill=as.factor(segment))) +
    geom_histogram(alpha=.5, position = "identity", bins = 100)
  #plot_data <- plot_data + facet_wrap(vars(karyotype, segment))

  #posterior distribution of tau obtained from the single segment inference
  inference_single_segment
  ggsave("./plots/inference_single_segments.png", width = 12, height = 12,  device = png)

  tau_single_inference <- summarized_results$tau_mean
  return(tau_single_inference)
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
    peaks[i,] <- get_clonal_peaks(karyo[i], purity)
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

