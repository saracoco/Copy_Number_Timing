library(rstan)
library(tidyr)
library(ggthemes)
library(matrixStats)
library(dplyr)


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

fit_model_selection_best_K = function(all_sim, karyo, purity=0.99, max_attempts=4, INIT=TRUE, simulation_params = simulation_params){

  set.seed(123)
  if (INIT==TRUE){
    #tau_single_inference <- fit_single_segments(all_sim, purity=purity) # uncomment to run initialization with single segment model
  }

  #MODEL SELECTION
  model_selection_tibble <- dplyr::tibble()
  #Regola empirica: In generale, si può considerare come massimo plausibile un numero di cluster pari a 
  # sqrt(n/2), dove n è il numero di punti dati.
  
  if (length(karyo) <= 15){
    k_max = (length(karyo)/2)-1
  } else { k_max = sqrt(length(karyo))  #does not converge with initialization (don't know why (ex #seg = 10 --> k_max = 5, with 5 does not converge))
  }
  
  for (K in 1:k_max) {
    input_data <- prepare_input_data(all_sim, karyo, K, purity)
    accepted_mutations = readRDS("results/accepted_mutations.rds") #potrei prenderli direttamente da input_data quando non farò model validation

    if (INIT==TRUE){
      # inits_chain <- get_init(tau_single_inference, K)     # uncomment to run initialization with single segment model
      inits_chain <- get_init_simpler(accepted_mutations, K)
      # print(paste0("These are the values used for initializing the model ", inits_chain)) # uncomment to check for initialization problems
      res <- fit_variational(input_data, max_attempts=max_attempts, initialization = inits_chain, INIT = TRUE)
    } else {
      res <- fit_variational(input_data, max_attempts=max_attempts, INIT = FALSE)
    }
        # Plot ELBO values over iterations
        output_files <- res$latent_dynamics_files()
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

          ############################################
          # SELECT THE BEST MODEL 
          stanfit <- rstan::read_stan_csv(res$output_files())
          S <- length(unique(accepted_mutations$segment_id))
          total_number_params <- K+((K-1)*S)+2                    # tau = K, w = K*S, phi, kappa (dirichlet reparametrization)
          N <- nrow(accepted_mutations)
          
          # mcmc_hist(fit_vb$draws("theta"))
          # fit_vb$summary()

          # Check log likelihood values 
          lp__ <- res$draws("lp__")
          cat (paste0(" mode of lp__ = unnormalized_log_lik = ", mode(lp__), "last value = ", lp__[100], " dimension of lp__ = ", dim(lp__),"\n"))
          lp_approx__ <- res$draws("lp_approx__")
          cat (paste0(" mode of lp_approx after rowSums = ", mode(rowSums(lp_approx__)), " dimension of lp_approx__ = ", dim(lp_approx__),"\n"))
          log_lik_matrix <- extract_log_lik(stanfit, parameter_name = "log_lik", merge_chains = TRUE) # merge chains potrebbe non servire
          # cat (paste0(" mode of log_lik after rosSums = likelihood from generated quantities", log_lik_matrix[1], "dimension of log_lik = ", dim(log_lik_matrix), "\n"))
          cat (paste0(" mode of log_lik after rosSums = likelihood from generated quantities", mode(rowSums(log_lik_matrix)),"last value = ", rowSums(log_lik_matrix)[100], "dimension of log_lik = ", dim(log_lik_matrix), "\n"))
          # lp_approx <- res$lp_approx() #res_tibble %>% dplyr::filter(param == "lp__") %>% pull(median)
          # lp <- res$lp() #res_tibble %>% dplyr::filter(param == "lp__") %>% pull(median)
          # lp__ <- as.matrix(stanfit, pars = "lp__")
          # lp__ <- res$lp()

          log_lik_total_per_sample <- rowSums(log_lik_matrix)
          L <- mode(log_lik_total_per_sample)

          BIC <- ((total_number_params * log(N)) - 2 * L) # %>% unname()
          AIC <- 2 * total_number_params - 2 * L

          # Generate names for w[sim_params_num, K] format
          names_weights <- outer(1:simulation_params$number_events, 1:K, 
                                FUN = function(i, j) paste0("w[", i, ",", j, "]"))
          # Convert to a vector for easier extraction
          draws = 1000
          names_weights <- as.vector(names_weights)
          w_ICL <- as.matrix(stanfit, pars = names_weights)
          dim(w_ICL) = c(draws,S*K)
          w_ICL <- apply(w_ICL, 2, mode) # mode over draws
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
          for (k in post ){
            post_k = post[k]
            log_post_k = log(post_k + 0.000001)
            post_k_entr = post_k * log_post_k * mut_per_seg
            post_k_entr = sum(post_k_entr)
            post_k_entr = -1 * (post_k_entr)
            res_entropy = res_entropy + post_k_entr
          }
          entropy = res_entropy
    
          ICL = BIC + entropy
          print(paste0("entropy",entropy))


          #LOO
          loo_result<-loo(log_lik_matrix)
          loo_value <- loo_result$estimates[3, "Estimate"]  # LOO-CV estimate

          #PSIS
          # log_ratios <- -1*(extract_log_lik(stanfit))
          # r_eff <- relative_eff(exp(-log_ratios))
          # psis_result <- psis(log_ratios, r_eff = r_eff)

          #WAIC
          #waic_result <- waic(log_lik_matrix)
          #waic_value <- waic_result$estimates[1, "Estimate"]# WAIC 


          model_selection_tibble <- dplyr::bind_rows(model_selection_tibble, dplyr::tibble(K = K, BIC = BIC, AIC = AIC, LOO = loo_value, Log_lik = L, ICL = ICL))
          

          # # prior predictive check without using fit_variational
          #   model_stan <- cmdstanr::cmdstan_model("../../CopyNumber/models/timing_mixed_simple.stan")

          #   fit_prior <- model_stan$sample(
          #       data = input_data,         # Data for prior predictive checks now the data of the simulation, is it wrong?
          #       iter_sampling = 1000,     # Number of prior samples
          #       chains = 4,               # Number of chains
          #       fixed_param = TRUE        # Sample from priors only
          #     )

        
              # # AGGIUNGI PHI, K, THETA
              # # PRIOR PLOT 
              # draws_df <- res$draws(format = "df")  # Extract draws as a data frame
              # print(draws_df)
              #       # Extract all columns related to tau and w
              # tau_columns <- grep("^tau_prior\\[", colnames(draws_df), value = TRUE)
              # w_columns <- grep("^w_prior\\[", colnames(draws_df), value = TRUE)

              # # Extract the draws for tau and w
              # tau_draws <- draws_df[, tau_columns]
              # w_draws <- draws_df[, w_columns]


              # tau_long <- pivot_longer(as.data.frame(tau_draws), cols = everything(), names_to = "tau_param", values_to = "tau_value")


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

                  
              ggsave(paste0("./plots/priors_tau_",K,".png"),  width = (8 + (simulation_params$number_events/2)), height = ( 8 + (simulation_params$number_events/2)), limitsize = FALSE, device = png, plot=areas_tau)


                  # names_w <- paste("w_prior[", 1:K, "]", sep = "")
                  intervals_weigths_per_tau <- list()
                    for (k in 1:K){
                          names_weights <- paste("w_prior[",1:simulation_params$number_events,",", k, "]", sep = "") 
                          intervals_weigths_per_tau[[k]] <- mcmc_areas_ridges(draws_matrix, pars = names_weights, point_est = "median", prob = 0.8)+
                                      labs(
                                        title =  str_wrap( paste0("Prior distributions of the weigths for tau ",k), width = 30 + K + sqrt(simulation_params$number_events)),
                                        subtitle = "with median and 80% and 95% intervals"
                                      ) +
                                      theme(plot.background = element_rect(fill = "white"),text = element_text(size = 15))

                    }
                    areas_w <- gridExtra::grid.arrange(grobs = intervals_weigths_per_tau, ncol=K) #add global title

        
              ggsave(paste0("./plots/priors_w_",K,".png"),  width = (15 + (simulation_params$number_events/2)), height = (8 + (simulation_params$number_events/2)), limitsize = FALSE, device = png, plot=areas_w)


              # # Plot tau
              # prior_tau <- ggplot(tau_long, aes(x = tau_value, fill = tau_param)) +
              #   geom_histogram(binwidth = 0.05, alpha = 0.7, position = "identity") +
              #   labs(title = "Prior Distributions for tau", x = "tau", y = "Frequency") +
              #   theme_minimal()



              # # Reshape w to long format
              # w_long <- pivot_longer(as.data.frame(w_draws), cols = everything(), names_to = "w_param", values_to = "w_value")

              # # Separate w_param into 'segment' and 'clock'
              # w_long <- separate(w_long, w_param, into = c("param", "segment", "clock"), sep = "\\[|,|\\]", convert = TRUE)

              # # Plot w[s,k] by segment and clock
              # prior_w <- ggplot(w_long, aes(x = w_value, fill = interaction(segment, clock))) +
              #   geom_histogram(binwidth = 0.05, alpha = 0.7, position = "identity") +
              #   labs(title = "Prior Distributions for w[s,k]", x = "w[s,k]", y = "Frequency") +
              #   theme_minimal() +
              #   facet_wrap(~ interaction(segment, clock), scales = "free")



        

          # prior_plot <- (prior_tau|prior_w )  +
          #   plot_layout(widths = c(8), heights = c(10)) +
          #   plot_annotation(
          #     title = paste0("Prior Predictive Check"),
          #     subtitle = paste0("Simulation with ", simulation_params$number_clocks," clocks, ", simulation_params$number_events, " segments, epsilon = ", simulation_params$epsilon, " purity = ", simulation_params$purity)
          #   ) & theme(text = element_text(size = 12), plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12), axis.text = element_text(size = 12 ), plot.caption = element_text(size = 10 ))
          
          # ggsave(paste0("./plots/priors_",K,".png"), width = 8, height = 5, limitsize = FALSE, device = png, plot=prior_plot)








          #instead of saving here save all together as soon as you can
          saveRDS(res, paste0("results/res",K,".rds"))
          saveRDS(input_data, paste0("results/input_data_",K,".rds"))


          p <- plotting(res,input_data, all_sim ,K, simulation_params)
          ggsave(paste0("./plots/plot_inference_",K,".png"), width = (12 + (simulation_params$number_events/2)), height = (16 + (simulation_params$number_events/2)), limitsize = FALSE, device = png, plot=p)
          

          plot_partition = plotting_cluster_partition(res, K, VALIDATION=TRUE)
          ggsave(paste0("./plots/inferred_partition_",K,".png"), width = 8, height = 5, limitsize = FALSE, device = png, plot=plot_partition)

  }
  
   best_K <- model_selection_tibble %>% dplyr::filter(BIC == min(BIC)) %>% pull(K)
   input_data <- prepare_input_data(all_sim, karyo, best_K, purity)
   
    
   if (INIT==TRUE){
    accepted_mutations = readRDS("results/accepted_mutations.rds")
    inits_chain <- get_init_simpler(accepted_mutations, best_K)
    #  inits_chain <- get_init(tau_single_inference, best_K)
     print(paste0("These are the values used for initializing the model ",inits_chain))
     res <- fit_variational(input_data, max_attempts=max_attempts, initialization = inits_chain, INIT=TRUE)
   } else {
     res <- fit_variational(input_data, max_attempts=max_attempts, INIT=FALSE)
   }
   

  
   p_best_K <- plotting(res,input_data, all_sim, best_K, simulation_params)
   ggsave(paste0("./plots/plot_inference_",best_K,"_best_K.png"), width = (12 + (simulation_params$number_events/2)), height = (16 + (simulation_params$number_events/2)), limitsize = FALSE, device = png, plot=p_best_K)
   
   p_elbo_iter <- plotting_elbo(k_max)
   ggsave(paste0("./elbo_vs_iterations_.png"), plot = p_elbo_iter)


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

#fit variational taking the best run 
fit_variational <- function(input_data, max_attempts = 5, initialization = NULL, INIT = TRUE, initial_iter = 10000, grad_samples = 10, elbo_samples = 100) {
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
            tol_rel_obj = 0.0001
          )
          print(res$init())

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
            tol_rel_obj = 0.0001
          )
        }


        # print(res$metadata())
        # print(res$cmdstan_summary()


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
        if (length(output_files) > 0) {
          invalid_log_output <- readLines(output_files[1])
          invalid_log_lines <- grep("Invalid", invalid_log_output, value = TRUE)
          if (length(invalid_log_lines) > 0) {
            message("Invalid log evaluations found:\n", paste(invalid_log_lines, collapse = "\n"))
          }
        }
        
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
    accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx], segment_id=all_sim$segment_name[accepted_idx], karyotype=all_sim$karyotype[accepted_idx], tau=all_sim$tau[accepted_idx] , segment_index=all_sim$segment_id[accepted_idx])
    
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
    peaks = peaks, #peaks_inference(karyo,purity),
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
get_init_simpler = function(accepted_mutations, K, phi=c(), kappa=5){

    accepted_mutations <- accepted_mutations[order(accepted_mutations$segment_index), ]

      karyotype <- accepted_mutations %>%
      group_by(segment_index) %>%
      summarise(karyotype = first(karyotype)) %>%
      pull(karyotype)

      seg <- accepted_mutations %>%
      group_by(segment_index) %>%
      summarise(karyotype = first(karyotype)) %>%
      pull(segment_index)

      peaks <- matrix(0, nrow = length(karyotype), ncol = 2)
      for (i in 1:length(karyotype)) {
        peaks[i,] <- get_clonal_peaks(karyotype[i], purity)
      }


      alpha = 0.05
      min_mutations_number = 2
      
      if (nrow(accepted_mutations) > 0) {
        # Check if mutation is inside CI
        probs <- c(alpha/2, 1 - alpha/2)
        
        DP <- accepted_mutations$DP
        NV <- accepted_mutations$NV
        
        
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
      proportions$tau_posterior[i] <- 3 * proportions$proportion_beta[i] / 
                                      (2 * proportions$proportion_beta[i] + proportions$proportion_alpha[i])
    } else {
      proportions$tau_posterior[i] <- 2 * proportions$proportion_beta[i] / 
                                      (2 * proportions$proportion_beta[i] + proportions$proportion_alpha[i])
    }
  }

  # Check if tau_posterior contains any NA values
  if (any(is.na(proportions$tau_posterior))) {
    stop("tau_posterior contains NA values, check the loop calculations.")
  }

  # Apply fuzzy c-means clustering
  res.fcm <- fcm(as.matrix(proportions$tau_posterior), centers = K)

  init_taus <- c(res.fcm$v)
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

    # filtering step

    all_sim_single = all_sim %>% filter(segment_id==segment_index)
    k = unique(all_sim$karyotype[all_sim$segment_id==segment_index])
    peaks_single <- get_clonal_peaks(k, purity)

    probs <- c(alpha/2, 1 - alpha/2)

    DP <- all_sim_single$DP
    NV <- all_sim_single$NV

    accepted_idx <- lapply(1:length(all_sim_single$DP), function(i) {
      for (p in peaks_single) {
        
        quantiles <- qbinom(probs, all_sim_single$DP[i], p)

        if ((all_sim_single$NV[i] >= quantiles[1]) && (all_sim_single$NV[i] <= quantiles[2])) {
          return(i)
        }
      }
    }) %>% unlist()

    # Get only good mutations
    accepted_mutations <- data.frame(DP = all_sim_single$DP[accepted_idx], NV = all_sim_single$NV[accepted_idx])



    input_data <- list(
      N = length(all_sim_single$NV[accepted_idx]),
      NV = all_sim_single$NV[accepted_idx], #all_sim$NV[all_sim$segment_id==segment_index],
      DP = all_sim_single$DP[accepted_idx], #all_sim$DP[all_sim$segment_id==segment_index],
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

    inference_results <- dplyr::bind_rows(inference_results, dplyr::tibble(tau = tau_posteriors, segment = unique(all_sim$segment_id[all_sim$segment_id==segment_index]), karyotype = k))
    summarized_results <- dplyr::bind_rows(summarized_results, dplyr::tibble(tau_low = tau_low, tau_mean = tau_mean, tau_high = tau_high, segment = unique(all_sim$segment_id[all_sim$segment_id==segment_index]), karyotype = k))

    inference_results_tot = c(inference_results_tot, inference_results)
    data_plot <- dplyr::bind_rows(data_plot, dplyr::tibble(tau = tau_posteriors, segment = unique(all_sim$segment_id[all_sim$segment_id==segment_index]), karyotype = k))

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




mode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}


