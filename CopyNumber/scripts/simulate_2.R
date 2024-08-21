library(patchwork)
library(loo)
library(bayesplot)
library(cmdstanr)
library(factoextra)
library(dplyr)


#setwd("C:/Users/sarac/CDS_git/Copy-Number-Timing/CopyNumber/")
#orfeo
setwd("/orfeo/cephfs/scratch/cdslab/scocomello/Copy-Number-Timing/CopyNumber")

original_dir <- getwd()

source("./CNTiming/R/simulate_functions.R")
source("./CNTiming/R/fitting_functions.R")
source("./CNTiming/R/plotting_functions.R")


S = 10
K = 3
INIT = FALSE
epsilon = 0.15
n_simulations = 20
purity = 0.99

vector_karyo <- c("2:0", "2:1", "2:2")
weights_karyo <- c(0.33, 0.33, 0.33)



for(i in 1:n_simulations){
  # Create a unique directory for each iteration
  iter_dir <- paste0("./simulation_2_iteration_", i)
  dir.create(iter_dir)
  setwd(iter_dir)
  dir.create(paste0("./plots"), showWarnings = TRUE)
  dir.create(paste0("./results"), showWarnings = FALSE)
  
  
  number_events = S
  number_clusters = K
  vector_tau = c()
  
  for (j in 1:K){
    vector_tau[j] = runif(1, 0)
    if (j != 1){
      vector_tau[j] = runif(1, 0)
      while (abs(vector_tau[j-1] - vector_tau[j]) < epsilon){
        vector_tau[j] = runif(1, 0)
      }
    }
  }
  weights_tau <- rep(1/K, K)
  
  data <- get_taus_karyo(number_events, vector_tau, vector_karyo, weights_tau, weights_karyo)
  all_sim = get_simulation(data$taus, data$karyo, purity)
  data_sim <- all_sim
  #add statistics on number of mutations from the simulation
  
  plot_data <- data_sim %>%
    ggplot(mapping = aes(x = NV / DP, fill = as.factor(j))) +
    geom_histogram(alpha = .5, position = "identity") +
    facet_wrap(vars(karyotype, tau, j))
  
  ggsave("./plots/original_data.png", plot = plot_data, width = 12, height = 10, device = "png")
  
  simulation_params <- list(
    vector_tau = vector_tau,
    vector_karyo = vector_karyo,
    weights_tau = weights_tau,
    weights_karyo = weights_karyo,
    taus = data$taus,
    karyo = data$karyo,
    purity = purity
  )
  
  
  #in "fit model selection best K" the plots for each K inference is directly saved 
  results <- fit_model_selection_best_K(data_sim, data$karyo, purity, INIT = INIT, simulation_params = simulation_params)
  
  
  
  
  
  
  
  model_selection <- results$model_selection_tibble
  best_K <- results$best_K
  
  
  bic_plot <- ggplot(data = model_selection, aes(x = K, y = BIC)) + 
    geom_line(aes(colour = "BIC Line")) + 
    geom_point(aes(colour = "BIC Point")) +
    geom_point(aes(x = best_K, y = BIC[K == best_K], colour = "Best K Point")) +
    geom_point(aes(x = K[BIC == min(BIC)], y = min(BIC), colour = "Minimum Point")) +
    scale_colour_manual(name = "Legend",
                        values = c("BIC Line" = "blue", 
                                   "BIC Point" = "blue", 
                                   "Best K Point" = "red",
                                   "Minimum Point" = "green"),)
  
  aic_plot <- ggplot(data = model_selection, aes(x = K, y = AIC)) + 
    geom_line(aes(colour = "AIC Line")) + 
    geom_point(aes(colour = "AIC Point")) +
    geom_point(aes(x = best_K, y = AIC[K == best_K], colour = "Best K Point")) +
    geom_point(aes(x = K[AIC == min(AIC)], y = min(AIC), colour = "Minimum Point")) +
    scale_colour_manual(name = "Legend",
                        values = c("AIC Line" = "blue", 
                                   "AIC Point" = "blue", 
                                   "Best K Point" = "red",
                                   "Minimum Point" = "green"))
  
  
  loo_plot <- ggplot(data = model_selection, aes(x = K, y = LOO)) + 
    geom_line(aes(colour = "LOO Line")) + 
    geom_point(aes(colour = "LOO Point")) +
    geom_point(aes(x = best_K, y = LOO[K == best_K], colour = "Best K Point")) +
    geom_point(aes(x = K[LOO == min(LOO)], y = min(LOO), colour = "Minimum Point")) +
    scale_colour_manual(name = "Legend",
                        values = c("LOO Line" = "blue", 
                                   "LOO Point" = "blue", 
                                   "Best K Point" = "red",
                                   "Minimum Point" = "green"))
  
  
  log_lik_plot <- ggplot(data = model_selection, aes(x = K, y = Log_lik)) + 
    geom_line(aes(colour = "Log likelihood Line")) + 
    geom_point(aes(colour = "Log likelihood Point")) +
    geom_point(aes(x = best_K, y = Log_lik[K == best_K], colour = "Best K Point")) +
    geom_point(aes(x = K[Log_lik == max(Log_lik)], y = max(Log_lik), colour = "Minimum Point")) +
    scale_colour_manual(name = "Legend",
                        values = c("Log likelihood Line" = "blue", 
                                   "Log likelihood Point" = "blue", 
                                   "Best K Point" = "red",
                                   "Minimum Point" = "green"))
  
  
  
  Subtitle <- vector("list", S)
  
  for (i in 1:S) {
    n_parameters <- i+(i*S)+2
    Subtitle[[i]] <- paste0(" ", i, ": ", n_parameters," parameters")
  }
  
  Subtitle <- paste(Subtitle, collapse = "\n")
  
  
  
  
  model_selection_plot <- (bic_plot | aic_plot) / (loo_plot|log_lik_plot) +
    plot_annotation(
      title = "Model selection graphs: score vs number of clusters" ,
      subtitle = paste0 ("Correspondence between number of clusters and number of parameters:  \n", Subtitle),
      caption = "caption"
    ) & theme(text = element_text(size = 8), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), axis.text = element_text(size = 8), plot.caption = element_text(size = 5))
  
  
  
  
  ggsave("./plots/model_selection_plot.png", plot = model_selection_plot, width = 12, height = 10, device = "png")
  saveRDS(model_selection, "model_selection.rds")
  
  
  
  setwd(original_dir)
  
}

