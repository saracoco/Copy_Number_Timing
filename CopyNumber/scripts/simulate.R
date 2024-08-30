library(patchwork)
library(loo)
library(bayesplot)
library(cmdstanr)
library(factoextra)
library(dplyr)
library(stringr) #for plotting add in the right script
library(fossil) #RI and ARI computation<


#setwd("C:/Users/sarac/CDS_git/Copy-Number-Timing/CopyNumber/")
#orfeo
setwd("/orfeo/cephfs/scratch/cdslab/scocomello/Copy_Number_Timing/CopyNumber")

original_dir <- getwd()

source("./CNTiming/R/simulate_functions.R")
source("./CNTiming/R/fitting_functions.R")
source("./CNTiming/R/plotting_functions.R")


self_name = "try"
new_dir = paste0("../",self_name) #relative path of the new created directory where to save the simulation results
dir.create(new_dir)

number_events = 6
number_clocks = 2

INIT = FALSE
epsilon = 0.20
n_simulations = 1
purity = 0.99

vector_karyo <- c("2:0", "2:1", "2:2")
weights_karyo <- c(0.33, 0.33, 0.33)

options(bitmapType='cairo')




for(i in 1:n_simulations){
  # Create a unique directory for each iteration
  iter_dir <- paste0("/simulation_iteration_", i)
  iter_dir <- paste0(new_dir,iter_dir)
  dir.create(iter_dir)
  setwd(iter_dir)
  dir.create(paste0("./plots"), showWarnings = TRUE)
  dir.create(paste0("./results"), showWarnings = FALSE)
  
  

  vector_tau = rep(0, number_clocks)
  
  for (j in 1:number_clocks){
    vector_tau[j] = runif(1, 0)
    if (j != 1){
      while (!all ( abs(vector_tau[1:j-1] - vector_tau[j]) > epsilon  )   ){
        vector_tau[j] = runif(1, 0)
      }
    }
  }
  weights_tau <- rep(1/number_clocks, number_clocks)
  
  data <- get_taus_karyo(number_events, vector_tau, vector_karyo, weights_tau, weights_karyo)
  all_sim = get_simulation(data$taus, data$karyo, purity)
  data_sim <- all_sim
  #add statistics on number of mutations from the simulation
  
  simulation_params <- list(
    vector_tau = vector_tau,
    vector_karyo = vector_karyo,
    weights_tau = weights_tau,
    weights_karyo = weights_karyo,
    taus = data$taus,
    karyo = data$karyo,
    purity = purity,
    number_events = number_events,
    number_clocks = number_clocks,
    epsilon = epsilon
  )

  data_sim_plot = data_sim %>% mutate (tau = round(tau, 2))
  plot_data <- data_sim_plot %>% 
    ggplot(mapping = aes(x = NV / DP, fill = as.factor(j))) +
    geom_histogram(alpha = .5, position = "identity") +
    labs(
      title = "Histogram of the VAF spectrum, per segment"
    )+
    facet_wrap(vars(karyotype, tau, j))
  
  ggsave("./plots/original_data.png", plot = plot_data, width = 12 +  simulation_params$number_events, height = 10 + simulation_params$number_events + (simulation_params$number_events/1.3), limitsize = FALSE,   device = png)
  #simulation_params can be substituted in relation with data_sim variables
  
  
  #in "fit model selection best K" the plots for each K inference is directly saved 
  results <- fit_model_selection_best_K(data_sim, data$karyo, purity, INIT = INIT, simulation_params = simulation_params)
  saveRDS(results, paste0("./results/results",K,"_",all_sim$j[1],".rds"))
  

  
  
  
  model_selection_plot = plotting_model_selection(results)
  model_selection_plot
  ggsave("./plots/model_selection_plot.png", plot = model_selection_plot, width = 12, height = 10,  device = png)
  
  model_selection <- results$model_selection_tibble
  saveRDS(model_selection, "./results/model_selection.rds")
  

  
  setwd(original_dir)
  
}

