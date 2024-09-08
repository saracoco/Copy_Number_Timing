library(patchwork)
library(loo)
library(bayesplot)
library(cmdstanr)
library(factoextra)
library(dplyr)
library(stringr) #for plotting add in the right script
library(fossil) #RI and ARI computation
library(gridExtra)



#setwd("C:/Users/sarac/CDS_git/Copy-Number-Timing/CopyNumber/")
#orfeo
setwd("/orfeo/cephfs/scratch/cdslab/scocomello/Copy_Number_Timing/CopyNumber")

original_dir <- getwd()

source("./CNTiming/R/simulate_functions.R")
source("./CNTiming/R/fitting_functions.R")
source("./CNTiming/R/plotting_functions.R")


self_name = "3"
new_dir = paste0("../",self_name) #relative path of the new created directory where to save the simulation results
dir.create(new_dir)

number_events = 10
number_clocks = 2

INIT = TRUE
epsilon = 0.20
n_simulations = 20
purity = 0.99

vector_karyo <- c("2:0", "2:1", "2:2")
weights_karyo <- c(0.33, 0.33, 0.33)


# get simulation parametes
coverage = 100 # average number of reads that align to a reference base
mu = 1e-4 # mutation rate
w = 1e-2 # cell division rate
l = 1e7 # length of the segment
time_interval = 20


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
  simulation_data_all_segments = get_simulation(data$taus, data$karyo, purity, coverage = 100) # the other parameters have default value assigned if none is specified
  simulation_data_all_segments <- simulation_data_all_segments[order(simulation_data_all_segments$segment_id), ]


  saveRDS(simulation_data_all_segments, "./results/all_sim_input_prepare_input_data.rds")



  
  Subtitle <- vector("list", (length(unique(simulation_data_all_segments$segment_id))+1))
  Subtitle[[1]]  <- paste0("Number of mutations per segment: ")
  num_mutations_all_segments <- c()

    for (i in seq_along(unique(simulation_data_all_segments$segment_id))) {
    segment <- unique(simulation_data_all_segments$segment_id)[i]
    num_mutations <- nrow(simulation_data_all_segments %>% filter(segment_id == segment))
    num_mutations_all_segments <- c(num_mutations_all_segments, num_mutations)
    Subtitle[[i+1]] <- paste0(segment, "=", num_mutations," ")
  }
  
  Subtitle <- paste(Subtitle, collapse = "   ")
  cat(Subtitle)

  mean_mut <- mean(num_mutations_all_segments)
  max_mut <- max(num_mutations_all_segments)
  min_mut <- min(num_mutations_all_segments)

  Subtitle_short <- paste0("Average number of mutations per segment: ", mean_mut, "  Minimum number of mutations per segment: ", min_mut, "  Maximum number of mutations per segment: ", max_mut )


  #add statistics on number of mutations from the simulation
  
  simulation_params <- list(
    vector_tau = vector_tau,
    vector_karyo = vector_karyo,
    weights_tau = weights_tau,
    weights_karyo = weights_karyo,
    taus = data$taus,
    karyo = data$karyo,
    purity = purity,
    number_events = number_events, # = nrow(vector-tau) / nrow(vector_karyo)
    number_clocks = number_clocks, # = unique(vector_tau)
    epsilon = epsilon
  )



  simulation_data_plot = simulation_data_all_segments %>% mutate (tau = round(tau, 2))
  plot_data <- simulation_data_plot %>% 
    ggplot(mapping = aes(x = NV / DP, fill = segment_name)) +
    geom_histogram(alpha = .5, position = "identity") +
    labs(
      title = "Distribution on the VAF for each segment in the simulated data",
      subtitle = paste0(Subtitle_short)
    )+
    facet_wrap(vars(karyotype, tau, segment_name), scales = "free_x", strip.position = "bottom") +
    theme_minimal() +
    theme(
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    plot.background = element_rect(fill = "white", color = NA),   # White plot background
    strip.background = element_rect(fill = "white", color = NA),  # White strip background
    strip.placement = "outside",   # Place facet labels outside
    axis.text.x = element_text(angle = 360, hjust = 1, color = "black", size = 8),  # Rotate and adjust x-axis text
    axis.ticks.x = element_line(color = "black"),  # Black x-axis ticks
    panel.spacing = unit(1, "lines"),  # Adjust space between facets
    strip.text.x = element_text(size = 10, color = "black"),  # Adjust and color strip text
    axis.line = element_line(color = "black"),  # Black axis lines
    axis.title.x = element_text(color = "black"),  # Black x-axis title
    axis.title.y = element_text(color = "black")   # Black y-axis title
  )


  #save plot of the simulated data in which we can see each single segment VAF distribution
  ggsave("./plots/simulation_data.png", plot = plot_data, width = 12 + simulation_params$number_events, height = 10 + simulation_params$number_events + (simulation_params$number_events/1.3), limitsize = FALSE,   device = png) 
  #simulation_params can be substituted in relation with simulation_data variables
  
  
  #in "fit model selection best K" the plots for each K inference is directly saved 
  results <- fit_model_selection_best_K(simulation_data_all_segments, data$karyo, purity, INIT = INIT, simulation_params = simulation_params)
  saveRDS(results, paste0("./results/results_simulation",i,".rds"))
  

  
  
  
  model_selection_plot = plotting_model_selection(results)
  model_selection_plot
  ggsave("./plots/model_selection_plot.png", plot = model_selection_plot, width = 12, height = 10,  device = png)
  
  model_selection <- results$model_selection_tibble
  saveRDS(model_selection, "./results/model_selection.rds")
  

  
  setwd(original_dir)
  
}

