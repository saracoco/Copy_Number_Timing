library(patchwork)
library(loo)
library(bayesplot)
library(cmdstanr)
library(factoextra)
library(dplyr)
library(stringr) #for plotting add in the right script
library(fossil) #RI and ARI computation
library(gridExtra)
library(ppclust)
library(tidyr)


set.seed(133)
tollerance = 0.0001
print(paste0("tolerance: ", tollerance))    

max_attempts = 2

#setwd("C:/Users/sarac/CDS_git/Copy-Number-Timing/CopyNumber/")
#orfeo

sim_list = c("6")
number_clocks_list = c(4)
print(paste0("number_clocks_list: ", number_clocks_list))
number_events_list = c(16)
epsilon_list = c(0.20)

for (i in (1:length(sim_list))) {

    # setwd("D:/scratch/Copy_Number_Timing/CopyNumber")
    setwd("/orfeo/cephfs/scratch/cdslab/scocomello/Copy_Number_Timing/CopyNumber")

    original_dir <- getwd()

    source("./CNTiming/R/simulate_functions.R")
    source("./CNtime/R/fitting_functions.R")
    source("./CNtime/R/plotting_functions.R")

    # simulate the 3 events happening at 3 different time points
    self_name = as.character(sim_list[i])
    new_dir = paste0("../",self_name) #relative path of the new created directory where to save the simulation results
    dir.create(new_dir)

    number_events = number_events_list[i]
    print(paste0("number of clocks: ", number_clocks ))
    number_clocks = number_clocks_list[i]
    print(paste0("number_events: ", number_events ))
    

    INIT = TRUE
    epsilon = epsilon_list[i]
    n_simulations = 20
    purity = 0.98

    # get simulation parametes
    coverage = 100 # average number of reads that align to a reference base
    mu = 1e-4 # mutation rate
    w = 1e-2 # cell division rate
    l = 1e7 # length of the segment
    time_interval = 5

    options(bitmapType='cairo')


    for(i in 1:n_simulations){
      print(paste0("n_simulations: ", n_simulations))

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
  
      vector_karyo <- c("2:0", "2:1", "2:2")
      weights_karyo <- c(0.33, 0.33, 0.33)
      data <- get_taus_karyo(number_events, vector_tau, vector_karyo, weights_tau, weights_karyo)

      taus = data$taus
      print(paste0("taus ", taus))

      karyo = data$karyo
      print(paste0("karyo ", karyo))

      simulation_data_all_segments = get_simulation(taus, karyo, purity, coverage = 100) # the other parameters have default value assigned if none is specified
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
        Subtitle_short <- paste0("Mutations per segment: average =", mean_mut, ",  min = ", min_mut, ", max = ", max_mut )


      # plot simulated data settings
      simulation_plot <- plot_simulate(simulation_data_all_segments)
      ggsave("./plots/simulation_plot.png", plot = simulation_plot, width = 10, height = 7,  device = png)



        #add statistics on number of mutations from the simulation
        simulation_params <- list(
          vector_tau = taus,
          vector_karyo = karyo,
          taus = taus,
          karyo = karyo,
          purity = purity,
          number_events = number_events, # = nrow(vector-tau) / nrow(vector_karyo)
          number_clocks = number_clocks # = unique(vector_tau)
        )
      
        #in fit model selection best K the plots for each K inference is directly saved 
        results <- fit_model_selection_best_K(simulation_data_all_segments, karyo, purity, INIT = INIT, max_attempts = max_attempts, tollerance = tollerance, compute_external_metric = TRUE )
        saveRDS(results, paste0("./results/results_simulation",i,".rds"))
        
        
        model_selection_plot = plotting_model_selection(results)
        model_selection_plot
        ggsave("./plots/model_selection_plot.png", plot = model_selection_plot, width = 12, height = 10,  device = png)
        
        model_selection <- results$model_selection_tibble
        saveRDS(model_selection, "./results/model_selection.rds")
        
        setwd(original_dir)

        print(paste0("directory at the end of simulation 1: ", original_dir))
      

    }


}



