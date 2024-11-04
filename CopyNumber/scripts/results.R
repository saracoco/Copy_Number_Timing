# setwd("scripts")


library(dplyr)
library(ggplot2)
if(!require(gtools)) install.packages("gtools")
library(gtools)


library(dplyr)
library(tidyr)

# Assuming you have these lists set up
simulation_folders <- list("../../01",  "../../02", "../../03", "../../7", "../../8", "../../9", "../../1",  "../../2", "../../3", "../../4", "../../5", "../../6")  # Add your simulation folder paths
simulation_sizes <- c(10, 20, 15, 10, 20, 15, 10, 20, 15, 10, 20, 15)  # Example of dataset sizes per simulation
true_clusters <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)  # Example of true K values per simulation

# Initialize an empty data frame to store results
final_results <- data.frame()

# Iterate over each simulation folder and setting
for (i in seq_along(simulation_folders)) {
  sim_dir <- simulation_folders[[i]]
  sim_size <- simulation_sizes[i]
  true_k <- true_clusters[i]
  
  # Get all iteration folders inside the current simulation directory
  iteration_folders <- list.dirs(sim_dir, recursive = FALSE, full.names = TRUE)
  iteration_folders <- iteration_folders[grepl("simulation_iteration_\\d+$", iteration_folders)]
  iteration_folders <- mixedsort(iteration_folders)

  
  # Loop through each simulation iteration
  for (iter_folder in iteration_folders[-length(iteration_folders)]) {
    
    # Define the path to the results folder inside the iteration folder
    results_dir <- file.path(iter_folder, "results")
    
    # Read in the RI and MAE for all K and the best K
    K_max <- length(list.files(results_dir, pattern = "RI_")) - 1  # Assuming RI files follow a pattern
    RI_best <- readRDS(file.path(results_dir, "RI_best_K.rds"))
    MAE_best <- readRDS(file.path(results_dir, "MAE_best_K.rds"))
    
    # Loop over all K values for RI and MAE
    for (K in 1:K_max) {
      RI_K <- readRDS(file.path(results_dir, paste0("RI_", K, ".rds")))
      MAE_K <- readRDS(file.path(results_dir, paste0("MAE_", K, ".rds")))
      
      # Add results to final_results data frame
      final_results <- final_results %>%
        bind_rows(data.frame(
          SimulationType = sim_size,  # or other factor you want to include
          TrueK = true_k,
          K = K,
          RI = RI_K,
          MAE = MAE_K,
          IsBestK = (K == K_max),  # Logical to mark best K
          Metric = "RI"
        )) %>%
        bind_rows(data.frame(
          SimulationType = sim_size,
          TrueK = true_k,
          K = K,
          RI = RI_best,  # Best K value
          MAE = MAE_best,
          IsBestK = TRUE,
          Metric = "MAE"
        ))
    }
  }
}

# Check structure of final_results
str(final_results)

# At this point, final_results contains all the data in long format with columns:
# SimulationType, TrueK, K, RI, MAE, IsBestK, and Metric


# Basic boxplot for RI or MAE with faceting
p <- ggplot(final_results, aes(x = as.factor(K), y = MAE)) +
  geom_boxplot() +
  facet_grid(SimulationType ~ TrueK) +  # You can replace TrueK or SimulationType with other factors
  labs(x = "K", y = "MAE", title = "RI Distribution by Simulation Settings") +
  ylim(0, 1) +
  theme_bw()

ggsave("MAE_all_sim.pdf")














library(dplyr)
library(ggplot2)
if(!require(gtools)) install.packages("gtools")
library(gtools)


library(dplyr)
library(tidyr)


simulation_dirs_list <- list("../../01",  "../../02", "../../03", "../../7", "../../8", "../../9", "../../1",  "../../2", "../../3", "../../4", "../../5", "../../6")  # Add your simulation folder paths
simulation_sizes <- c(10, 20, 15, 10, 20, 15, 10, 20, 15, 10, 20, 15)   # Example of dataset sizes per simulation
true_K_list <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)  # Example of true K values per simulation

library(tidyverse)



# Function to get best K for each model selection criterion
get_best_K <- function(model_selection_table) {
  best_K_BIC <- which.min(model_selection_table$BIC)
  best_K_AIC <- which.min(model_selection_table$AIC)
  best_K_LOO <- which.min(model_selection_table$LOO)
  best_K_ICL <- which.min(model_selection_table$ICL)
  
  return(c(BIC = best_K_BIC, AIC = best_K_AIC, LOO = best_K_LOO, ICL = best_K_ICL))
}


# Initialize an empty data frame for results
best_K_results <- data.frame(Simulation_Setting = character(),
                             Iteration = integer(),
                             TrueK = integer(),
                             BIC = integer(),
                             AIC = integer(),
                             LOO = integer(),
                             ICL = integer())

# Assuming you have a list of true Ks for each simulation type
# and that simulation_dirs is a list of directories for each simulation setting
for (sim_setting in 1:length(simulation_dirs_list)) {
  sim_dir <- simulation_dirs_list[[sim_setting]]
  true_K <- true_K_list[sim_setting]  # corresponding true K for each simulation setting
  
  iteration_folders <- list.dirs(sim_dir, recursive = FALSE, full.names = TRUE)
  iteration_folders <- iteration_folders[grepl("simulation_iteration_\\d+$", iteration_folders)]
  iteration_folders <- mixedsort(iteration_folders)


  for (iteration in 1:(length(iteration_folders)-1)) {  # Assuming up to 20 iterations (adjust if necessary)
    iter_dir <- file.path(sim_dir, paste0("simulation_iteration_", iteration), "results")
    
    if (dir.exists(iter_dir)) {
      # Load the model_selection table for this iteration
      model_selection_table <- readRDS(file.path(iter_dir, "model_selection.rds"))
      
      # Get the best K for each criterion
      best_Ks <- get_best_K(model_selection_table)
      
      # Add results to the dataframe
      best_K_results <- rbind(best_K_results, 
                              data.frame(Simulation_Setting = sim_setting, 
                                         Iteration = iteration,
                                         TrueK = true_K,
                                         BIC = best_Ks['BIC'], 
                                         AIC = best_Ks['AIC'], 
                                         LOO = best_Ks['LOO'], 
                                         ICL = best_Ks['ICL']))
    }
  }
}





best_K_results <- best_K_results %>%
  mutate(
    Match_BIC = ifelse(TrueK == BIC, 1, 0),
    Match_AIC = ifelse(TrueK == AIC, 1, 0),
    Match_LOO = ifelse(TrueK == LOO, 1, 0),
    Match_ICL = ifelse(TrueK == ICL, 1, 0)
  )



  # Summary of correct predictions per criterion for each simulation setting
match_summary_by_setting <- best_K_results %>%
  group_by(Simulation_Setting) %>%
  summarize(
    BIC_correct = mean(Match_BIC) * 100,
    AIC_correct = mean(Match_AIC) * 100,
    LOO_correct = mean(Match_LOO) * 100,
    ICL_correct = mean(Match_ICL) * 100
  )




print(match_summary_by_setting)


long_data <- match_summary_by_setting %>%
  pivot_longer(cols = c(BIC_correct, AIC_correct, LOO_correct, ICL_correct),
               names_to = "Criterion", values_to = "Correct_Prediction")

match_summary_long <- match_summary_by_setting %>%
  pivot_longer(cols = ends_with("_correct"), names_to = "Criterion", values_to = "Correct_Match") %>%
  mutate(Criterion = str_replace(Criterion, "_correct", ""))

# Plot the results
ggplot(match_summary_long, aes(x = Simulation_Setting, y = Correct_Match, fill = Criterion)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0, 100) +
  labs(title = "Percentage of Correct K Predictions by Model Selection Criteria", 
       x = "Simulation Setting", 
       y = "Percentage of Correct Predictions (%)") +
  theme_minimal() +
  facet_wrap(~ Criterion)


ggsave("compare_MS_scores.pdf")                               





























  # Initialize an empty data frame for MAE and RI results per criterion
mae_ri_results <- data.frame(Simulation_Setting = character(),
                             Iteration = integer(),
                             Criterion = character(),
                             K = integer(),
                             MAE = numeric(),
                             RI = numeric())

for (sim_setting in 1:length(simulation_dirs_list)) {
  sim_dir <- simulation_dirs_list[[sim_setting]]
  
  iteration_folders <- list.dirs(sim_dir, recursive = FALSE, full.names = TRUE)
  iteration_folders <- iteration_folders[grepl("simulation_iteration_\\d+$", iteration_folders)]
  iteration_folders <- mixedsort(iteration_folders)


  for (iteration in 1:(length(iteration_folders)-1)) {# Assuming up to 20 iterations (adjust if necessary)
    iter_dir <- file.path(sim_dir, paste0("simulation_iteration_", iteration), "results")
    
    if (dir.exists(iter_dir)) {
      # Load the model_selection table for this iteration
      model_selection_table <- readRDS(file.path(iter_dir, "model_selection.rds"))
      
      # Get the best K for each criterion
      best_Ks <- get_best_K(model_selection_table)
      
      # For each criterion, load the corresponding MAE and RI for the best K
      for (criterion in names(best_Ks)) {
        best_k <- best_Ks[criterion]
        
        # Load the corresponding MAE and RI for the best K
        mae <- readRDS(file.path(iter_dir, paste0("MAE_", best_k, ".rds")))
        ri  <- readRDS(file.path(iter_dir, paste0("RI_", best_k, ".rds")))
        
        # Add results to the dataframe
        mae_ri_results <- rbind(mae_ri_results, 
                                data.frame(Simulation_Setting = sim_setting,
                                           Iteration = iteration,
                                           Criterion = criterion,
                                           K = best_k,
                                           MAE = mae,
                                           RI = ri))
      }
    }
  }
}

# View the compiled MAE and RI results
print(mae_ri_results)

  # Reshape the data for combined plotting
mae_ri_results <- mae_ri_results %>%
  pivot_longer(cols = c("MAE", "RI"), names_to = "Metric", values_to = "Value")

# SEPARATE PLOTS
# Split the data by simulation settings into two subsets
mae_ri_results_1to4 <- mae_ri_results %>% filter(Simulation_Setting %in% 1:4)
mae_ri_results_5to8 <- mae_ri_results %>% filter(Simulation_Setting %in% 5:8)
mae_ri_results_9to12 <- mae_ri_results %>% filter(Simulation_Setting %in% 9:12)


# Plot for Simulation Settings 1-4
p1 <- ggplot(mae_ri_results_1to4, aes(x = Criterion, y = Value, fill = Criterion)) +
  geom_boxplot() +
  facet_grid(Metric ~ Simulation_Setting, scales = "free_y") +
  labs(title = "MAE and RI Distributions (Simulations 1-4) for Best K Selected by Each Criterion",
       x = "Model Selection Criterion", y = "Value") +
  theme_minimal()

ggsave("compare_criteria_performances_1to4.pdf", plot = p1, width = 10, height = 8)

# Plot for Simulation Settings 5-8
p2 <- ggplot(mae_ri_results_5to8, aes(x = Criterion, y = Value, fill = Criterion)) +
  geom_boxplot() +
  facet_grid(Metric ~ Simulation_Setting, scales = "free_y") +
  labs(title = "MAE and RI Distributions (Simulations 5-8) for Best K Selected by Each Criterion",
       x = "Model Selection Criterion", y = "Value") +
  theme_minimal()

ggsave("compare_criteria_performances_5to8.pdf", plot = p2, width = 10, height = 8)

# Plot for Simulation Settings 5-8
p3 <- ggplot(mae_ri_results_9to12, aes(x = Criterion, y = Value, fill = Criterion)) +
  geom_boxplot() +
  facet_grid(Metric ~ Simulation_Setting, scales = "free_y") +
  labs(title = "MAE and RI Distributions (Simulations 9-12) for Best K Selected by Each Criterion",
       x = "Model Selection Criterion", y = "Value") +
  theme_minimal()

ggsave("compare_criteria_performances_9to12.pdf", plot = p3, width = 10, height = 8)























# Plot the MAE distribution
ggplot(mae_ri_results, aes(x = Criterion, y = MAE, fill = Criterion)) +
  geom_boxplot() +
  facet_wrap(~ Simulation_Setting, scales = "free_y") +
  labs(title = "MAE Distribution for Best K Selected by Each Criterion",
       x = "Model Selection Criterion", y = "MAE") +
  theme_minimal()



  # Plot the RI distribution
ggplot(mae_ri_results, aes(x = Criterion, y = RI, fill = Criterion)) +
  geom_boxplot() +
  facet_wrap(~ Simulation_Setting, scales = "free_y") +
  labs(title = "RI Distribution for Best K Selected by Each Criterion",
       x = "Model Selection Criterion", y = "RI") +
  theme_minimal()


  # Reshape the data for combined plotting
mae_ri_long <- mae_ri_results %>%
  pivot_longer(cols = c("MAE", "RI"), names_to = "Metric", values_to = "Value")

# Plot combined MAE and RI
p <- ggplot(mae_ri_long, aes(x = Criterion, y = Value, fill = Criterion)) +
  geom_boxplot() +
  facet_grid(Metric ~ Simulation_Setting, scales = "free_y") +
  labs(title = "MAE and RI Distributions for Best K Selected by Each Criterion",
       x = "Model Selection Criterion", y = "Value") +
  theme_minimal()

ggsave("compare_criteria_performances.pdf", width=20, height=10)                               
























