library(dplyr)
library(ggplot2)
library(patchwork)
library(loo)
library(bayesplot)
library(cmdstanr)
library(factoextra)
library(ppclust)

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

set.seed(123)

tollerance = 0.01


source("./CNTiming/R/simulate_functions.R")
source("./CNtime/R/fitting_functions.R")
source("./CNtime/R/plotting_functions.R")



self_name = "UPN04_HLA"
new_dir = paste0("../",self_name) #relative path of the new created directory where to save the simulation results
dir.create(new_dir)


# LOAD DATA #
UPN04 <- readRDS("../../Data/extra_cnloh/alpha_beta/UPN04/mutations.rds")
UPN05 <- readRDS("../../Data/extra_cnloh/alpha_beta/UPN05/mutations.rds")

UPN_4_NV = UPN04 %>% filter(timing_classification %in% c("alpha private", "beta"),
                            PASS == TRUE)
UPN_5_NV = UPN05 %>% filter(timing_classification  %in% c("alpha private", "beta"), 
                            chr == "chr1",
                            PASS == TRUE)
UPN04_hla <- readRDS("../../Data/alpha_beta/UPN04/mutations.rds")
UPN05_hla <- readRDS("../../Data/alpha_beta/UPN05/mutations.rds")

UPN_4_hla_NV = UPN04_hla %>% ungroup %>% filter(timing_classification %in% c("alpha"), 
                                        PASS == TRUE)
UPN_5_hla_NV = UPN05_hla %>% filter(timing_classification %in% c("alpha", "beta"), 
                            PASS == TRUE)


data_hla <- list(UPN04 = UPN_4_NV, UPN05 = UPN_5_NV, UPN04_hla = UPN04_hla, UPN05_hla = UPN05_hla)
names <- c("UPN04","UPN04_hla")
name_segments <- c("segment 1", "segment 2")


options(bitmapType='cairo')



data <- dplyr::tibble()
karyo_all <- c()

for(i in 1:length(names)){

    data_single <- data_hla[[names[i]]]
    data_single <- data_single %>% rename(karyo = segment.REL, 
                        DP = DP.REL,
                        NV = NV.REL) %>%
                  mutate(segment_name = paste0(name_segments[i]),
                        segment_id = i,
                        karyo = as.character(karyo),
                        karyotype = karyo,
                        tau = karyo)
  karyo <- data_single$karyo[1]
  
  data <- dplyr::bind_rows(data, data_single)
  karyo_all <- c(karyo_all, karyo)
}

n_simulations = 1

for(i in 1:n_simulations){
      # Create a unique directory for each iteration
      iter_dir <- paste0("/simulation_iteration_", i)
      iter_dir <- paste0(new_dir,iter_dir)
      dir.create(iter_dir)
      setwd(iter_dir)
      dir.create(paste0("./plots"), showWarnings = TRUE)
      dir.create(paste0("./results"), showWarnings = FALSE)
      
      results <- fit_model_selection_best_K(data, karyo=karyo_all, max_attempts=4, purity=0.98, INIT=TRUE,  tollerance = tollerance)
      # results$model_selection_tibble

      saveRDS(results, paste0("./results/results_simulation",i,".rds"))

      model_selection_plot = plotting_model_selection(results)
      model_selection_plot
      ggsave("./plots/model_selection_plot.png", plot = model_selection_plot, width = 12, height = 10,  device = png)
      
      model_selection <- results$model_selection_tibble
      saveRDS(model_selection, "./results/model_selection.rds")
}



