library(dplyr)
library(ggplot2)
library(patchwork)
library(loo)
library(bayesplot)
library(cmdstanr)
library(factoextra)
library(ppclust)


setwd("C:/Users/sarac/CDS/CopyNumber")
source("./CNTiming/R/simulate_functions.R")
source("./CNTiming/R/fitting_functions.R")
source("./CNTiming/R/plotting_functions.R")


# LOAD DATA #
setwd("E:/scratch/CDS_ORFEO/Timing_CDS")
UPN04_extra <- readRDS("Data/extra_cnloh/alpha_beta/UPN04/mutations.rds")
UPN05_extra <- readRDS("Data/extra_cnloh/alpha_beta/UPN05/mutations.rds")
UPN04_alpha_beta <- readRDS("Data/alpha_beta/UPN04/mutations.rds")
UPN05_alpha_beta <- readRDS("Data/alpha_beta/UPN05/mutations.rds")
# FILTERING #
UPN04_extra_NV = UPN04_extra %>% filter(timing_classification %in% c("alpha private", "beta"),PASS == TRUE)
UPN05_extra_NV = UPN05_extra %>% filter(timing_classification  %in% c("alpha private", "beta"),chr == "chr1",PASS == TRUE)
UPN04_alpha_beta_NV = UPN04_alpha_beta %>% ungroup %>% filter(timing_classification %in% c("alpha"),PASS == TRUE)
UPN05_alpha_beta_NV = UPN05_alpha_beta %>% filter(timing_classification %in% c("alpha", "beta"),PASS == TRUE)


data_lsh <- list(UPN04 = UPN04_extra_NV, UPN05 = UPN05_extra_NV, UPN04_LSH = UPN04_alpha_beta_NV, UPN05_LSH = UPN05_alpha_beta_NV)
#names <- c("UPN04","UPN04_LSH", "UPN05", "UPN05_LSH")
names <- c("UPN04","UPN04_LSH")
#names <- c("UPN05", "UPN05_LSH")


setwd("E:/scratch/CDS_ORFEO/Timing_CDS/initialization")


data <- dplyr::tibble()
karyo_all <- c()

for(i in 1:length(names)){
  
  data_single <- data_lsh[[names[i]]]

  data_single <- data_single %>% rename(karyo = segment.REL, 
                        DP = DP.REL,
                        NV = NV.REL) %>%
                  mutate(j = paste0(names[i]),
                        segment_id = i,
                        karyo = as.character(karyo),
                        karyotype = karyo)
  
  karyo <- data_single$karyo[1]
  
  data <- dplyr::bind_rows(data, data_single)
  karyo_all <- c(karyo_all, karyo)
}



##### INFERENCE MODEL SEELECTION MULTIPLE K ##################################################################################
results <- fit_model_selection_best_K(data, karyo=karyo_all, purity=0.98, INIT=FALSE)

results$model_selection_tibble

p <- plotting(results$res_best_K,results$input_data,results$best_K)
ggsave(paste0("./plots/plot_best_K_",names[1],".png"), width = 12, height = 16, device = "png", plot=p)




##### SINGLE MODELS RESULTS #######################################################################################à

#UPN04
input_data_1 <- readRDS("./results/input_data1_UPN04.rds")
input_data_2 <- readRDS("./results/input_data2_UPN04.rds")

res_1 <- readRDS("./results/res1_UPN04.rds")
res_2 <- readRDS("./results/res2_UPN04.rds")


#UPN05
input_data_1 <- readRDS("./results/input_data1_UPN05.rds")
input_data_2 <- readRDS("./results/input_data2_UPN05.rds")

res_1 <- readRDS("./results/res1_UPN05.rds")
res_2 <- readRDS("./results/res2_UPN05.rds")

#### Plot accepted mutations #########################################################################################


accepted_mutations_UPN04 <- readRDS("./results/input_data1_UPN04.rds")
accepted_mutations_UPN05 <- readRDS("./results/input_data1_UPN05.rds")


#UPN04
accepted_mutations_UPN04 <- as.data.frame(accepted_mutations_UPN04)
df <- filter(accepted_mutations_UPN04, seg_assignment == 1)
hist_data_extra <- df %>% 
  ggplot(mapping = aes(x=(NV/DP))) +
  geom_histogram(binwidth=0.01, , fill="blue", alpha = 0.5, position = "identity") +
  theme(legend.position = 'top') +
  xlim(0, 1) +
  ylim(0, 5) +
  ggtitle("UPN04 extra", data_lsh["UPN04_LSH"]$UPN04_LSH$segment.REL[1])


df <- filter(accepted_mutations_UPN04, seg_assignment == 2)
hist_data_LSH <- df %>% 
  ggplot(mapping = aes(x=(NV/DP), fill=timing_classification)) +
  geom_histogram(binwidth=0.01, , fill="blue", alpha = 0.5, position = "identity") +
  theme(legend.position = 'top') +
  xlim(0, 1) +
  ylim(0, 5) +
  ggtitle("UPN04 LSH", data_lsh["UPN04"]$UPN04$segment.REL[1])


plot_mutations_04 <- (hist_data_extra + hist_data_LSH) + 
  plot_layout(widths = c(6), heights = c(10)) +
  plot_annotation(
    title = 'Accepted mutations after inference on UPN04 LSH and UPN04 extra CN event ',
    subtitle = " ", #Extra event is on chr
    caption = "" #caption
  ) & theme(text = element_text(size = 8), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), axis.text = element_text(size = 8), plot.caption = element_text(size = 5))
plot_mutations_04

ggsave(paste0("./plots/accepted_mutations_UPN04.png"), width = 14, height = 8, device = "png", plot=plot_mutations_04)






#UPN05
accepted_mutations_UPN05 <- as.data.frame(accepted_mutations_UPN05)
df <- filter(accepted_mutations_UPN04, seg_assignment == 1)
hist_data_extra <- df %>% 
  ggplot(mapping = aes(x=(NV/DP))) +
  geom_histogram(binwidth=0.01, , fill="blue", alpha = 0.5, position = "identity") +
  theme(legend.position = 'top') +
  xlim(0, 1.5) +
  ylim(0, 5) +
  ggtitle("UPN05 extra", data_lsh["UPN05_LSH"]$UPN05_LSH$segment.REL[1])


df <- filter(accepted_mutations_UPN05, seg_assignment == 2)
hist_data_LSH <- df %>% 
  ggplot(mapping = aes(x=(NV/DP), fill=timing_classification)) +
  geom_histogram(binwidth=0.01, , fill="blue", alpha = 0.5, position = "identity") +
  theme(legend.position = 'top') +
  xlim(0, 1.5) +
  ylim(0, 5) +
  ggtitle("UPN05 LSH", data_lsh["UPN05"]$UPN05$segment.REL[1])


plot_mutations_05 <- (hist_data_extra + hist_data_LSH) + 
  plot_layout(widths = c(6), heights = c(10)) +
  plot_annotation(
    title = 'Accepted mutations after inference on UPN05 LSH and UPN05 extra CN event ',
    subtitle = " ", #Extra event is on chr
    caption = "" #caption
  ) & theme(text = element_text(size = 8), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), axis.text = element_text(size = 8), plot.caption = element_text(size = 5))
plot_mutations_05

ggsave(paste0("./plots/accepted_mutations_UPN05.png"), width = 14, height = 8, device = "png", plot=plot_mutations_05)









#### INFERENCE FOR A SINGLE k ##################################################################################

input_data_1 <- prepare_input_data(data, karyo_all, K=1, purity=0.98)
res_1 <- fit_variational(input_data_1, max_attempts = 10, initialization = inits_chain1, INIT = FALSE, initial_iter = 10000)

p_1 <- plotting(res_1,input_data_1,1)
ggsave(paste0("./plots/plot_inference",names[1],"_1.png"), width = 12, height = 16, device = "png", plot=p_1)



input_data_2 <- prepare_input_data(data, karyo_all, K=2, purity=0.98)
res_2 <- fit_variational(input_data_2, max_attempts = 10, initialization = inits_chain1, INIT = FALSE, initial_iter = 10000)

p_2 <- plotting(res_2,input_data_2,2)
ggsave(paste0("./plots/plot_inference",names[1],"_2.png"), width = 12, height = 16, device = "png", plot=p_2)


#########################################################################################################à
