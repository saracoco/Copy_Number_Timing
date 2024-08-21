library(dplyr)
library(ggplot2)
library(patchwork)
library(loo)
library(bayesplot)
library(cmdstanr)
library(factoextra)
library(ppclust)

setwd("C:/Users/sarac/CDS_git/Copy-Number-Timing/CopyNumber/")
source("./CNTiming/R/simulate_functions.R")
source("./CNTiming/R/fitting_functions.R")
source("./CNTiming/R/plotting_functions.R")


number_events=6
vector_tau<-c(0.01,0.5,0.9)
vector_karyo<-c("2:0", "2:1")
weignths_tau<-c(0.33,0.33,0.33) #vedere se il modello ha trovato le giuste proporzioni --> comp tau proportions
weights_karyo<-c(0.5,0.5)
purity = 0.99

data <- get_taus_karyo(number_events,vector_tau,vector_karyo,weignths_tau,weights_karyo)
all_sim = get_simulation(data$taus,data$karyo, purity=.99)
data_sim <- all_sim


# data_sim %>% 
#   ggplot(mapping=aes(x=NV/DP,col=karyotype))+
#   geom_histogram()+
#   facet_grid(~j +tau)

#change j with "segment_name" segment_id with segment_index
plot_data <- data_sim %>%
  ggplot(mapping = aes(x=NV/DP, fill=as.factor(j))) +
  geom_histogram(alpha=.5, position = "identity")
plot_data <- plot_data + facet_wrap(vars(karyotype, tau, j))
plot_data
ggsave("./plots/data.png", plot=plot_data, width = 12, height = 10, device = "png")







# directly
results <- fit_model_selection_best_K(data_sim, data$karyo, purity, INIT = FALSE)



results$model_selection_tibble
plot_best_K <- plotting(results$res_best_K,results$input_data,results$best_K)
ggsave("./plots/plot_inferenc_best_K.png", width = 12, height = 16, device = "png", plot = plot_best_K)







# fit with a specific number of K 
input_data_1 <- prepare_input_data(data_sim, data$karyo, K=1, purity=0.99)
res_1 <- fit_variational(input_data_1, max_attempts = 10, initialization = inits_chain1, INIT = FALSE, initial_iter = 10000)

p_1 <- plotting(res_1,input_data_1,1)



input_data_2 <- prepare_input_data(data_sim, data$karyo, K=2, purity=0.99)
res_2 <- fit_variational(input_data_2, max_attempts = 10, initialization = inits_chain1, INIT = FALSE, initial_iter = 10000)

p_2 <- plotting(res_2,input_data_2,2)









# FIRST OPTION INITIALIZATION: inizializzazione con cmeans sul rapporto alpha beta o solo alpha e beta come due variabili 
# model selection with initialization (c-MEANS) fit model selection INIT = TRUE ########################################## 

# indirectly
tau_single_inference <- fit_single_segments(data_sim, purity=0.99)
purity=0.99
K = 2
input_data <- prepare_input_data(all_sim, data$karyo, K, purity)
inits_chain1 <- get_init(tau_single_inference, K)
fit_result <- fit_variational(input_data, max_attempts = 10, initialization = inits_chain1, INIT = TRUE, initial_iter = 10000)

