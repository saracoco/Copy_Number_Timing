library(patchwork) # plotting
library(loo)        # model selection
library(bayesplot)  # plotting
library(cmdstanr)   # inference
library(factoextra)  
library(dplyr)     # simulation + inference
library(stringr) # for plotting add in the right script
library(fossil) # RI and ARI computation
library(gridExtra) # plotting
library(ppclust) # initialization step 
library(tidyr)

set.seed(133)
tolerance = 0.01
print(paste0("tolerance: ", tolerance))
max_attempts = 2

sim_list = c(0)
number_clocks_list = c(3)
number_events_list = c(8 )
epsilon_list = c(0.20)

setwd("../")
original_dir <- getwd()

source("./CNTiming/R/simulate_functions.R")
source("./CNtime/R/fitting_functions.R")
source("./CNtime/R/plotting_functions.R")


i = 1
self_name = as.character(sim_list[i])
new_dir = paste0("../",self_name) #relative path of the new created directory where to save the simulation results
dir.create(new_dir)

number_events = number_events_list[i]
number_clocks = number_clocks_list[i]
print(paste0("number of clocks: ", number_clocks) )

INIT = TRUE
epsilon = epsilon_list[i]
n_simulations = 1
purity = 0.98

vector_karyo <- c("2:0", "2:1", "2:2")
weights_karyo <- c(0.33, 0.33, 0.33)

# get simulation parametes
coverage = 100 # average number of reads that align to a reference base
mu = 1e-4 # mutation rate
w = 1e-2 # cell division rate
l = 1e7 # length of the segment
time_interval = 7

options(bitmapType='cairo')

i=1
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
data_simulation <- get_taus_karyo(number_events, vector_tau, vector_karyo, weights_tau, weights_karyo)
simulation_data_all_segments = get_simulation(data_simulation$taus, data_simulation$karyo, purity, time_interval, l, mu, w, coverage) # the other parameters have default value assigned if none is specified
data <- simulation_data_all_segments[order(simulation_data_all_segments$segment_id), ]

simulation_params <- list(
  vector_tau = vector_tau,
  vector_karyo = vector_karyo,
  weights_tau = weights_tau,
  weights_karyo = weights_karyo,
  taus = data_simulation$taus,
  karyo = data_simulation$karyo,
  purity = purity,
  number_events = number_events, # = nrow(vector-tau) / nrow(vector_karyo)
  number_clocks = number_clocks, # = unique(vector_tau)
  epsilon = epsilon
)

input_for_fit = list(data = data, purity = purity, max_attempts = max_attempts, INIT = INIT, tolerance = tolerance)
saveRDS(input_for_fit, paste0("./results/input_for_fit.rds"))
# input_for_fit = readRDS("../../0/simulation_iteration_1/results/input_for_fit.rds")
# results <- fit_model_selection_best_K(data, input_for_fit$purity, input_for_fit$max_attempts, input_for_fit$INIT, input_for_fit$tolerance )
results <- fit_model_selection_best_K(data, purity, max_attempts, INIT, tolerance)

# now results is a list of lenght k_max of results form the variational method with the summary statistics of the draws from the approximate posterior

input <- readRDS("./results/input_.rds")
accepted_mutations <- input$accepted_mutations
results_model_selection <- model_selection(data, results, accepted_mutations, compute_external_metric = FALSE)

best_K <- results_model_selection$best_K
model_selection_tibble <- results_model_selection$model_selection_tibble

# input_list <- readRDS("./results/input_list.rds")
plotting_fit_all_k(results, accepted_mutations, data)


entropy <- results_model_selection$entropy_list

# model_selection_plot <- plotting_model_selection(model_selection_tibble, best_K, entropy_per_segment_matrix_norm, entropy_per_segment_matrix, accepted_mutations)
model_selection_plot <- plotting_model_selection(model_selection_tibble, best_K, accepted_mutations, entropy)
ggsave("./plots/model_selection_plot.png", plot = model_selection_plot, width = 12, height = 10,  device = png)


      

# setwd(original_dir)