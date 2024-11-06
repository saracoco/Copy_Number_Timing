library(tidyverse)
source("scripts/utils.R")
source("scripts/fit.R")
source("scripts/plot.R")
source("scripts/beta_mixture.R")
source("scripts/COMPARE_TIMING_PLOT.R")


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



data_all <-readRDS("data/sara_dataset.rds")

data <- data_all[[3]]

#purity 1 0.46
#purity 2 0.504
#purity 3 0.611




purity <- data$metadata$purity
print(paste0("purity = ", purity))
possible_k = c("2:1", "2:2", "2:0")
segments=data$cna
mutations = data$mutations

data_all <- dplyr::tibble()
segments <- segments %>%
    drop_na(Major, minor)
n_segments <- nrow(segments)
j=0


for (segment_idx in 1:n_segments) { #n_segments   i=5 106 16
    if (segment_idx != 16){  # chance wrt the specific PCAWG segment
    print(segment_idx)

    # Segments
    segment <- segments[segment_idx, ]
    chr <- segment$chr
    print(chr)
    print(segment)

    segment_id <- paste(chr, segment$from, segment$to, sep = "_")

    # Get karyotype
    Major <- segment$Major
    minor <- segment$minor

    k <- paste(Major, minor, sep=':')

    peaks <- get_clonal_peaks(k, purity)

    if (k %in% possible_k) {
         print("YES")
        # Get info for mutations
        segment_mutations <- mutations %>%
        filter(chr == segment$chr,from > segment$from, to < segment$to) %>%
        drop_na(DP)
            
        DP <- segment_mutations$DP
        NV <- segment_mutations$NV



                    
        if (nrow(segment_mutations)>2){


              j = j+1
            # Get only good mutations
            
            
            
            cli::cli_alert_info("Adding segment with index {.val {segment_idx}}")
        
            data_single <- data.frame(DP = DP, NV = NV, segment_id = j, tau = 0 )

            # # model and input data
            # input_data <- list(
            #   N = nrow(accepted_mutations),
            #   NV = accepted_mutations$NV,
            #   DP = accepted_mutations$DP,
            #   peaks = peaks,
            #   beta_dispersion = beta_binomial_disp
            # )

            data_all <- dplyr::bind_rows(data_all, dplyr::tibble(DP = DP, NV = NV, karyotype = k, segment_name = paste0("segment ",j), segment_id = j, tau = 0, segment_name_real = segment_idx))

        }
        
    
    }

    }
}






set.seed(133)
tollerance = 0.001
print(paste0("tolerance: ", tollerance))    

max_attempts = 2

#setwd("C:/Users/sarac/CDS_git/Copy-Number-Timing/CopyNumber/")
#orfeo

sim_list = c("WHICH_SAMPLE_3")

i = 1
setwd("/orfeo/cephfs/scratch/cdslab/scocomello/Copy_Number_Timing/CopyNumber")

original_dir <- getwd()

source("./CNTiming/R/simulate_functions.R")
source("./CNtime/R/fitting_functions.R")
source("./CNtime/R/plotting_functions.R")

# simulate the 3 events happening at 3 different time points
self_name = sim_list
new_dir = paste0("../",self_name) #relative path of the new created directory where to save the simulation results
dir.create(new_dir)

INIT = TRUE
options(bitmapType='cairo')

n_simulations = 1

# Create a unique directory for each iteration
iter_dir <- paste0("/simulation_iteration_", i)
iter_dir <- paste0(new_dir,iter_dir)
dir.create(iter_dir)
setwd(iter_dir)
dir.create(paste0("./plots"), showWarnings = TRUE)
dir.create(paste0("./results"), showWarnings = FALSE)

    
      
# vector_tau
taus = data_all$tau

# in karyo there should be a list of karyotypes that is as long as the number of segments
# karyo = data_all$karyotype
data <- data_all[order(data_all$segment_id), ]
karyo <- data %>%
group_by(segment_id) %>%
summarise(karyotype = first(karyotype)) %>%
pull(karyotype)


simulation_data_all_segments <- data_all[order(data_all$segment_id), ]
saveRDS(data_all, "./results/all_sim_input_prepare_input_data.rds")

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



#in fit model selection best K the plots for each K inference is directly saved 
results <- fit_model_selection_best_K(simulation_data_all_segments, karyo, purity, INIT = INIT, max_attempts = max_attempts, tollerance = tollerance )
saveRDS(results, paste0("./results/results_simulation",i,".rds"))


model_selection_plot = plotting_model_selection(results)
model_selection_plot
ggsave("./plots/model_selection_plot.png", plot = model_selection_plot, width = 12, height = 10,  device = png)

model_selection <- results$model_selection_tibble
saveRDS(model_selection, "./results/model_selection.rds")

