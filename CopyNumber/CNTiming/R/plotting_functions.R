
#' plotting Function
#'
#' This function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param res result
#' @param input_data data
#' @param K number of mixture components
#' @keywords input
#' @export
#' @examples
#' plotting()


plotting <- function(res, input_data, all_sim, K, simulation_params){


  ################################################## fit results ######################################################################
  draws <- res$draws(format = "matrix")

  color_scheme_set("blue")
  #manage to get the K from the model fit directly rather than as input
  names <- paste("tau[", 1:K, "]", sep = "")
  areas_tau <- mcmc_areas(
    draws,
    pars = names,
    prob = 0.8, # 80% intervals
    prob_outer = 0.95, # 99%
    point_est = "median"
  )+
    labs(
      title = "Approximate Posterior distributions",
      subtitle = "With median and 80% and 95% intervals"
    )+
    xlim(0, 1) # + scale_x_continuous(breaks = c(1:5), labels = c("A", "B", "C", "D", "E"))

  color_scheme_set("blue")
  intervals_weigths_per_tau <- list()
  for (k in 1:K){
        names_weights <- paste("w[",1:input_data$S,",", k, "]", sep = "") 
        intervals_weigths_per_tau[[k]] <- mcmc_intervals(draws, pars = names_weights, point_est = "median", prob = 0.8, prob_outer = 0.95)+
                    labs(
                      title =  str_wrap( paste0("Posterior distributions of the weigths for tau ",k), width = 30 + K + sqrt(input_data$S)),
                      subtitle = "With median and 80% and 95% intervals"
                    )
  }
  intervals_weigths_per_tau <- gridExtra::grid.arrange(grobs = intervals_weigths_per_tau, ncol=K) #add global title


  #################################################### posterior predictive checkn#############################################
  
  stanfit <- rstan::read_stan_csv(res$output_files())
  y_rep <- as.matrix(stanfit, pars = "NV_pred")
  y = input_data$NV
  #distribution of replicated data vs real data
  ppc <- ppc_dens_overlay(
    y = input_data$NV,
    yrep = y_rep) +
    labs(
      title = "Posterior Predictive checks with 10000 iterations",
      subtitle = "Density distribution"
    )
  #empirical cumulative distribution
  ecdf_compare <-ppc_ecdf_overlay(y,y_rep)
  #predictive intervals vs observed values
  intervals_compare <- ppc_intervals(y,y_rep) +
    labs(
      subtitle = " Observed data vs Predicted values for each segment",
    )
  #Compare statistics
  mean_compare <- ppc_stat(y,y_rep, stat = 'mean')
  max_compare <- ppc_stat(y,y_rep, stat = 'max')
  min_compare <- ppc_stat(y,y_rep, stat = 'min')
  median_compare <- ppc_stat(y,y_rep, stat = 'median')

  #######################################################################################################

  #############################à prepare plot for final plot ##############################################
  accepted_mutations <- readRDS("results/accepted_mutations.rds") #Give as input do not call it from inside the function for the package
  
  Subtitle <- vector("list", (length(unique(accepted_mutations$segment_id))+1))
  Subtitle[[1]]  <- paste0("Number of mutations per segment: ")
  num_mutations_all <- c()
    for (i in seq_along(unique(accepted_mutations$segment_id))) {
    segment <- unique(accepted_mutations$segment_id)[i]
    num_mutations_single <- nrow(accepted_mutations %>% filter(segment_id == segment))
    num_mutations_all <- c(num_mutations_all, num_mutations_single)
    Subtitle[[i+1]] <- paste0(segment, "=", num_mutations_single," ")
  }
  Subtitle <- paste(Subtitle, collapse = "   ")
  cat(Subtitle)

  mean_mut <- mean(num_mutations_all)
  max_mut <- max(num_mutations_all)
  min_mut <- min(num_mutations_all)
  Subtitle_short <- paste0("Mutations per segment: average =", mean_mut, ",  min = ", min_mut, ", max = ", max_mut )
  
  accepted_mutations = accepted_mutations %>% mutate (tau = round(tau, 2))
  plot_filtered_data <- accepted_mutations %>%
  ggplot(mapping = aes(x = NV / DP, fill = segment_id)) +
  geom_histogram(alpha = .5, position = "identity") +
  labs(x = "VAF")+
  labs(
    title = paste0( str_wrap("Histogram of the VAF spectrum, per segment, resulting from the simulation (only the data used in the inference after the filtering step are plotted here)", width = 90 + K + (simulation_params$number_events) ) ),
    subtitle = paste0( str_wrap(Subtitle_short, width = 90 + K + (simulation_params$number_events) ) )
  )+
  facet_wrap(vars(karyotype, tau))


  ################################# plot data with accepted mutations modified with peaks #######################################à
  vline_positions <- data.frame(
  segment_id = paste0("segment ",1:input_data$S),
  peaks_1 = input_data$peaks[,1],
  peaks_2 = input_data$peaks[,2]
  )
  accepted_mutations <- merge(accepted_mutations, vline_positions, by = "segment_id")

  #plot filter data to save separately and compare vaf distribution before and after the filtering step
  plot_filtered_data_complete <- accepted_mutations %>%
  ggplot(mapping = aes(x = NV / DP, fill = segment_id)) +
  geom_histogram(alpha = 0.6, position = "identity", color = "black") + # Add a black border to bars
  facet_wrap(vars(karyotype, tau, segment_id), scales = "free_x", strip.position = "bottom") +
  geom_vline(aes(xintercept = peaks_1), color = "grey", linetype = "dashed") +  # Dashed lines for peaks
  geom_vline(aes(xintercept = peaks_2), color = "grey", linetype = "dashed") + 
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),   
    strip.background = element_rect(fill = "white", color = NA),  
    strip.placement = "outside",   
    axis.text.x = element_text(angle = 360, hjust = 1, color = "black", size = 8),  
    axis.ticks.x = element_line(color = "black"),  
    panel.grid.major = element_line(color = "grey90"), # Subtle grid lines
    panel.spacing = unit(1, "lines"),  
    strip.text.x = element_text(size = 10, color = "black"),  
    axis.line = element_line(color = "black"),  
    axis.title.x = element_text(color = "black", size = 12),  # Bigger font for axis titles
    axis.title.y = element_text(color = "black", size = 12),   
    plot.title = element_text(size = 16, face = "bold"),  # Bold title for better readability
    plot.subtitle = element_text(size = 12, face = "italic")  # Italic subtitle
  ) +
  labs(
    x = "VAF",
    y = "Count",
    title = "Histogram of the VAF spectrum, per segment",
    subtitle = paste0("Only the data used in the inference after the filtering step are plotted \n", Subtitle_short)
  )+
    xlim(0, 1)

  ggsave("./plots/plot_filtered_data_complete.png", plot = plot_filtered_data_complete, width = 8 + simulation_params$number_events, height = 6 + simulation_params$number_events + (simulation_params$number_events/1.3), limitsize = FALSE,   device = png) 


######################################### number of mutations  per segment taking the simulation data before filtering ###############################################
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

  vline_positions <- data.frame(
  segment_name = paste0("segment ",1:input_data$S),
  peaks_1 = input_data$peaks[,1],
  peaks_2 = input_data$peaks[,2]
  )
  all_sim_ <- merge(all_sim, vline_positions, by = "segment_name")

#### PLOT WITH ALL MUTATIONS BEFORE FILTERING, IF NEEDED ##############################################################s
  # simulation_data_plot = all_sim_ %>% mutate (tau = round(tau, 2))
  # plot_data <- simulation_data_plot %>% 
  #   ggplot(mapping = aes(x = NV / DP, fill = segment_name)) +
  #   geom_histogram(alpha = .5, position = "identity") +
  #   labs(
  #     title = "Distribution on the VAF for each segment in the simulated data",
  #     subtitle = paste0(Subtitle_short)
  #   )+
  #   facet_wrap(vars(karyotype, tau, segment_name), scales = "free_x", strip.position = "bottom") +
  #   geom_vline(aes(xintercept = peaks_1), color = "grey") +
  #   geom_vline(aes(xintercept = peaks_2), color = "grey") +
  #   theme_minimal() +
  #   theme(
  #   panel.background = element_rect(fill = "white", color = NA),  # White panel background
  #   plot.background = element_rect(fill = "white", color = NA),   # White plot background
  #   strip.background = element_rect(fill = "white", color = NA),  # White strip background
  #   strip.placement = "outside",   # Place facet labels outside
  #   axis.text.x = element_text(angle = 360, hjust = 1, color = "black", size = 8),  # Rotate and adjust x-axis text
  #   axis.ticks.x = element_line(color = "black"),  # Black x-axis ticks
  #   panel.spacing = unit(1, "lines"),  # Adjust space between facets
  #   strip.text.x = element_text(size = 10, color = "black"),  # Adjust and color strip text
  #   axis.line = element_line(color = "black"),  # Black axis lines
  #   axis.title.x = element_text(color = "black"),  # Black x-axis title
  #   axis.title.y = element_text(color = "black")   # Black y-axis title
  # )
  
  # #save plot of the simulated data in which we can see each single segment VAF distribution
  # ggsave("./plots/simulation_data.png", plot = plot_data, width = 12 + simulation_params$number_events, height = 10 + simulation_params$number_events + (simulation_params$number_events/1.3), limitsize = FALSE,   device = png) 
  # #simulation_params can be substituted in relation with simulation_data variables

##########################################################################################################à
  
  
# per ora commentato perché da errore
  # hist_data <- all_sim %>%
  #   count(tau = round(tau, 2)) %>%
  #   filter(n > 0)

  # taus <- all_sim %>%
  #   ggplot(mapping = aes(x = tau)) +
  #   geom_histogram(alpha = .5, position = "identity", fill = "skyblue4") +
  #   labs(x = "tau")+
  #   labs(
  #     title = paste0( str_wrap("Distribution of tau in the simulated data", width = 40 + sqrt(simulation_params$number_events)) ),
  #     subtitle = " " 
  #   )+
  #   xlim(0,1) +
  #   stat_bin( geom="text", aes(label = ifelse(..count.. > 0, ..count.., "")), vjust=-1) +
  #   stat_bin(geom = "text", aes(label = ifelse(..count.. > 0, paste0("bold(", sprintf("%.2f", ..x..), ")"), "")), vjust=2, parse = TRUE)  +
  #   geom_hline(data = hist_data, aes(yintercept = n), linetype = "dashed", color = "grey")
  
    
  segment_counts <- all_sim %>%
    group_by(tau) %>%
    summarise(num_segments = n_distinct(segment_id)) %>%
    ungroup()
  
  # plot for validation tau parameters of the inference ##############################################################
  tau_segments_plot <- ggplot(segment_counts, aes(x = factor(round(tau, 2)), y = num_segments, fill = factor(tau))) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(x = "Tau", y = "Number of Segments") +
    labs(
      title = paste0( str_wrap("Number of Unique Segments Associated with Each Tau Value", width = 40 + K + (simulation_params$number_events)) ),
      subtitle = " "
    ) +
    theme_minimal()
      
  segment_counts <- all_sim %>%
    group_by(karyotype) %>%
    summarise(num_segments = n_distinct(segment_id)) %>%
    ungroup()
  
  karyo_segments_plot <- ggplot(segment_counts, aes(x = factor(karyotype), y = num_segments, fill = factor(karyotype))) +
    scale_fill_manual(values = c("skyblue", "orange", "lightgreen")) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(x = "Karyotype", y = "Number of Segments") +
    labs(
      title = paste0( str_wrap("Number of Unique Segments Associated with Each Karyotype configuration", width = 40 + (simulation_params$number_events)) ),
      subtitle = " "
    ) +
    theme_minimal()


  # EXTERNAL METRIC 
    #obtain score of simulation accuracy EXTERNAL METRICS (with known ground truth)
    #MAE evaluated on the tau values associated to each segment compare original with predicted/arrigned
    #RI compare the original partition of segments with respect to the tau they are associated with and the partition obtained from the tau assigned through the model
    accepted_mutations <- readRDS("results/accepted_mutations.rds") #Reload as I modify it for visualization purposes. Give as input do not call it from inside the function for the package
    
    names_tau <- paste("tau[", 1:K, "]", sep = "")
    tau_inferred <- res$draws(names_tau, format = "matrix")
    tau_inferred_map <- lapply(1:ncol(tau_inferred), function(i) {mode(tau_inferred[,i])} ) %>% unlist() 

    all_differences = c()
    identity_matrix_RI = c()
    for (i in seq_along(unique(accepted_mutations$segment_id))) {
        segment <- unique(accepted_mutations$segment_id)[i] #potrei mettere direttamente i 
        seg_modified <- segment %>%  str_split(" ")
        segment_number <- seg_modified[[1]][2]
        tau_original <- unique(accepted_mutations %>% filter(segment_id == segment)%>% select(tau)) %>% pull(tau)
        
        names_weights <- paste("w[",segment_number,",", 1:K, "]", sep = "")  #regex_pars = c("w")
        weights_inferred <- res$draws(names_weights, format = "matrix")
        weights_inferred_median <- lapply(1:ncol(weights_inferred), function(i) {median(weights_inferred[,i])} ) %>% unlist() 

        tau_index_assigned <- which.max(weights_inferred_median)
        tau_inferred_assigned =  tau_inferred_map[which.max(weights_inferred_median)]

        difference <- abs(tau_inferred_assigned - tau_original) # MAE - meglio MSE?
        all_differences <- c(all_differences, difference)
        identity_matrix_RI = c(identity_matrix_RI, tau_index_assigned) #extract the vector of taus (unordered) to which the ordered segments are assigned
    }

    MAE <- mean(all_differences)
    saveRDS(MAE, paste0("results/MAE_",K,".rds"))
    real_assignment <- accepted_mutations %>%
                                  group_by(segment_id) %>%
                                  summarize(tau = first(tau)) %>%
                                  arrange(match(segment_id, unique(accepted_mutations$segment_id))) %>%
                                  pull(tau) # here real_assignment is the true tau associated to the segment
    real_assignment <- rank(real_assignment, ties.method = "min") # here is the "class" to which it is associated --> the first event, the second ecc

    model_assignment <- identity_matrix_RI        #i,2,3 ecc 
    RI <- rand.index(real_assignment,model_assignment)
    #ARI <- adj.rand.index(real_assignment,model_assignment)  #NaN
    saveRDS(RI, paste0("results/RI_",K,".rds"))
    #saveRDS(ARI, paste0("results/ARI_",K,".rds"))

    MAE = round(MAE, 3)
    RI = round(RI, 3)

    final_plot <- (tau_segments_plot|karyo_segments_plot ) / plot_filtered_data /  (areas_tau | ppc) / intervals_weigths_per_tau / intervals_compare / (mean_compare|max_compare|min_compare|median_compare) +
      plot_layout(widths = c(8, 6, 6, 8, 8, 8), heights = c(10 + sqrt(simulation_params$number_events) , 15 + simulation_params$number_events + (simulation_params$number_events/1.3) , 15 + simulation_params$number_events + (simulation_params$number_events/1.3) + K, 15 + simulation_params$number_events + (simulation_params$number_events/1.3) + K*2, 8 + (simulation_params$number_events/1.3),  8 + (simulation_params$number_events/2))) +
      plot_annotation(
        title = paste0("Simulation with ", simulation_params$number_clocks," clocks, ", simulation_params$number_events, " segments, epsilon = ", simulation_params$epsilon, " purity = ", simulation_params$purity ),
        subtitle = " ",
        caption = paste0("MAE (tau real vs tau predicted) = ",MAE,"   RI (of real partition Vs model partition) = ",RI)
      ) & theme(text = element_text(size = 12+sqrt(simulation_params$number_events)), plot.title = element_text(size = 15+sqrt(simulation_params$number_events)), plot.subtitle = element_text(size = 12+sqrt(simulation_params$number_events)), axis.text = element_text(size = 12 + sqrt(simulation_params$number_events)), plot.caption = element_text(size = 10 + sqrt(simulation_params$number_events)))
    
    # ggsave("./plots/plots_inference.pdf", plot = final_plot, width = 12, height = 24)
    # ggsave("./plots/plots_inference.png", width = 12, height = 24, device = "png")
    
  return(final_plot)
}



library(ggplot2)
library(dplyr)
library(tidyr)
#justify why I am using the median of tau MAP 


#' plotting_cluster_partition Function
#'
#' This function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param res result
#' @param K number of mixture components
#' @param VALIDATION boolean
#' @export
#' @examples
#' plotting_cluster_partition()
 
plotting_cluster_partition <- function(res, K, VALIDATION = FALSE){

    accepted_mutations = readRDS("results/accepted_mutations.rds") # Adjust input as needed
    
    # Prepare data as for the MAE computation 
    names_tau <- paste("tau[", 1:K, "]", sep = "")
    tau_inferred <- res$draws(names_tau, format = "matrix")
    tau_inferred_map <- lapply(1:ncol(tau_inferred), function(i) {mode(tau_inferred[,i])} ) %>% unlist() 

    identity_matrix_RI = c()

    # Loop through each segment and determine its tau assignment
    for (i in seq_along(unique(accepted_mutations$segment_id))) {
        segment <- unique(accepted_mutations$segment_id)[i]
        seg_modified <- str_split(segment, " ")
        segment_number <- seg_modified[[1]][2]

        names_weights <- paste("w[",segment_number,",", 1:K, "]", sep = "")  
        weights_inferred <- res$draws(names_weights, format = "matrix")
        weights_inferred_map <- lapply(1:ncol(weights_inferred), function(i) {mode(weights_inferred[,i])} ) %>% unlist()

        tau_index_assigned <- which.max(weights_inferred_map)
        identity_matrix_RI <- c(identity_matrix_RI, tau_index_assigned)
    }

    # Prepare the data for ggplot
    segments <- unique(accepted_mutations$segment_id)
    plot_data <- data.frame(
        segment = segments,
        tau_index_assigned = identity_matrix_RI,
        tau_inferred_map = tau_inferred_map[identity_matrix_RI]
    )

    # Prepare identity matrix data
    identity_matrix_df <- data.frame(
        segment = rep(segments, each = K),
        tau_index = rep(1:K, times = length(segments)),
        value = as.numeric(rep(identity_matrix_RI, each = K) == rep(1:K, times = length(segments))),
        tau_inferred_map = rep(tau_inferred_map, times = length(segments))  # Correct the length mismatch
    )

    # Plot identity matrix with RColorBrewer palette
    library(RColorBrewer)

    # Plot identity matrix with RColorBrewer palette
    p <- ggplot(identity_matrix_df, aes(x = tau_index, y = segment)) +  # Change x to tau_index instead of tau_inferred_map
        geom_tile(aes(fill = factor(value)), color = "white") +
        scale_fill_brewer(palette = "Set1") + # Use a predefined palette
        labs(x = "Tau Index", y = "Segment", fill = "Membership") +
        theme_minimal(base_size = 14) +
        theme(
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            legend.background = element_rect(fill = "white"),
            legend.box.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = 12, angle = 0, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.ticks = element_line(color = "grey70"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(color = "black"),
            plot.subtitle = element_text(color = "black")
        ) +
        scale_x_continuous(breaks = 1:K)  # Set x-axis to tau indices

    return(p) # Return the ggplot object
}





#' plotting_model_selection Function
#'
#' This function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param results result
#' @export
#' @examples
#' plotting_model_selection()

plotting_model_selection <- function(results){

  model_selection <- results$model_selection_tibble
  best_K <- results$best_K
  k_max = nrow(model_selection)

# BIC plot
bic_plot <- ggplot(data = model_selection, aes(x = K, y = BIC)) + 
    geom_line(aes(colour = "BIC Line"), size = 1) + 
    geom_point(aes(colour = "BIC Point"), size = 3) +
    geom_point(aes(x = best_K, y = BIC[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[BIC == min(BIC)], y = min(BIC), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("BIC Line" = "steelblue", 
                                   "BIC Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  # Tema minimalista
    theme(legend.position = "top",  # Sposta la legenda sopra il grafico
          legend.title = element_blank(),  # Rimuove il titolo della legenda
          panel.grid.minor = element_blank(),  # Rimuove le griglie minori
          panel.grid.major.x = element_blank(),  # Rimuove le griglie verticali
          axis.title.x = element_text(margin = margin(t = 10)),  # Spazio per l'asse X
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))  # Spazio per l'asse Y


# AIC Plot
aic_plot <- ggplot(data = model_selection, aes(x = K, y = AIC)) + 
    geom_line(aes(colour = "AIC Line"), size = 1) + 
    geom_point(aes(colour = "AIC Point"), size = 3) +
    geom_point(aes(x = best_K, y = AIC[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[AIC == min(AIC)], y = min(AIC), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("AIC Line" = "steelblue", 
                                   "AIC Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  
    theme(legend.position = "top",  
          legend.title = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.grid.major.x = element_blank(),  
          axis.title.x = element_text(margin = margin(t = 10)),  
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))

# LOO Plot
loo_plot <- ggplot(data = model_selection, aes(x = K, y = LOO)) + 
    geom_line(aes(colour = "LOO Line"), size = 1) + 
    geom_point(aes(colour = "LOO Point"), size = 3) +
    geom_point(aes(x = best_K, y = LOO[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[LOO == min(LOO)], y = min(LOO), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("LOO Line" = "steelblue", 
                                   "LOO Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  
    theme(legend.position = "top",  
          legend.title = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.grid.major.x = element_blank(),  
          axis.title.x = element_text(margin = margin(t = 10)),  
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))

# Log Likelihood Plot
log_lik_plot <- ggplot(data = model_selection, aes(x = K, y = Log_lik)) + 
    geom_line(aes(colour = "Log likelihood Line"), size = 1) + 
    geom_point(aes(colour = "Log likelihood Point"), size = 3) +
    geom_point(aes(x = best_K, y = Log_lik[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[Log_lik == max(Log_lik)], y = max(Log_lik), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("Log likelihood Line" = "steelblue", 
                                   "Log likelihood Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  
    theme(legend.position = "top",  
          legend.title = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.grid.major.x = element_blank(),  
          axis.title.x = element_text(margin = margin(t = 10)),  
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))



# ICL plot
ICL_plot <- ggplot(data = model_selection, aes(x = K, y = ICL)) + 
    geom_line(aes(colour = "ICL Line"), size = 1) + 
    geom_point(aes(colour = "ICL Point"), size = 3) +
    geom_point(aes(x = best_K, y = ICL[K == best_K], colour = "Best K Point"), size = 4) +
    geom_point(aes(x = K[ICL == min(ICL)], y = min(ICL), colour = "Minimum Point"), size = 4) +
    scale_colour_manual(name = "Legend",
                        values = c("ICL Line" = "steelblue", 
                                   "ICL Point" = "steelblue", 
                                   "Best K Point" = "firebrick",
                                   "Minimum Point" = "forestgreen")) +
    theme_minimal(base_size = 15) +  # Tema minimalista
    theme(legend.position = "top",  # Sposta la legenda sopra il grafico
          legend.title = element_blank(),  # Rimuove il titolo della legenda
          panel.grid.minor = element_blank(),  # Rimuove le griglie minori
          panel.grid.major.x = element_blank(),  # Rimuove le griglie verticali
          axis.title.x = element_text(margin = margin(t = 10)),  # Spazio per l'asse X
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.margin = margin(5, 20, 5, 5))  # Spazio per l'asse Y


# Update the BIC plot
bic_plot <- bic_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(BIC, 2)), nudge_y = -10, size = 3, fill = "white", color = "black") 
# Update the AIC plot
aic_plot <- aic_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(AIC, 2)), nudge_y = -2, size = 3, fill = "white", color = "black")
# Update the LOO plot
loo_plot <- loo_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(LOO, 2)), nudge_y = 50, size = 3, fill = "white", color = "black") 
# Update the Log Likelihood plot
log_lik_plot <- log_lik_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(Log_lik, 2)), nudge_y = -0.01, size = 3, fill = "white", color = "black") 
# Update the ICL plot
ICL_plot <- ICL_plot +
  scale_x_continuous(breaks = seq(min(model_selection$K), max(model_selection$K), by = 1)) +
  geom_label(aes(label = round(ICL, 2)), nudge_y = -10, size = 3, fill = "white", color = "black") 


  
  S <- length(unique(results$accepted_mutations$segment_id))
  Subtitle <- rep(NA, times=k_max)
  
  for (k in 1:k_max) {
    n_parameters <- k+(k*S)+2
    Subtitle[[k]] <- paste0(" ", k, ": ", n_parameters," parameters")
  }
  
  Subtitle <- paste(Subtitle, collapse = "\n")
  
  

# Combinazione dei grafici con layout e annotazioni
model_selection_plot <- (bic_plot | aic_plot) / (loo_plot | log_lik_plot) /(ICL_plot) +
  plot_annotation(
    title = "Model Selection Graphs: Score vs Number of Clusters",
    subtitle = paste0("Correspondence between number of clusters and number of parameters: \n", Subtitle),
    caption = "Source: Your Data",
    theme = theme_minimal(base_size = 12) +  # Tema minimalista con base 12
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
        plot.caption = element_text(size = 8, hjust = 0.5, margin = margin(t = 10)),
        plot.background = element_rect(fill = "white", color = NA),  # Sfondo bianco senza bordo
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray80", size = 0.2),  # Griglie sottili e discrete
        panel.grid.minor = element_blank()
      )
  ) & 
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(size = 12, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 10, face = "italic", margin = margin(b = 10)),
    axis.text = element_text(size = 8),
    plot.caption = element_text(size = 8)
  )



  return(model_selection_plot)
}





#' plotting_elbo Function
#'
#' This function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param K number of clusters
#' @export
#' @examples
#' plotting_elbo()

plotting_elbo <- function(K){

  elbo_vs_iter <- list()
  for (k in 1:K){
        elbo_vs_iter[[k]] <- readRDS(paste0("./elbo_vs_iterations_",k,".rds"))+
                    labs(
                      title =  str_wrap( paste0("Elbo VS iterations with ",k," components"))  # , width = 30 + K*2
                    )
  }
  elbo_vs_iter_plot <- gridExtra::grid.arrange(grobs = elbo_vs_iter, ncol=K) #add global title
  return(elbo_vs_iter_plot)

}



#' plotting_fit Function
#'
#' For non validation setting, this function plot the essential information to be seen in order to evaluate the performance of the tool out of a simulation where real data or simulated data to be inferred are known.
#' @param res results
#' @param input_data input_data
#' @param all_sim all_sim
#' @param K K
#' @export
#' @examples
#' plotting_fit()

plotting_fit <- function(res, input_data, all_sim, K){

  draws <- res$draws(format = "matrix")
  names <- paste("tau[", 1:K, "]", sep = "")

  areas_tau <- mcmc_areas(
    draws,
    pars = names,
    prob = 0.8, # 80% intervals
    prob_outer = 0.95, # 99%
    point_est = "median"
  )+
    labs(
      title = "Approximate Posterior distributions",
      subtitle = "With median and 80% and 95% intervals"
    )+
    xlim(0, 1) # + scale_x_continuous(breaks = c(1:5), labels = c("A", "B", "C", "D", "E"))

  color_scheme_set("blue")
  intervals_weigths_per_tau <- list()
  for (k in 1:K){
        names_weights <- paste("w[",1:input_data$S,",", k, "]", sep = "") 
        intervals_weigths_per_tau[[k]] <- mcmc_intervals(draws, pars = names_weights, point_est = "median", prob = 0.8, prob_outer = 0.95)+
                    labs(
                      title =  str_wrap( paste0("Posterior distributions of the weigths for tau ",k), width = 30 + K + sqrt(input_data$S)),
                      subtitle = "With median and 80% and 95% intervals"
                    )
  }
  intervals_weigths_per_tau <- gridExtra::grid.arrange(grobs = intervals_weigths_per_tau, ncol=K) #add global title

  #not using it, just need to add the general title to avoid repetitions
  intervals <- mcmc_intervals(draws, regex_pars = c("w"), point_est = "median", prob = 0.8, prob_outer = 0.95)+
    labs(
      title = "Posterior distributions of the weigths for each tau inferred",
      subtitle = "With median and 80% and 95% intervals"
    )
  #

  #posterior predictive check
  stanfit <- rstan::read_stan_csv(res$output_files())
  y_rep <- as.matrix(stanfit, pars = "NV_pred")
  y = input_data$NV
  #distribution of replicated data vs real data
  ppc <- ppc_dens_overlay(
    y = input_data$NV,
    yrep = y_rep) +
    labs(
      title = "Posterior Predictive checks with 10000 iterations",
      subtitle = "Density distribution"
    )
  #empirical cumulative distribution
  ecdf_compare <-ppc_ecdf_overlay(y,y_rep)

  #predictive intervals vs observed values
  intervals_compare <- ppc_intervals(y,y_rep) +
    labs(
      subtitle = " Observed data vs Predicted values for each segment",
    )
  mean_compare <- ppc_stat(y,y_rep, stat = 'mean')
  max_compare <- ppc_stat(y,y_rep, stat = 'max')
  min_compare <- ppc_stat(y,y_rep, stat = 'min')
  median_compare <- ppc_stat(y,y_rep, stat = 'median')
  
  accepted_mutations <- readRDS("results/accepted_mutations.rds") #Give as input do not call it from inside the function for the package
  
  Subtitle <- vector("list", (length(unique(accepted_mutations$segment_id))+1))
  Subtitle[[1]]  <- paste0("Number of mutations per segment: ")
  num_mutations_all <- c()
    for (i in seq_along(unique(accepted_mutations$segment_id))) {
    segment <- unique(accepted_mutations$segment_id)[i]
    num_mutations_single <- nrow(accepted_mutations %>% filter(segment_id == segment))
    num_mutations_all <- c(num_mutations_all, num_mutations_single)
    Subtitle[[i+1]] <- paste0(segment, "=", num_mutations_single," ")
  }
  Subtitle <- paste(Subtitle, collapse = "   ")
  cat(Subtitle)

  mean_mut <- mean(num_mutations_all)
  max_mut <- max(num_mutations_all)
  min_mut <- min(num_mutations_all)

  Subtitle_short <- paste0("Mutations per segment: average =", mean_mut, ",  min = ", min_mut, ", max = ", max_mut )
  
  print(paste0("accepted mutation plot filtered data", dim(accepted_mutations)))


  plot_filtered_data <- accepted_mutations %>%
  ggplot(mapping = aes(x = NV / DP, fill = segment_id)) +
  geom_histogram(alpha = .5, position = "identity", bins=100) +
  labs(x = "VAF")+
  labs(
    title = paste0( str_wrap("Histogram of the VAF spectrum, per segment, resulting from the simulation (only the data used in the inference after the filtering step are plotted here)", width = 90 + K + (input_data$S) ) ),
    subtitle = paste0( str_wrap(Subtitle_short, width = 90 + K + (input_data$S) ) )
  )+
  facet_wrap(vars(karyotype))

  vline_positions <- data.frame(
  segment_id = paste0("segment ",1:input_data$S),
  peaks_1 = input_data$peaks[,1],
  peaks_2 = input_data$peaks[,2]
  )
  accepted_mutations <- merge(accepted_mutations, vline_positions, by = "segment_id")

  # print(paste0("accepted mutations complete", dim(accepted_mutations)))

  # #plot filter data to save separately and compare vaf distribution before and after the filtering step, I can remove maybe later
  # plot_filtered_data_complete <- accepted_mutations %>%
  # ggplot(mapping = aes(x = NV / DP, fill = segment_id)) +
  # geom_histogram(alpha = 0.6, position = "identity", color = "black") + # Add a black border to bars
  # facet_wrap(vars(karyotype, segment_id), scales = "free_x", strip.position = "bottom") +
  # geom_vline(aes(xintercept = peaks_1), color = "grey", linetype = "dashed") +  # Dashed lines for peaks
  # geom_vline(aes(xintercept = peaks_2), color = "grey", linetype = "dashed") + 
  # theme_minimal() +
  # theme(
  #   panel.background = element_rect(fill = "white", color = NA),  
  #   plot.background = element_rect(fill = "white", color = NA),   
  #   strip.background = element_rect(fill = "white", color = NA),  
  #   strip.placement = "outside",   
  #   axis.text.x = element_text(angle = 360, hjust = 1, color = "black", size = 8),  
  #   axis.ticks.x = element_line(color = "black"),  
  #   panel.grid.major = element_line(color = "grey90"), # Subtle grid lines
  #   panel.spacing = unit(1, "lines"),  
  #   strip.text.x = element_text(size = 10, color = "black"),  
  #   axis.line = element_line(color = "black"),  
  #   axis.title.x = element_text(color = "black", size = 12),  # Bigger font for axis titles
  #   axis.title.y = element_text(color = "black", size = 12),   
  #   plot.title = element_text(size = 16, face = "bold"),  # Bold title for better readability
  #   plot.subtitle = element_text(size = 12, face = "italic")  # Italic subtitle
  # ) +
  # labs(
  #   x = "VAF",
  #   y = "Count",
  #   title = "Histogram of the VAF spectrum, per segment",
  #   subtitle = "Only the data used in the inference after the filtering step are plotted"
  # )+
  #   xlim(0, 1)

  # ggsave("./plots/plot_filtered_data_complete.png", plot = plot_filtered_data_complete, width = 8 + input_data$S, height = 6 + input_data$S + (input_data$S/1.3), limitsize = FALSE,   device = png) 



  Subtitle <- vector("list", (length(unique(all_sim$segment_id))+1))
  Subtitle[[1]]  <- paste0("Number of mutations per segment: ")
  num_mutations_all_segments <- c()

    for (i in seq_along(unique(all_sim$segment_id))) {
    segment <- unique(all_sim$segment_id)[i]
    num_mutations <- nrow(all_sim %>% filter(segment_id == segment))
    num_mutations_all_segments <- c(num_mutations_all_segments, num_mutations)
    Subtitle[[i+1]] <- paste0(segment, "=", num_mutations," ")
  }
  
  Subtitle <- paste(Subtitle, collapse = "   ")
  cat(Subtitle)

  mean_mut <- mean(num_mutations_all_segments)
  max_mut <- max(num_mutations_all_segments)
  min_mut <- min(num_mutations_all_segments)

  Subtitle_short <- paste0("Average number of mutations per segment: ", mean_mut, "  Minimum number of mutations per segment: ", min_mut, "  Maximum number of mutations per segment: ", max_mut )


  vline_positions <- data.frame(
  segment_name = paste0("segment ",1:input_data$S),
  peaks_1 = input_data$peaks[,1],
  peaks_2 = input_data$peaks[,2]
  )
  all_sim_ <- merge(all_sim, vline_positions, by = "segment_name")

  # plot_data <- simulation_data_plot %>% 
  #   ggplot(mapping = aes(x = NV / DP, fill = segment_name)) +
  #   geom_histogram(alpha = .5, position = "identity") +
  #   facet_wrap(vars(karyotype, segment_name), scales = "free_x", strip.position = "bottom") +
  #   labs(
  #     title = "Distribution on the VAF for each segment in the simulated data",
  #     subtitle = paste0(Subtitle_short)
  #   )+
  #   geom_vline(aes(xintercept = peaks_1), color = "grey") +
  #   geom_vline(aes(xintercept = peaks_2), color = "grey") +
  #   theme_minimal() +
  #   theme(
  #   panel.background = element_rect(fill = "white", color = NA),  # White panel background
  #   plot.background = element_rect(fill = "white", color = NA),   # White plot background
  #   strip.background = element_rect(fill = "white", color = NA),  # White strip background
  #   strip.placement = "outside",   # Place facet labels outside
  #   axis.text.x = element_text(angle = 360, hjust = 1, color = "black", size = 8),  # Rotate and adjust x-axis text
  #   axis.ticks.x = element_line(color = "black"),  # Black x-axis ticks
  #   panel.spacing = unit(1, "lines"),  # Adjust space between facets
  #   strip.text.x = element_text(size = 10, color = "black"),  # Adjust and color strip text
  #   axis.line = element_line(color = "black"),  # Black axis lines
  #   axis.title.x = element_text(color = "black"),  # Black x-axis title
  #   axis.title.y = element_text(color = "black")   # Black y-axis title
  # )
  
  #save plot of the simulated data in which we can see each single segment VAF distribution
  # ggsave("./plots/simulation_data.png", plot = plot_data, width = 12 + input_data$S, height = 10 + input_data$S + (simulation_params/1.3), limitsize = FALSE,   device = png) 
  #simulation_params can be substituted in relation with simulation_data variable   


   accepted_mutations <- readRDS("results/accepted_mutations.rds") #Reload as I modify it for visualization purposes. Give as input do not call it from inside the function for the package
    
    print(paste0("accepted_mutations ", accepted_mutations))

    names_tau <- paste("tau[", 1:K, "]", sep = "")
    tau_inferred <- res$draws(names_tau, format = "matrix")
    tau_inferred_map <- lapply(1:ncol(tau_inferred), function(i) {mode(tau_inferred[,i])} ) %>% unlist() 

    all_differences = c()
    identity_matrix_RI = c()
    for (i in seq_along(unique(accepted_mutations$segment_id))) {
        segment <- unique(accepted_mutations$segment_id)[i] #potrei mettere direttamente i 
        seg_modified <- segment %>%  str_split(" ")
        segment_number <- seg_modified[[1]][2]
      
                print(paste0("segment_number",segment_number))

        names_weights <- paste("w[",segment_number,",", 1:K, "]", sep = "")  #regex_pars = c("w")
        weights_inferred <- res$draws(names_weights, format = "matrix")
        weights_inferred_median <- lapply(1:ncol(weights_inferred), function(i) {median(weights_inferred[,i])} ) %>% unlist() 

        tau_index_assigned <- which.max(weights_inferred_median)
        tau_inferred_assigned =  tau_inferred_map[which.max(weights_inferred_median)]
        identity_matrix_RI = c(identity_matrix_RI, tau_index_assigned) #extract the vector of taus (unordered) to which the ordered segments are assigned
    }

    model_assignment <- identity_matrix_RI        #i,2,3 ecc 

    final_plot <-plot_filtered_data /  (areas_tau | ppc) / intervals_weigths_per_tau / intervals_compare / (mean_compare|max_compare|min_compare|median_compare) +
      plot_layout(widths = c( 6, 6, 8, 8, 8), heights = c(15 + input_data$S + (input_data$S/1.3) , 15 + input_data$S + (input_data$S/1.3) + K, 15 + input_data$S + (input_data$S/1.3) + K*2, 8 + (input_data$S/1.3),  8 + (input_data$S/2))) +
      plot_annotation(
        title = paste0("Number of events", input_data$S),
        subtitle = " ",
        caption = ""
      ) & theme(text = element_text(size = 12+sqrt(input_data$S)), plot.title = element_text(size = 15+sqrt(input_data$S)), plot.subtitle = element_text(size = 12+sqrt(input_data$S)), axis.text = element_text(size = 12 + sqrt(input_data$S)), plot.caption = element_text(size = 10 + sqrt(input_data$S)))
    

  return(final_plot)
}
