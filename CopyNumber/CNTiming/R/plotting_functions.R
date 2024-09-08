
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

  draws <- res$draws(format = "matrix")

  #color_scheme_set("red")
  #manage to get the K from the model fit directly rather than as input
  names <- paste("tau[", 1:K, "]", sep = "")

  areas_tau <- mcmc_areas(
    draws,
    pars = names,
    prob = 0.8, # 80% intervals
    prob_outer = 0.95, # 99%
    point_est = "mean"
  )+
    labs(
      title = "Approximate Posterior distributions",
      subtitle = "with mean and 80% and 95% intervals"
    )+
    xlim(0, 1) # + scale_x_continuous(breaks = c(1:5), labels = c("A", "B", "C", "D", "E"))



  color_scheme_set("blue")
  intervals_weigths_per_tau <- list()
  for (k in 1:K){
        names_weights <- paste("w[",1:simulation_params$number_events,",", k, "]", sep = "") 
        intervals_weigths_per_tau[[k]] <- mcmc_intervals(draws, pars = names_weights, point_est = "mean", prob = 0.8, prob_outer = 0.95)+
                    labs(
                      title =  str_wrap( paste0("Posterior distributions of the weigths for tau ",k), width = 30 + K + sqrt(simulation_params$number_events)),
                      subtitle = "with mean and 80% and 95% intervals"
                    )
  }
  intervals_weigths_per_tau <- gridExtra::grid.arrange(grobs = intervals_weigths_per_tau, ncol=K) #add global title

  
  #not using it, just need to add the general title to avoid repetitions
  intervals <- mcmc_intervals(draws, regex_pars = c("w"), point_est = "mean", prob = 0.8, prob_outer = 0.95)+
    labs(
      title = "Posterior distributions of the weigths for each tau inferred",
      subtitle = "with mean and 80% and 95% intervals"
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

  #Compare statistics
  #compute the estimated Bayesian p-value
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

  Subtitle_short <- paste0("Average number of mutations per segment: ", mean_mut, "  Minimum number of mutations per segment: ", min_mut, "  Maximum number of mutations per segment: ", max_mut )


  
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





  vline_positions <- data.frame(
  segment_id = paste0("segment ",1:input_data$S),
  peaks_1 = input_data$peaks[,1],
  peaks_2 = input_data$peaks[,2]
  )
  accepted_mutations <- merge(accepted_mutations, vline_positions, by = "segment_id")


  #plot filter data to save separately and compare vaf distribution before and after the filtering step, I can remove maybe later
  plot_filtered_data_complete <- accepted_mutations %>%
  ggplot(mapping = aes(x = NV / DP, fill = segment_id)) +
  geom_histogram(alpha = .5, position = "identity") +
  facet_wrap(vars(karyotype, tau, segment_id), scales = "free_x", strip.position = "bottom") +
  geom_vline(aes(xintercept = peaks_1), color = "grey") +
  geom_vline(aes(xintercept = peaks_2), color = "grey") +
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
  )+
  labs(x = "VAF")+
  labs(
    title = paste0( str_wrap("Histogram of the VAF spectrum, per segment, resulting from the simulation (only the data used in the inference after the filtering step are plotted here)", width = 90 + K + (simulation_params$number_events) ) ),
    subtitle = paste0( Subtitle_short )
  )

  ggsave("./plots/plot_filtered_data_complete.png", plot = plot_filtered_data_complete, width = 8 + simulation_params$number_events, height = 6 + simulation_params$number_events + (simulation_params$number_events/1.3), limitsize = FALSE,   device = png) 











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
  all_sim <- merge(all_sim, vline_positions, by = "segment_name")

  simulation_data_plot = all_sim %>% mutate (tau = round(tau, 2))
  plot_data <- simulation_data_plot %>% 
    ggplot(mapping = aes(x = NV / DP, fill = segment_name)) +
    geom_histogram(alpha = .5, position = "identity") +
    labs(
      title = "Distribution on the VAF for each segment in the simulated data",
      subtitle = paste0(Subtitle_short)
    )+
    facet_wrap(vars(karyotype, tau, segment_name), scales = "free_x", strip.position = "bottom") +
    geom_vline(aes(xintercept = peaks_1), color = "grey") +
    geom_vline(aes(xintercept = peaks_2), color = "grey") +
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
 


  
  





  hist_data <- all_sim %>%
    count(tau = round(tau, 2)) %>%
    filter(n > 0)

  taus <- all_sim %>%
    ggplot(mapping = aes(x = tau)) +
    geom_histogram(alpha = .5, position = "identity", fill = "skyblue4") +
    labs(x = "tau")+
    labs(
      title = paste0( str_wrap("Distribution of tau in the simulated data", width = 40 + sqrt(simulation_params$number_events)) ),
      subtitle = " " 
    )+
    xlim(0,1) +
    stat_bin( geom="text", aes(label = ifelse(..count.. > 0, ..count.., "")), vjust=-1) +
    stat_bin(geom = "text", aes(label = ifelse(..count.. > 0, paste0("bold(", sprintf("%.2f", ..x..), ")"), "")), vjust=2, parse = TRUE)  +
    geom_hline(data = hist_data, aes(yintercept = n), linetype = "dashed", color = "grey")
  
    

  segment_counts <- all_sim %>%
    group_by(tau) %>%
    summarise(num_segments = n_distinct(segment_id)) %>%
    ungroup()
  
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
  
    




    #obtain score of simulation accuracy    EXTERNAL METRICS (with known ground truth)
    #MAE evaluated on the tau values associated to each segment compare original with predicted/arrigned
    #RI compare the original partition of segments with respect to the tau they are associated with and the partition obtained from the tau assigned through the model
    accepted_mutations <- readRDS("results/accepted_mutations.rds") #Reload as I modify it for visualization purposes. Give as input do not call it from inside the function for the package
    
    names_tau <- paste("tau[", 1:K, "]", sep = "")
    tau_inferred <- res$draws(names_tau, format = "matrix")
    tau_inferred_median <- lapply(1:ncol(tau_inferred), function(i) {median(tau_inferred[,i])} ) %>% unlist() 

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
        tau_inferred_assigned =  tau_inferred_median[which.max(weights_inferred_median)]

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
                                  pull(tau)
    model_assignment <- identity_matrix_RI
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
    #



  return(final_plot)
}















plotting_model_selection <- function(results){

  model_selection <- results$model_selection_tibble
  best_K <- results$best_K
  k_max = nrow(model_selection)
  
  
  bic_plot <- ggplot(data = model_selection, aes(x = K, y = BIC)) + 
    geom_line(aes(colour = "BIC Line")) + 
    geom_point(aes(colour = "BIC Point")) +
    geom_point(aes(x = best_K, y = BIC[K == best_K], colour = "Best K Point")) +
    geom_point(aes(x = K[BIC == min(BIC)], y = min(BIC), colour = "Minimum Point")) +
    scale_colour_manual(name = "Legend",
                        values = c("BIC Line" = "blue", 
                                   "BIC Point" = "blue", 
                                   "Best K Point" = "red",
                                   "Minimum Point" = "green"),)
  
  aic_plot <- ggplot(data = model_selection, aes(x = K, y = AIC)) + 
    geom_line(aes(colour = "AIC Line")) + 
    geom_point(aes(colour = "AIC Point")) +
    geom_point(aes(x = best_K, y = AIC[K == best_K], colour = "Best K Point")) +
    geom_point(aes(x = K[AIC == min(AIC)], y = min(AIC), colour = "Minimum Point")) +
    scale_colour_manual(name = "Legend",
                        values = c("AIC Line" = "blue", 
                                   "AIC Point" = "blue", 
                                   "Best K Point" = "red",
                                   "Minimum Point" = "green"))
  

  loo_plot <- ggplot(data = model_selection, aes(x = K, y = LOO)) + 
    geom_line(aes(colour = "LOO Line")) + 
    geom_point(aes(colour = "LOO Point")) +
    geom_point(aes(x = best_K, y = LOO[K == best_K], colour = "Best K Point")) +
    geom_point(aes(x = K[LOO == min(LOO)], y = min(LOO), colour = "Minimum Point")) +
    scale_colour_manual(name = "Legend",
                        values = c("LOO Line" = "blue", 
                                   "LOO Point" = "blue", 
                                   "Best K Point" = "red",
                                   "Minimum Point" = "green"))
  
  
  log_lik_plot <- ggplot(data = model_selection, aes(x = K, y = Log_lik)) + 
    geom_line(aes(colour = "Log likelihood Line")) + 
    geom_point(aes(colour = "Log likelihood Point")) +
    geom_point(aes(x = best_K, y = Log_lik[K == best_K], colour = "Best K Point")) +
    geom_point(aes(x = K[Log_lik == max(Log_lik)], y = max(Log_lik), colour = "Minimum Point")) +
    scale_colour_manual(name = "Legend",
                        values = c("Log likelihood Line" = "blue", 
                                   "Log likelihood Point" = "blue", 
                                   "Best K Point" = "red",
                                   "Minimum Point" = "green"))
  
  S <- length(unique(results$accepted_mutations$segment_id))
  Subtitle <- rep(NA, times=k_max)
  
  for (k in 1:k_max) {
    n_parameters <- k+(k*S)+2
    Subtitle[[k]] <- paste0(" ", k, ": ", n_parameters," parameters")
  }
  
  Subtitle <- paste(Subtitle, collapse = "\n")
  
  
  
  model_selection_plot <- (bic_plot | aic_plot) / (loo_plot|log_lik_plot) +
  plot_annotation(
    title = "Model selection graphs: score vs number of clusters" ,
    subtitle = paste0 ("Correspondence between number of clusters and number of parameters:  \n", Subtitle),
    caption = " "
  ) & theme(text = element_text(size = 8), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), axis.text = element_text(size = 8), plot.caption = element_text(size = 8))
  
  return(model_selection_plot)
}

