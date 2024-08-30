
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
  intervals <- mcmc_intervals(draws, regex_pars = c("w"), point_est = "mean", prob = 0.8, prob_outer = 0.95)+
    labs(
      title = "Posterior distributions",
      subtitle = "with mean and 80% and 95% intervals"
    )


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

    for (i in seq_along(unique(accepted_mutations$segment_id))) {
    segment <- unique(accepted_mutations$segment_id)[i]
    num_mutations <- nrow(accepted_mutations %>% filter(segment_id == segment))
    Subtitle[[i+1]] <- paste0(segment, "=", num_mutations," ")
  }
  
  Subtitle <- paste(Subtitle, collapse = "   ")
  

  mean_mut <- mean(num_mutations)
  max_mut <- max(num_mutations)
  min_mut <- min(num_mutations)

  Subtitle_short <- paste0("Average number of mutations per segment: ", mean_mut, "  Minimum number of mutations per segment: ", min_mut, "  Maximum number of mutations per segment: ", max_mut )


  
  accepted_mutations = accepted_mutations %>% mutate (tau = round(tau, 2))
    plot_filtered_data <- accepted_mutations %>%
    ggplot(mapping = aes(x = NV / DP, fill = as.factor(segment_id))) +
    geom_histogram(alpha = .5, position = "identity") +
    labs(x = "VAF")+
    labs(
      title = "Histogram of the VAF spectrum, per segment, resulting from the simulation (only the data used in the inference after the filtering step are plotted here)",
      subtitle = paste0( str_wrap(Subtitle_short, width = 80 ) )
    )+
    facet_wrap(vars(karyotype, tau))

  
  
  
  
  
  
  hist_data <- all_sim %>%
    count(tau = round(tau, 2)) %>%
    filter(n > 0)

  taus <- all_sim %>%
    ggplot(mapping = aes(x = tau)) +
    geom_histogram(alpha = .5, position = "identity", fill = "skyblue4") +
    labs(x = "tau")+
    labs(
      title = paste0( str_wrap("Distribution of tau in the simulated data", width = 40) ),
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
      title = paste0( str_wrap("Number of Unique Segments Associated with Each Tau Value", width = 40) ),
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
      title = "Number of Unique Segments Associated with Each Karyotype configuration",
      subtitle = " "
    ) +
    theme_minimal()
  
    




    #obtain score of simulation accuracy 
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


    final_plot <- (tau_segments_plot|karyo_segments_plot ) / plot_filtered_data /  (areas_tau | intervals) / ppc / intervals_compare / (mean_compare|max_compare|min_compare|median_compare) +
      plot_layout(widths = c(8, 6, 6, 8, 8, 8), heights = c(10 + sqrt(simulation_params$number_events) , 15 + simulation_params$number_events + (simulation_params$number_events/1.3) , 15 + simulation_params$number_events + (simulation_params$number_events/1.3) + K, 8 + (simulation_params$number_events/1.3), 8 + (simulation_params$number_events/1.3),  8 + (simulation_params$number_events/2))) +
      plot_annotation(
        title = paste0("Simulation with ", simulation_params$number_clocks," clocks, ", simulation_params$number_events, " segments, epsilon = ", simulation_params$epsilon, " purity = ", simulation_params$purity ),
        subtitle = " ",
        caption = " "
      ) & theme(text = element_text(size = 14+sqrt(simulation_params$number_events)), plot.title = element_text(size = 18+sqrt(simulation_params$number_events)), plot.subtitle = element_text(size = 14+sqrt(simulation_params$number_events)), axis.text = element_text(size = 14 + sqrt(simulation_params$number_events)), plot.caption = element_text(size = 8))
    
    
    
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

