
#' plotting Function
#'
#' This function obtains the peaks to be used for CN timing in the input data.
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
    xlim(0, 1)

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

  
   
  
  accepted_mutations <- readRDS("results/accepted_mutations.rds") #giveas input do not call it from inside the function for the package
  
  Subtitle <- vector("list", (length(unique(accepted_mutations$segment_id))+1))
  Subtitle[[1]]  <- paste0("Number of mutations per segment: ")

    for (i in seq_along(unique(accepted_mutations$segment_id))) {
    segment <- unique(accepted_mutations$segment_id)[i]
    num_mutations <- nrow(accepted_mutations %>% filter(segment_id == segment))
    Subtitle[[i+1]] <- paste0(segment, "= ", num_mutations," ")
  }
  
  Subtitle <- paste(Subtitle, collapse = "   ")
  
  
  accepted_mutations = accepted_mutations %>% mutate (tau = round(tau, 2))
    plot_filtered_data <- accepted_mutations %>%
    ggplot(mapping = aes(x = NV / DP, fill = as.factor(segment_id))) +
    geom_histogram(alpha = .5, position = "identity") +
    labs(x = "VAF")+
    labs(
      title = "Histogram for the VAF spectrum, per segment, resulting from the simulation (only the data used in the inference after the filtering step are plotted here)",
      subtitle = paste0(Subtitle)
    )+
    facet_wrap(vars(karyotype, tau, segment_id))

  
  
  
  
  
  
  hist_data <- all_sim %>%
    count(tau = round(tau, 2)) %>%
    filter(n > 0)

  taus <- all_sim %>%
    ggplot(mapping = aes(x = tau)) +
    geom_histogram(alpha = .5, position = "identity", fill = "skyblue4") +
    labs(x = "tau")+
    labs(
      title = "Distribution of tau in the simulated data",
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
      title = "Number of Unique Segments Associated with Each Tau Value",
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
  
    

    final_plot <- (tau_segments_plot|karyo_segments_plot ) / plot_filtered_data /  (areas_tau | intervals) / (ppc | intervals_compare) / (mean_compare|max_compare|min_compare|median_compare) +
      plot_layout(widths = c(8, 6, 6, 8, 8), heights = c(8, 15 + simulation_params$number_events + (simulation_params$number_events/2) , 15 + simulation_params$number_events + (simulation_params$number_events/2) , 8 (simulation_params$number_events/2), 8 + (simulation_params$number_events/2))) +
      plot_annotation(
        title = paste0("Simulation with ", simulation_params$number_clocks," clocks, ", simulation_params$number_events, " segments, epsilon = ", simulation_params$epsilon, " purity = ", simulation_params$purity ),
        subtitle = " ",
        caption = " "
      ) & theme(text = element_text(size = 8), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), axis.text = element_text(size = 8), plot.caption = element_text(size = 5))
    
    
    
    # ggsave("./plots/plots_inference.pdf", plot = final_plot, width = 12, height = 24)
    # ggsave("./plots/plots_inference.png", width = 12, height = 24, device = "png")
    #



  return(final_plot)
}
