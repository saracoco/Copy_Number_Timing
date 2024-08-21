//
//

data {
  int N; //N = nrow(accepted_mutations)
  array[N] int NV;
  array[N] int DP;

  array[2] real peaks;
  // real beta_dispersion;
  // real purity;
  // int multiplicities[n_groups];
}

parameters {
  simplex[2] omega; //weigths of the binomial mixture for the mutations, parametrized by (#reads, expected peak)
}

model {
  vector[2] contributions;
  // priors
  omega ~ dirichlet(rep_vector(2.0, 2)); //alpha = 2.0 

  // likelihood
  for (i in 1:N) {
    for (k in 1:2) {
      contributions[k] = log(omega[k]) + binomial_lpmf(NV[i] | DP[i], peaks[k]);
    }
    target += log_sum_exp(contributions);
  }
}


generated quantities {
  array[N] int comp_binomial; // Binomial component identifier
  vector[N] y_rep;  // Predicted NV

  vector[N] log_lik; // log-likelihood for each data point

  vector[2] contributions;
  for (i in 1:N) {
      for (k in 1:2) {
      contributions[k] = log(omega[k]) + binomial_lpmf(NV[i] | DP[i], peaks[k]);
    }
    log_lik[i] = log_sum_exp(contributions);
    
    // Generate quantities for posterior predictive checks
    comp_binomial[i] = categorical_rng(omega);
    y_rep[i] = binomial_rng(DP[i], peaks[comp_binomial[i]]);
  }
}
