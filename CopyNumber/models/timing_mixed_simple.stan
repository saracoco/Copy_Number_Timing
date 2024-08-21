// reparametrization with dirichlet alpha being mean value of w and nuber of observations of the prior
data {
  // input the karyotyoe to specify the formula for theta
  //
  int S; // number of segments
  int K; // clocks
  int N; // total number of mutations
  array[S] int karyotype; // list of karyotype associated wit the segments
  array[N] int seg_assignment; // segment_id assignment to the mutations
  array[S,2] real peaks;  //for each segment, S vectors of dim 2
  array[N] int NV;     // for all the segments
  array[N] int DP;
}


parameters {
  array[S] simplex[K] w;  // simplex[K] w[N] mixing proportions for each segment group;
  vector<lower=0,upper=1>[K] tau;     //clocks   ordered[K] tau  vector<lower=0,upper=1>[K]
  //vector[K] alpha;
  simplex[K] phi;
  real<lower=0> kappa;
}

transformed parameters{

  array[S,K,2] real<lower=0,upper=1> theta; //binomial mixing proportions // array[S,K] simplex[2] theta;

  for (s in 1:S){
   for (k in 1:K){
      if (karyotype[s] == 1) {
        theta[s,k,1] = (3 - 2*tau[k]) / (3 - tau[k]);      // 2:1
        theta[s,k,2] = tau[k] / (3 - tau[k]);
      }
      else {
        theta[s,k,1] = (2 - 2*tau[k]) / (2 - tau[k]);      // 2:0 - 2:2
        theta[s,k,2] = tau[k]/(2 - tau[k]);
      }
    }
  }

  vector[K] alpha = kappa * phi;

}



model {
  vector[K*2] contributions;

  // priors

  //alpha ~
  phi ~ dirichlet(rep_vector(1.0, K));;
  kappa ~ normal(5,2.5);
  // phi = expected value of w, kappa (minus K) = concentrazione della distribuzione / strength of the prior mean measured in number of prior observations.

  for (s in 1:S){
      w[s] ~ dirichlet(alpha);
  }
  tau ~ beta(2,2);



  //likelihood
      for (i in 1:N) {
        int c = 1;
        for (k in 1:K) {
         for (j in 1:2) {
          contributions[c] = log(w[seg_assignment[i],k]) + log(theta[seg_assignment[i],k,j]) + binomial_lpmf(NV[i] | DP[i], peaks[seg_assignment[i],j]);
          c += 1;
        }
      }
      target += log_sum_exp(contributions);
      //print(target)
    }
}


generated quantities {
  array[N] int comp_tau;    // Component identifier "comp"
  array[N] int comp_binomial; // Binomial component identifier
  vector[N] NV_pred;  // Predicted NV
  vector[2] theta_vector;

  vector[N] log_lik; // log-likelihood for each data point

  for (i in 1:N) {
    vector[K*2] contributions;
    int c = 1;
    for (k in 1:K) {
      for (j in 1:2) {
        contributions[c] = log(w[seg_assignment[i],k]) + log(theta[seg_assignment[i],k,j]) + binomial_lpmf(NV[i] | DP[i], peaks[seg_assignment[i],j]);
        c += 1;
      }
    }
    log_lik[i] = log_sum_exp(contributions);
    
    // Generate quantities for posterior predictive checks
    comp_tau[i] = categorical_rng(w[seg_assignment[i]]);
    theta_vector = to_vector(theta[seg_assignment[i], comp_tau[i]]);
    comp_binomial[i] = categorical_rng(theta_vector);
    NV_pred[i] = binomial_rng(DP[i], peaks[seg_assignment[i], comp_binomial[i]]);
  }
}

