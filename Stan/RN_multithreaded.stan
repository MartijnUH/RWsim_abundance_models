// Royle-Nichols (RN) model - multithreading

functions {
  /**
   * Return log probability of RN-model for a site
   * 
   * @param det      Detection status 
   * @param J        Number of surveys/ replicates
   * @param occ      Occupancy status 
   * @param max_n    Maximum site-abundance
   * @param lambda   Mean site-abundance
   * @param p        Per Capita detection probability
   *
   * @return         Log probability
   */
   
  vector rn_logp(int[] det, int occ, int max_n,
               real lambda, real p) {
    
    vector[max_n - occ + 1] lp; 
    if(occ == 0){
      lp[1] = poisson_log_lpmf(0 | log(lambda)) + 1;
    }
    else lp[1] = poisson_log_lpmf(1 | log(lambda)) + binomial_lpmf(det | 1, p);
    
    for (j in 2:(max_n - occ + 1)){
      lp[j] = poisson_log_lpmf(occ + j - 1 | log(lambda))
      + binomial_lpmf(det | 1, 1 - (1 - p)^(occ + j - 1) );
    }
    return lp;
  }
  
  real rn_lpmf(int[] det, int occ, int max_n,
               real lambda, real p) {
    
    vector[max_n - occ + 1] lp; 
    if(occ == 0){
      lp[1] = poisson_log_lpmf(0 | log(lambda)) + 1;
    }
    else lp[1] = poisson_log_lpmf(1 | log(lambda)) + binomial_lpmf(det | 1, p);
    
    for (j in 2:(max_n - occ + 1)){
      lp[j] = poisson_log_lpmf(occ + j - 1 | log(lambda))
      + binomial_lpmf(det | 1, 1 - (1 - p)^(occ + j - 1) );
    }
    return log_sum_exp(lp);
  }

  real partial_sum(int[] site,
                   int start, int end,
                   int[,] det, int[] occ,
                   int max_n, real lambda, real p) {
    real lp = 0;

    for (m in start:end)
      lp = lp + rn_lpmf(det[m] | occ[m], max_n, lambda, p);
    return lp;
  }
}

data {
  int<lower = 0> R;              // number of sites
  int<lower = 0> T;              // number of replications
  int<lower = 0> K;              // upper limit of the abundance
  int<lower = 0> y[R, T];        // detection/ non-detection for each site and replication

}

transformed data {
  int site[R] = rep_array(0, R); // Dummy for site index
  int<lower=0,upper=1> z[R];     // at least one detection
  for (i in 1:R) {
    z[i] = 0;
    if (sum(y[i]))
      z[i] = 1;
  }    
}

parameters {
  real<lower=0> lambda;          // Mean site-abundance
  real<lower=0, upper=1> p;      // Detection probability
}

model {
  int grainsize = 1;
  
  // Priors
  lambda ~ cauchy(0, 10);
  // A flat prior [0, 1] is implicitly used on p.
  
  // Partial likelihood
  target += reduce_sum(partial_sum, site, grainsize, y, z, K, lambda, p);
}

generated quantities {
  vector[R] log_lik;                // Pointwise log likelihood
  real<lower=0,upper=1> psi;        // Mean site-occupancy
  real<lower=0,upper=1> P;          // Species detection probability 
  int totalN;                       // Total population size
  array[R, T] int y_new;
  real fit = 0;                     // Test statistic actual data
  real fit_new = 0;                 // Test statistic new data
  real chat;                        // C-statistic (degree of overdisperion)
  int<lower=0, upper=1> pval = 0;   // Bayesian p-value for test statistic
  
  P = 1 - (1 - p)^poisson_rng(lambda);
  psi = 1 - exp(-lambda);
  
  {
    array[R] int N;                 // Abundance
    array[R, T] real eval;          // Expected values
    matrix[R, T] E;
    matrix[R, T] E_new;
    
    // Initialize counter, N, E and E_new
    N = rep_array(0, R);
    E = rep_matrix(0, R, T);
    E_new = rep_matrix(0, R, T);
    
    for (i in 1 : R) {
      
      // pointwise log likelihood contributions
      log_lik[i] = rn_lpmf(y[i,] | z[i], K, lambda, p);
          
      // latent N's
      vector[K - z[i] + 1] lp;
      lp = rn_logp(y[i,], z[i], K, lambda, p);
      N[i] = categorical_rng(softmax(lp)) - 1;
          
      for (t in 1 : T) {
        // expected detections
        eval[i, t] = 1 - (1 - p)^N[i];
        
        // Assess model fit using Chi-squared discrepancy
        // Compute fit statistic E for observed data
        E[i, t] = square(y[i, t] - eval[i, t]) / (eval[i, t] + 0.5);
        // Generate replicate data and
        // Compute fit statistic E_new for replicate data
        y_new[i, t] = bernoulli_rng(1 - (1 - p)^N[i]);
        E_new[i, t] = square(y_new[i, t] - eval[i, t])/ (eval[i, t] + 0.5);
      }
      // Summed fit statistics
      fit = fit + sum(E[i]);
      fit_new = fit_new + sum(E_new[i]);
    }
    
    // Bayesian P-values
    chat = fit/ (fit_new + 1e-9);
    if (fit_new > fit) pval = 1;
    
    // Total pop. size across all sampled sites
    totalN = sum(N[1 : R]);
  }
}
