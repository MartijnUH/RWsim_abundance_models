// N-mixture model - multithreading

functions {
  
  /**
   * Return log probability of Poisson Binomial Mixture
   *
   * @param y          Count
   * @param n          Population size
   * @param lambda     Poisson mean
   * @param p          Detection probability
   *
   * @return Log probability
   */
  real poisbin_lpmf(array[] int y, int n, real lambda, real p) {
    if (max(y) > n) {
      return negative_infinity();
    }
    return poisson_log_lpmf(n | log(lambda)) + binomial_lpmf(y | n, p);
  }
  
  /**
   * Return log probability of N-mixture model for a site
   * 
   * @param count    Count in a site
   * @param max_n    Maximum site-abundance
   * @param lambda   Mean site-abundance
   * @param p        Per capita detection probability
   *
   * @return         Log probability
   */
   
  vector pb_logp(int[] count, int max_n,
                      real lambda, real p) {
                 
    int c_max = max(count);
    vector[max_n + 1] lp;

    for (k in 0:(c_max - 1))
      lp[k + 1] = negative_infinity();
    for (k in c_max:max_n) 
      lp[k + 1] = poisbin_lpmf(count| k, lambda, p);
    return lp;
  }
  
  real pb_lpmf(int[] count, int max_n,
                      real lambda, real p) {
                        
    vector[max_n + 1] lp;
    lp = pb_logp(count, max_n, lambda, p);
    
    return log_sum_exp(lp);
                        
  }
  
  real partial_sum(int[] site,
                   int start, int end,
                   int[, ] count,
                   int max_n, real lambda, real p) {
    real lp = 0;

    for (m in start:end)
      lp = lp + pb_lpmf(count[m] | max_n, lambda, p);
    return lp;
  }
}

data {
  int<lower = 0> R;               // Number of sites
  int<lower = 0> T;               // Number of replications
  int<lower = 0> K;               // Upper limit of the abundance
  //int<lower = 0> N[R];            // Abundance for each site
  int<lower = 0> y[R, T];         // Count for each site and replication
}

transformed data {
  int site[R] = rep_array(0, R);  // Dummy for site index
}

parameters {
  real<lower = 0> lambda;         // Mean site-abundance
  real<lower = 0, upper = 1> p;   // Per capita detection probability
}

model {
  int grainsize = 1;

  // Priors
  lambda ~ cauchy(0, 10);
  // A flat prior [0, 1] is implicitly used on p.
  
  // Partial Likelihood
  target += reduce_sum(partial_sum, site, grainsize, y, K, lambda, p);
}

generated quantities {
  vector[R] log_lik;                // Pointwise log likelihood
  real<lower=0,upper=1> psi;        // Mean site-occupancy
  real<lower=0,upper=1> P;          // Species detection probability 
  int totalN;                       // Total population size
  array[R, T] int y_new;            // New data points
  real fit = 0;                     // Test statistic actual data
  real fit_new = 0;                 // Test statistic new data
  real chat;                        // C-statistic (degree of overdisperion)
  int<lower=0, upper=1> pval = 0;   // Bayesian p-value for test statistic
  
  P = 1 - (1 - p)^poisson_rng(lambda);
  psi = 1 - exp(-lambda);
  
  {
    array[R] int N;                 // Abundance
    array[R] real eval;             // Expected values
    matrix[R, T] E;
    matrix[R, T] E_new;
    
    // Initialize counter, N, E and E_new
    N = rep_array(0, R);
    E = rep_matrix(0, R, T);
    E_new = rep_matrix(0, R, T);
    
    for (i in 1 : R) {
      
      // pointwise log likelihood contributions
      log_lik[i] = pb_lpmf(y[i,] | K, lambda, p);
          
      // latent N's and expected counts
      vector[K + 1] lp;
      lp = pb_logp(y[i,], K, lambda, p);
      N[i] = categorical_rng(softmax(lp)) - 1;
      eval[i] = p * N[i];
          
      for (t in 1 : T) {
        // Assess model fit using Chi-squared discrepancy
        // Compute fit statistic E for observed data
        E[i, t] = square(y[i, t] - eval[i]) / (eval[i] + 0.5);
        // Generate replicate data and
        // Compute fit statistic E_new for replicate data
        y_new[i, t] = binomial_rng(N[i], p);
        E_new[i, t] = square(y_new[i, t] - eval[i])/ (eval[i] + 0.5);
      }
      // Summed fit statistics
      fit = fit + sum(E[i,]);
      fit_new = fit_new + sum(E_new[i,]);
    }
    
    // Bayesian P-values
    chat = fit/ fit_new;
    if (fit_new > fit) pval = 1;

    // Total pop. size across all sampled sites
    totalN = sum(N[1 : R]);
  }
}

