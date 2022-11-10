// Zero-Inflated Poisson-Binomial mixture model - multithreading

functions{
  
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
   * Return log probability of a zero inflated N-mixture model for a site
   * 
   * @param count        Count in a site
   * @param max_n        Maximum site-abundance
   * @param log lambda   Mean site-abundance (log scale)
   * @param p            Mean per capita detection probability
   *
   * @return         Log probability
   */
   
  real zipb_lpmf(int[] count, int max_n,
               real lambda, real psi, real p) {
                 
    int c_max = max(count);
    vector[max_n + 1] lp;
    
    if(c_max){  // If non-zero max count
    
      // Probability that the site is occupied
      lp = rep_vector(bernoulli_lpmf(1 | psi), max_n + 1);
      
      // Probability of k individuals, when 'count' are counted
      for (k in 0:(c_max - 1)){
        lp[k + 1] = lp[k + 1] + negative_infinity();
      }
      for (k in c_max:max_n){
        lp[k + 1] = lp[k + 1] + poisbin_lpmf(count| k, lambda, p);
        
      }
    } else {  // If max count is zero
      
      vector[max_n + 1] lp2[2];
      // Probability that the site is unoccupied
      lp2[1,] = rep_vector(bernoulli_lpmf(0 | psi), max_n + 1);
      // Probability that the site is occupied
      lp2[2,] = rep_vector(bernoulli_lpmf(1 | psi), max_n + 1);
      
      // Probability of k individuals, when 'count' are counted
      for (k in 0:(c_max - 1)){
        lp2[2, k + 1] = lp2[2, k + 1] + negative_infinity();
      }
      
      for (k in c_max:max_n){
        lp2[2, k + 1] = lp2[2, k + 1] + poisbin_lpmf(count| k, lambda, p);
        
      }
      
      // Sum the log probabilities for 'unoccupied' and 'occupied' states
      for (k in 0:max_n){
        lp[k + 1] = log_sum_exp(lp2[1:2, k + 1]);
      }
    }
    return log_sum_exp(lp);
  }

  real partial_sum(int[] site,
                   int start, int end,
                   int[,] count, int max_n, 
                   real lambda, real psi, real p) {
    real lp = 0;

    for (m in start:end){
      lp = lp + zipb_lpmf(count[m] | max_n, lambda, psi, p);
    }
    return lp;
  }
}


data {
  
  int<lower = 0> R;               // Number of sites
  int<lower = 0> T;               // Number of surveys (replicates)
  int<lower = 0> K;               // Upper limit of the abundance
  int<lower = 0> y[R, T];         // Count for each site and replication

}

transformed data {
  int site[R] = rep_array(0, R);  // Dummy for site index
}

parameters {
  real<lower=0> lambda;           // Mean site-abundance
  real<lower=0, upper=1> psi;     // Mean occupancy probability
  real<lower=0, upper=1> p;       // Mean pc detection probability
 
}

model {
  int grainsize = 1;

  // Priors
  lambda ~ cauchy(0, 10);
  // Flat priors [0, 1] are implicitly used on p and psi.
  
  // Partial Likelihood
  target += reduce_sum(partial_sum, site, grainsize, y, K, lambda, psi, p);
}

generated quantities {
  vector[R] log_lik;                // Pointwise log likelihood
  real<lower=0,upper=1> P;          // Species detection probability 
  int totalN;                       // Total population size
  array[R, T] int y_new;
  real fit = 0;                     // Test statistic actual data
  real fit_new = 0;                 // Test statistic new data
  real chat;                        // C-statistic (degree of overdisperion)
  int<lower=0, upper=1> pval = 0;   // Bayesian p-value for test statistic
  
  P = 1 - (1 - p)^poisson_rng(lambda);
  
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
      log_lik[i] = zipb_lpmf(y[i,] | K, lambda, psi, p);
          
      // latent N's and expected counts
      vector[K + 1] lp;
      
      for (n in 0 : K) {
        lp[n + 1] = poisbin_lpmf(y[i,] | K, lambda, p);
      }
      
      N[i] = categorical_rng(softmax(lp)) - 1;
          
      real log_p_unobs; // Log of prob. site is suitable
      if (max(y[i,]) == 0) {
        // Unobserved
        log_p_unobs = log(psi) + binomial_lpmf(0 | N[i], p);
        if (bernoulli_rng(exp(log_p_unobs)) == 0) {
          // Site is not suitable
          N[i] = 0;
        }
      }
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
      fit = fit + sum(E[i]);
      fit_new = fit_new + sum(E_new[i]);
    }
    // Bayesian P-values
    chat = fit/ fit_new;
    if (fit_new > fit) pval = 1;

    // Total pop. size across all sampled sites
    totalN = sum(N[1 : R]);
  }
}
