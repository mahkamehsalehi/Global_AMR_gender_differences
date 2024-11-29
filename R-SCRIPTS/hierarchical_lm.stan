data{
  
  int<lower=0> N; // Total sample size
  int<lower=0> G; // Number of predictors
  
  vector[N] y; // response
  matrix[N, G] x; // predictors
  
  int n_groups;
  int<lower=0, upper=n_groups> ind[N]; // Indicator denoting group
  
}

parameters {
  
  vector[n_groups] intercept;
  matrix[G, n_groups] coef;
  
  real<lower=0> sigma;
  
  real mu_coef;
  real<lower=0> sigma_coef;
}

model {
  
  for(i in 1:N) {
    int my_group = ind[i];
    
    y[i] ~ normal( intercept[my_group] + dot_product(coef[, my_group], x[i, ]), sigma);
    
  }
  
  // Priors
  intercept ~ normal(0, 1);
  
  for(i in 1:G) {
    for(j in 1:n_groups) {
      coef[i, j] ~ normal(mu_coef, sigma_coef);
    }
  }
  
  mu_coef ~ normal(0, 1);
  sigma_coef ~ gamma(2, 1);
  
  sigma ~ gamma(2, 1);
}