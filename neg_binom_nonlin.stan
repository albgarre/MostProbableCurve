

data {
  int<lower=0> n;
  int count[n];
  vector[n] t;
}
parameters {
  real logN0;
  real logD;
  real<lower=1e-8> theta;
  real logp;
}
transformed parameters {

  real<lower=0> D;    
  real<lower=0> N0;
  real<lower=0> p;
  D = 10^logD;
  N0 = 10^logN0;
  p = 10^logp;
}
model {
  
  // Likelihood
  vector[n] mu;
  for (i in 1:n)
    mu[i] = N0 * 10 ^ ( - ((t[i] / D)^p));

  count ~ neg_binomial_2(mu, theta);
  
  // Priors
  
  logD ~ normal(1, 1);
  logN0 ~ normal(8, 1);
  theta ~ gamma(1, 1);
  logp ~ normal(0, 1);

}




