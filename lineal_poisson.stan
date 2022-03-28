
data {
  int<lower=0> n;
  int count[n];
  vector[n] t;
}


parameters {
  real N0;
  real beta;
}

model {

  count ~ poisson(N0 - beta*t);
  
  // Priors
  
  beta ~ normal(.05, 1);
  N0 ~ normal(1e8, 100);

}







