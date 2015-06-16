data{
  int degree;         #polynomial degree. See Efron 2006 Size, ....
  int K;              #Number of bins. See Efron 2006 Size, ....
  int <lower=0> y[K]; #Number of total counts in each bin
  real z[K,degree];   #z-score^i for bins i in 1:degree
}
parameters{
  real beta[degree];           #Regression coefs. See Efron 2006 Size, ....
  real<lower=0> sigma_sq_beta;
  real mu;                     #The mean term
  real epsilon[K];             #Added over-dispersion parameter.
  real<lower=0> sigma_sq_epsilon;
}
transformed parameters{
  real<lower=0> sigma_beta;
  real<lower=0> sigma_epsilon;

  sigma_beta <- sqrt(sigma_sq_beta);
  sigma_epsilon <- sqrt(sigma_sq_epsilon);
}
model{
  sigma_sq_beta ~ inv_gamma(0.001,0.001);
  mu ~ normal(0,sigma_beta);
  beta ~ normal(0,sigma_beta);
  epsilon ~ normal(0,sigma_epsilon);
  
  for(i in 1:K) {
    y[i] ~ poisson_log(mu + dot_product(beta,z[i]) + epsilon[i]);    
  }  
}
generated quantities{
  vector[K] lambda;
  for(i in 1:K) {
  	lambda[i] <- exp(mu + dot_product(beta,z[i]) + epsilon[i] );
  }
}