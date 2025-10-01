data {
  
  // data
  int<lower=0> n; 
  int<lower=0, upper=1> l[n];
  int<lower=0, upper=1> a_tilde[n];
  int<lower=0, upper=1> y[n];
  
  // sensitivity point mass priors
  real<lower=0, upper=1> xi1;
  real<lower=0, upper=1> xi2;
}

parameters {
  vector[3] eta; //outcome model parameters
  vector[2] gamma; // treatment model parameters
  real<lower=0, upper=1> theta; // covariate model parameter
}


model {
  
  // priors contributions
  eta ~ normal(0,3);
  gamma ~ normal(0,3);
  theta ~ beta(1,1);
  
  // covariate model likelihood contribution
  l ~ bernoulli(theta);
  
  // mixture term likelihood contribution
  for(i in 1:n){
    target += log_sum_exp( 
      // term for a=1
      bernoulli_lpmf(a_tilde[i] | xi1) 
      + bernoulli_lpmf(y[i] | inv_logit( eta[1] + eta[2]*l[i] + eta[3]))
      + bernoulli_lpmf(1 | inv_logit( gamma[1] + gamma[2]*l[i]) ),   
      // term for a=0
      bernoulli_lpmf(a_tilde[i] | xi2)
      + bernoulli_lpmf(y[i] | inv_logit( eta[1] + eta[2]*l[i]) )
      + bernoulli_lpmf(0 | inv_logit( gamma[1] + gamma[2]*l[i]) ) 
    );
  }
  
}

generated quantities {
  
  // compute average treatment effect using standardization formula
  real ATE;
  
  ATE = theta*(inv_logit(eta[1] + eta[2] + eta[3]) - inv_logit(eta[1] + eta[2])) +
    (1-theta)*(inv_logit( eta[1] + eta[3]) - inv_logit( eta[1] ) );
  
}
