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
  vector[3] eta;
  vector[2] gamma;
  real<lower=0, upper=1> theta;
}

transformed parameters{
  
  real<lower=0, upper=1> PY_a1[n];
  real<lower=0, upper=1> PY_a0[n];
  real<lower=0, upper=1> Pa1[n];

  for(i in 1:n){
    PY_a1[i] = inv_logit( eta[1] + eta[2]*l[i] + eta[3]);
    PY_a0[i] = inv_logit( eta[1] + eta[2]*l[i]);
    Pa1[i] = inv_logit( gamma[1] + gamma[2]*l[i]) ;
  }
  
  
}

model {
  eta ~ normal(0,3);
  gamma ~ normal(0,3);
  theta ~ beta(1,1);

  for(i in 1:n){
    target += log_sum_exp( 
                  bernoulli_lpmf(a_tilde[i] | xi1) 
                + bernoulli_lpmf(y[i] | PY_a1[i])
                + bernoulli_lpmf(1 | Pa1[i]),   
      
                  bernoulli_lpmf(a_tilde[i] | xi2)
                + bernoulli_lpmf(y[i] | PY_a0[i])
                + bernoulli_lpmf(0 | Pa1[i]) 
              );
              
      target += bernoulli_lpmf(l[i] | theta);
  }
  
}

generated quantities {
  
  real ATE;
  
  ATE = theta*(inv_logit(eta[1] + eta[2] + eta[3]) - inv_logit(eta[1] + eta[2])) +
        (1-theta)*(inv_logit( eta[1] + eta[3]) - inv_logit( eta[1] ) );
    
}
