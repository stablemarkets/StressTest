data {
  int<lower=0> n; 
  int<lower=0, upper=1> a[n]; 
  int<lower=0, upper=1> l[n]; 
  int<lower=0, upper=1> y[n]; 
  real xi1;
  real xi2;
}

parameters {
  
  real beta0;
  real beta1;
  real beta2;
  real u[n];
  
  real gam0;
  real gam1;
  
  real<lower=0, upper=1> theta1;
}

model {
  
  theta1 ~ beta(2,2);
  u ~ normal(0,1);
  beta1 ~ normal(0,1);
  
  for(i in 1:n){
    y[i] ~ bernoulli( inv_logit( beta0 + beta1*a[i] + beta2*l[i] + xi1 * u[i] ) );
    a[i] ~ bernoulli( inv_logit( gam0 + gam1*l[i] + xi2* u[i] ) );
    l[i] ~ bernoulli(theta1);
  }
  
}


generated quantities {
  
  int J = 500; // number of MC iterations
  real u_draw;
  real causal_effect = 0;
  
  for(j in 1:J)  {
    
    u_draw = normal_rng(0,1);
    causal_effect += (1.0/J)*( ( inv_logit(beta0 + beta1 *1 + beta2 * 1 + xi1 * u_draw) - 
                       inv_logit(beta0 + beta1*0 + beta2 * 1 + xi1*u_draw ) )*theta1 + 
                       ( inv_logit(beta0 + beta1 *1 + xi1 * u_draw) - 
                       inv_logit(beta0 + beta1*0 + xi1*u_draw ) )*(1-theta1) );
  }
}
