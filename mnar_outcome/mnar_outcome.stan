
data {
  int<lower=0> n;
  int<lower=0> n_o;
  int<lower=0> n_m;
  
  int y_o[n_o];
  int a_o[n_o];
  int l_o[n_o];
  
  int a_m[n_m];
  int l_m[n_m];
  
  real xi2;
  real xi3;
}

parameters {
  real<lower=0, upper=1> theta;
  
  real eta[3];
  real xi[2];
}


model{
  
  eta ~ normal(0,3);
  xi ~ normal(0,3);
  
  
  for(i in 1:n_o){
    y_o[i] ~ bernoulli(inv_logit( eta[1] + eta[2]*l_o[i] + eta[3]*a_o[i]));
    l_o[i] ~ bernoulli(theta);
    target += bernoulli_lpmf(0 | inv_logit(xi[1] + xi[2]*a_o[i] + xi2*y_o[i] + xi3*y_o[i]*a_o[i] ) );
  }

  for(i in 1:n_m){
    
  target += log_sum_exp( 
    
    bernoulli_lpmf(1 | inv_logit(xi[1] + xi[2]*a_m[i] + xi2*1.0 + xi3*1.0*a_m[i])) 
    + bernoulli_lpmf(1 | inv_logit( eta[1] + eta[2]*l_m[i] + eta[3]*a_m[i])),
                
    bernoulli_lpmf(1 | inv_logit(xi[1] + xi[2]*a_m[i] + xi2*0 + xi3*0*a_m[i])) 
    + bernoulli_lpmf(0| inv_logit( eta[1] + eta[2]*l_m[i] + eta[3]*a_m[i]))
                
  );
            
  target += bernoulli_lpmf(l_m[i] | theta);
    
  }
  
}

generated quantities {
  
  real ATE;
  
    ATE = theta*(inv_logit(eta[1] + eta[2] + eta[3]) - inv_logit(eta[1] + eta[2])) +
        (1-theta)*(inv_logit( eta[1] + eta[3]) - inv_logit( eta[1] ) );
  
}

