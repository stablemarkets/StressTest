data{
    int<lower=0> K; 
    int<lower=0> n; 
    int<lower=0> n_o; 
    int<lower=0> n_m; 
    
    // grid of values for computing regression line
    int<lower=0> N_test;
    real l_test[N_test];
    
    // observed data for those with/without missing outcomes
    int a_o[n_o];
    real l_o[n_o];
    real y_o[n_o];
    int a_m[n_m];
    real l_m[n_m];

    // specified value of sensitivity parameter
    real xi3;
}

parameters {
  real<lower=0> ax;
  real<lower=0> bx;
  real cx;
  real<lower=0> dx;
  real<lower=0> a;
  real<lower=0> b;
  
  real b0m;
  real b1m;
  real b2m;

  real g0m;
  real g1m;
  
  real<lower=0> sigma0;
  real<lower=0> sigma1;
  real<lower=0> sigma2;

  real<lower=0> g_sigma0;
  real<lower=0> g_sigma1;
  
  vector[K] mu_x;
  vector<lower=0>[K] phi_x;
  
  vector[K] beta0;
  vector[K] beta1;
  vector[K] beta2;
  
  vector[K] gam0;
  vector[K] gam1;
  
  vector<lower=0>[K] tau_mix;
  real <lower=0,upper=1> q[K];
  real <lower=0> alpha;
  
  real y_m[n_m];
  
  real xi1;
  real xi2;
}

transformed parameters{
   simplex[K] w;

   w[1] = q[1];
   for(j in 2:(K-1) ){
     w[j] = q[j]*(1-q[j-1])*w[j-1]/q[j-1];
   }
   w[K] = 1 - sum(w[1:(K-1)]);
   
}

model {
  alpha ~ gamma(1,1);
  
  // runs a bit faster and better results with these hyperpriors
  a ~ normal(0,1);
  b ~ normal(0,1);
  
  ax ~ normal(0,1);
  bx ~ normal(0,1);
  
  b0m ~ normal(0,1);
  b1m ~ normal(0,1);
  b2m ~ normal(0,1);
  
  cx ~ normal(0,1);
  dx ~ normal(0,1);
  
  sigma2 ~ normal(0,1);
  sigma1 ~ normal(0,1);
  sigma0 ~ normal(0,1);
  
  g_sigma1 ~ normal(0,1);
  g_sigma0 ~ normal(0,1);
  
  xi1 ~ normal(0,1);
  xi2 ~ normal(0,1);

  for(j in 1:K){
    q[j] ~ beta(1, alpha);
    
    mu_x[j] ~ normal(cx,dx);
    phi_x[j] ~ normal(ax,bx);
    
    beta2[j] ~ normal( b2m, sigma2 );
    beta1[j] ~ normal( b1m, sigma1 );
    beta0[j] ~ normal( b0m, sigma0 );
    
    gam1[j] ~ normal( g1m, g_sigma1 );
    gam0[j] ~ normal( g0m, g_sigma0 );
    
    tau_mix[j] ~ inv_gamma ( a, b );
  }
  
  for(i in 1:n_o){
    vector[K] lps = log(w);
    for(k in 1:K){
      lps[k] += normal_lpdf(y_o[i] | beta0[k] + beta1[k]*l_o[i] + beta2[k]*a_o[i], tau_mix[k]);
      lps[k] += bernoulli_lpmf(a_o[i] | inv_logit(gam0[k] + gam1[k]*l_o[i]) );
      lps[k] += normal_lpdf(l_o[i] | mu_x[k], phi_x[k]);
    }
    target += log_sum_exp(lps);
    target += bernoulli_lpmf(0 | inv_logit(xi1 + xi2*a_o[i] + xi3*y_o[i]*a_o[i] ));
  }

  for(i in 1:n_m){
    vector[K] lps = log(w);
    for(k in 1:K){
      lps[k] += normal_lpdf(y_m[i] | beta0[k] + beta1[k]*l_m[i] + beta2[k]*a_m[i], tau_mix[k]);
      lps[k] += bernoulli_lpmf(a_m[i] | inv_logit(gam0[k] + gam1[k]*l_m[i] ) );
      lps[k] += normal_lpdf(l_m[i] | mu_x[k], phi_x[k]);
    }
    target += log_sum_exp(lps);
    target += bernoulli_lpmf(1 | inv_logit(xi1 + xi2*a_m[i] + xi3*y_m[i]*a_m[i] ));
  }


}

generated quantities{
  
  real cond_mean_y_a1[N_test];
  real cond_mean_y_a0[N_test];
  
  // E[Y | A=1, L=l]
  for( i in 1:N_test){
    vector[K] wght;
    vector[K] cond_mod;
    for(k in 1:K){
      wght[k] = w[k]*exp(normal_lpdf(l_test[i] | mu_x[k], phi_x[k] ) + 
                         bernoulli_lpmf(1 | inv_logit(gam0[k] + gam1[k]*l_test[i] ) ) );
      cond_mod[k] = (beta0[k] + beta1[k]*l_test[i] + beta2[k])*wght[k];
    }
    cond_mean_y_a1[i] = sum(cond_mod) / sum(wght) ;
  }
  
  // E[Y | A=0, L=l]
  for( i in 1:N_test){
    vector[K] wght;
    vector[K] cond_mod;
    for(k in 1:K){
      wght[k] = w[k]*exp(normal_lpdf(l_test[i] | mu_x[k], phi_x[k] ) + 
                         bernoulli_lpmf(0 | inv_logit(gam0[k] + gam1[k]*l_test[i] ) ) );
      cond_mod[k] = (beta0[k] + beta1[k]*l_test[i])*wght[k];
    }
    cond_mean_y_a0[i] = sum(cond_mod) / sum(wght) ;
  }
  
  
}


