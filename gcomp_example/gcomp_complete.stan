
data {
  // declare known quantities (data)
  int<lower=0> n; int y[n]; int a[n]; int l[n];
}

parameters {
  // declare unknown quantities (parameters)
  real<lower=0, upper=1> theta;
  real eta[3];
}


model{
  // specify prior contributions to log unnormalized posterior
  theta ~ beta(1,1);
  eta ~ normal(0,3);

  // specify data model contributions to log unnormalized posterior
  for(i in 1:n){
    y[i] ~ bernoulli(inv_logit( eta[1] + eta[2]*l[i] + eta[3]*a[i]));
    l[i] ~ bernoulli(theta);
  }
}

generated quantities {
  // perform standardization with each 
  // posterior draw of parameter to get ATE
  real ATE;
  
  ATE = theta*( inv_logit(eta[1]+eta[2]+eta[3]) 
                - inv_logit(eta[1]+eta[2])) +
    (1-theta)*( inv_logit(eta[1]+eta[3]) 
                - inv_logit(eta[1]) );
  
}

