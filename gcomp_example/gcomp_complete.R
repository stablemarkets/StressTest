library(rstan)

###--- Simulate Some Data ---###
set.seed(3)
n = 1000
l = rbinom(n, 1, .5)
a = rbinom(n, 1,plogis(1 - 2*l))
y = rbinom(n, 1, plogis(1 - 1*l + 1.5*a ))

###--- Compute True ATE ---###

true_ate = (plogis(1 - 1*1 + 1.5*1 ) - plogis(1 - 1*1 + 1.5*0 ))*.5 + 
           (plogis(1 - 1*0 + 1.5*1 ) - plogis(1 - 1*0 + 1.5*0 ))*.5

true_ate

### --- Create Stan Data --- ###
stan_data = list(n=n, y=y, a=a, l=l)
 
### --- Load the Stan Model --- ###
mod = stan_model("gcomp_complete.stan")

## Obtain 1000 posterior draws (after 1000 warm-up draws)
res = sampling(mod, seed=31, data = stan_data,
               chains = 1, iter = 2000, warmup = 1000)

## extract and print posterior summary of ATE
draws = summary(res, pars=c('ATE') )$summary[1,c(1,4,8)]
