library(rstan)
library(latex2exp)

###--- Simulate Some Data ---###
set.seed(1)
n = 1000
l = rbinom(n, 1, .5)
a = rbinom(n, 1,plogis(1 - 2*l))
y = rbinom(n, 1, plogis(1 - 1*l + 1.5*a ))
delta = rbinom(n, 1, plogis(2 + 1*a - 0*y + 1*a*y ) )

# among a=1

## delta = rbinom(n, 1, plogis(2 + 1 + 1*y ) )
## patients with remission more likely to be missing. 
## i.e. observed treated subjects more likely to have (y=0)
## might make it seem like treatment a=1 doesn't lower remission. 

true_ate = (plogis(1 - 1*1 + 1.5*1 ) - plogis(1 - 1*1 + 1.5*0 ))*.5 + 
           (plogis(1 - 1*0 + 1.5*1 ) - plogis(1 - 1*0 + 1.5*0 ))*.5

true_ate

### --- Create Stan Data --- ###
y_o = y[delta==0]

a_o = a[delta==0]
a_m = a[delta==1]

l_o = l[delta==0]
l_m = l[delta==1]

n_o = sum(1-delta)
n_m = n - n_o

stan_data = list(n=n, n_o=n_o, n_m = n_m,
                 l_o=l_o, l_m=l_m,
                 a_o=a_o, a_m=a_m,
                 y_o=y_o,
                 xi2=0, xi3=0)
 
## stan model implementing sensitivty
mod = stan_model("mnar_outcome.stan")

res = sampling(mod, seed=31, data = stan_data,
               chains = 1, iter = 2000, warmup = 1000)

draws = summary(res, pars=c('ATE','eta','xi','theta') )$summary[1,c(1,4,8)]
draws

## repeat the Bayesian analysis for different values of xi1, xi2
xi3v = seq(0, 3, by=.25) # different xi2 values
nvals = length(xi3v)

res = matrix(NA, nrow=nvals, ncol=3)

for(k in 1:nvals){
  
  stan_data = list(n=n, n_o=n_o, n_m = n_m, 
                   l_o=l_o, l_m=l_m,
                   a_o=a_o, a_m=a_m, 
                   y_o=y_o, 
                   xi2 = 0 , xi3 = xi3v[k])
  
  sampling_res = sampling(mod, seed=31, data = stan_data, 
                          chains = 1, iter = 1000, warmup = 500)
  
  draws = summary(sampling_res, pars=c('ATE'))
  res[k,] = draws$summary[c(1,4,8)] ## keep inference
}


png('mnar.png', width = 500, height = 500)

par(mar = c(5, 6, 5, 2))  

plot(xi3v, res[,1], ylim=c(-1,1), pch=20,
     cex.lab=1.5, cex.axis=1.5, lwd=1.5, cex=1.5, cex.main=1.5,
     main = 'Posterior Inference for ATE \nUnder Different MAR Violations',
     ylab = TeX("Posterior Mean & Credible Intervals for ATE"), 
     xlab = TeX("$\\xi_3$"))
segments(xi3v, res[,2], xi3v, res[,3], lwd=1.5)
abline(h=0, col='red',lty=2, lwd=1.5)

dev.off()


