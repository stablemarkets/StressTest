#_______________________________________________________________________________
## Load Packaged _______________________________________________________________
#_______________________________________________________________________________

library(rstan)
library(latex2exp)
library(ggplot2)
library(reshape2) 

#_______________________________________________________________________________
## Simulate some data __________________________________________________________
#_______________________________________________________________________________

set.seed(1)

n = 300 ## 300 subjects
u = rnorm(n, 0, 1) ## unmeasured confounder, higher u = sicker
l = rbinom(n, 1, 0.5) ## measured confounder
a = rbinom(n, 1, plogis(1 + 1*l + 2*u ) ) ## sicker patients likely get A=1
y = rbinom(n,1 , plogis(1 + 1*l - 2*u + 0*a ) ) ## sicker patients less likely Y=1


## stan model implementing sensitivty method described in slides.
mod = stan_model("unmeasured_confounding.stan")

#_______________________________________________________________________________
## Estimate causal effect assuming no unmeasured confounding ___________________
#_______________________________________________________________________________

## frequentist standardization (over empirical distribution of l )
glm_res = glm(y ~ a + l, family=binomial(link='logit'))
summary(glm_res)
mean(l)* ( plogis( sum( glm_res$coefficients*c(1,1,1) ) ) - plogis( sum( glm_res$coefficients*c(1,0,1) ) ) ) +
  (1-mean(l))* ( plogis( sum( glm_res$coefficients*c(1,1,0) ) ) - plogis( sum( glm_res$coefficients*c(1,0,0) ) ) )

## Bayesian standardization (over empirical distribution of l. ideally should 
##   standardize over Bayesian bootstrap as in Oganisian & Roy 2021 but keeping
##   things simple for this course.)
stan_data = list(n=n, a=a, l=l, y=y, 
                 xi1 = 0, xi2 = 0 ## note xi1=xi2=0 implies to no unmeasured confounding
)
sampling_res = sampling(mod, seed=1,data = stan_data, chains = 1, iter = 2000, warmup = 1000)
summary(sampling_res, pars=c('causal_effect'))

#_______________________________________________________________________________
## Sensitivity analysis ________________________________________________________
#_______________________________________________________________________________
set.seed(11230)

## repeat the Bayesian analysis for different values of xi1, xi2
xi1v = seq(0, -3, length.out=10) # different xi1 values
xi2v = seq(0,  3, length.out=10) # different xi2 values

res = matrix(NA, nrow=10, ncol=10)

for( j in 1:10){
  for(k in 1:10){
    stan_data = list(n=n, a=a, l=l, y=y, xi1 = xi1v[j] , xi2 = xi2v[k] )
    
    sampling_res = sampling(mod, seed=1,data = stan_data, chains = 1, iter = 2000, warmup = 1000)
    
    draws = summary(sampling_res, pars=c('causal_effect'))
    res[j,k] = draws$summary[8] ## keep upper bound of 95% credible interval
  }
}

#_______________________________________________________________________________
## Plot Results ________________________________________________________________
#_______________________________________________________________________________
colnames(res) = xi2v
rownames(res) = xi1v
df <- melt(res)
colnames(df) = c("Y", "X", "value")

threshold = 0
# Create a new column to indicate above/below threshold
df$highlight = ifelse(df$value > threshold, "Above", "Below")

ggplot(df, aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  geom_text(data = subset(df, value > threshold), aes(label = round(value, 2)), color = "black") +
  scale_fill_gradient2(low = "skyblue", mid = "white", high = "darkred", midpoint = threshold) +
  theme_minimal() +
  labs(title = TeX("Upper limit of 95% CrI across priors. At $\\xi_1=\\xi_2=0$, -.31 [-.39,-.23]"), 
       x=TeX("$\\xi_1$"), y=TeX("$\\xi_2") ) +
  theme(
    axis.title.x = element_text(size = 16),  # x-axis label
    axis.title.y = element_text(size = 16),  # y-axis label
    axis.text.x = element_text(size = 14),   # x-axis tick labels
    axis.text.y = element_text(size = 14)    # y-axis tick labels
  )