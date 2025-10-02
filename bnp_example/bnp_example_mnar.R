library(rstan)
library(scales)

#####-------------------------------------------------------------------------##
#####------------------         Analysis Under MAR     -----------------------##
#####-------------------------------------------------------------------------##
n <- 200 ## sample size
set.seed(3)

# simulate some training data
l = seq(1,10*pi, length.out = n) + rnorm(n) # confounder
a = rbinom(n, 1, plogis( -2 + .1*l) )
y = rnorm( n = length(l), sin(.5*l) + 1.5*a, .07*l)

d = data.frame(y= as.numeric(scale(y)),
               a=a, 
               l = as.numeric(scale(l)))

d$delta = rbinom(n, 1, ifelse( -.5 < d$l & d$l < .5, .9, .01 ))

#plot(d$l,d$y, pch=20, col=ifelse(d$delta==1,'red', 'black'))

## load stan model
model = stan_model("bnp_example_mnar.stan")

n_m = sum(d$delta)
n_o = n - n_m

grid = seq(-1.75, 1.75, length.out=1000)

#####-------------------------------------------------------------------------##
#####------------------         Analysis Under MAR     -----------------------##
#####-------------------------------------------------------------------------##

## data to pass to Stan
stan_data = list(n=n, n_o = n_o, n_m=n_m,
                 y_o = d$y[d$delta==0], 
                 a_o = d$a[d$delta==0], 
                 l_o = d$l[d$delta==0], 
                 a_m = d$a[d$delta==1], 
                 l_m = d$l[d$delta==1],
                 
                 ## sensitivity parameters
                 xi3 = 0,
                 
                 ## for plotting regression line
                 l_test = grid,
                 N_test = length(grid),
                 
                 P = 10)

## run Stan model
#fit = sampling(model, seed=2,data = stan_data, iter = 2000, warmup=1000, chains = 1)
fit = sampling(model, seed=3,data = stan_data, iter = 1000, warmup=500, chains = 1)

## extract posterior draws of E[Y|X] for each X in test set.
post1 = extract(fit, pars= c('cond_mean_y_a1') )$cond_mean_y_a1
post0 = extract(fit, pars= c('cond_mean_y_a0') )$cond_mean_y_a0
post_ym = extract(fit, pars= c('y_m') )$y_m

#####-------------------------------------------------------------------------##
#####------------------        Analysis Under MNAR     -----------------------##
#####-------------------------------------------------------------------------##

## data to pass to Stan
stan_data = list(n=n, n_o = n_o, n_m=n_m,
                 y_o = d$y[d$delta==0], 
                 a_o = d$a[d$delta==0], 
                 l_o = d$l[d$delta==0], 
                 a_m = d$a[d$delta==1], 
                 l_m = d$l[d$delta==1],
                 
                 ## sensitivity parameters
                 xi3 = -1,
                 
                 ## for plotting regression line
                 l_test = grid,
                 N_test = length(grid),
                 
                 K = 10)

## run Stan model
#fit = sampling(model, seed=2,data = stan_data, iter = 2000, warmup=1000, chains = 1)
fit_mnar = sampling(model, seed=3,data = stan_data, iter = 1000, warmup=500, chains = 1)

## extract posterior draws of E[Y|X] for each X in test set.
post1_mnar = extract(fit_mnar, pars= c('cond_mean_y_a1') )$cond_mean_y_a1
post0_mnar = extract(fit_mnar, pars= c('cond_mean_y_a0') )$cond_mean_y_a0
post_ym_mnar = extract(fit_mnar, pars= c('y_m') )$y_m

#####-------------------------------------------------------------------------##
#####------------------        Plot Results      -----------------------------##
#####-------------------------------------------------------------------------##

png(filename = 'bnp_example.png', width = 700, height = 700)
par(mfrow=c(2,2))

plot(d$l[d$delta==0], d$y[d$delta==0], xlim=c(-2,2), pch=20,col=ifelse(d$a==1, 'red', 'darkblue'),
     xlab='L', ylab='Y',ylim=c(-10,4), 
     main='Posterior Draws of Missing Outcomes \n(MAR)')
for(j in 1:500){ 
  points(d$l[d$delta==1 & d$a==1], post_ym[0+j, d$a[d$delta==1]==1 ], 
         col= alpha('pink',alpha=.1), pch=20, cex=.8)
}

for(j in 1:500){ 
  points(d$l[d$delta==1 & d$a==0], post_ym[0+j, d$a[d$delta==1]==0 ], 
         col=alpha('lightblue',alpha=.1), pch=20, cex=.8)
}
legend('topleft', legend= c('Treated, A=1','Untreated, A=0'), pch=c(20,20),col=c('red','darkblue') )

points(d$l[d$delta==1], y= rep(-10.51, sum(d$delta)), pch=3,col=ifelse(d$a==1, 'red', 'darkblue'))
points(d$l[d$delta==0], d$y[d$delta==0], pch=20,col=ifelse(d$a==1, 'red', 'darkblue'))


plot(d$l[d$delta==0], d$y[d$delta==0], xlim=c(-2,2), pch=20,col=ifelse(d$a==1, 'red', 'darkblue'),
     xlab='L', ylab='Y', ylim=c(-10,4), main='Posterior Draws of Regression E[Y|A=a,L=l], a=0,1 \n(MAR)')
points(d$l[d$delta==1], y= rep(-10.51, sum(d$delta)), pch=3,col=ifelse(d$a==1, 'red', 'darkblue'))

for(j in 1:50){
  lines(grid, post1[100+j,], col=alpha('pink',alpha=.5), pch=20)
}
lines(grid, apply(post1,2,median), col='red', lwd=2) 

for(j in 1:50){
  lines(grid, post0[100+j,], col=alpha('lightblue',alpha=.5), pch=20)
}
lines(grid, apply(post0,2,median), col='darkblue', lwd=2) 

points(d$l[d$delta==0], d$y[d$delta==0], pch=20,col=ifelse(d$a==1, 'red', 'darkblue'))


plot(d$l[d$delta==0], d$y[d$delta==0], xlim=c(-2,2), pch=20,col=ifelse(d$a==1, 'red', 'darkblue'),
     xlab='L', ylab='Y',ylim=c(-10,4), 
     main='Posterior Draws of Missing Outcomes \n(MNAR)')
for(j in 1:500){ 
  points(d$l[d$delta==1 & d$a==1], post_ym_mnar[0+j, d$a[d$delta==1]==1 ], 
         col= alpha('pink',alpha=.1), pch=20, cex=.8)
}

for(j in 1:500){ 
  points(d$l[d$delta==1 & d$a==0], post_ym_mnar[0+j, d$a[d$delta==1]==0 ], 
         col=alpha('lightblue',alpha=.1), pch=20, cex=.8)
}
legend('topleft', legend= c('Treated, A=1','Untreated, A=0'), pch=c(20,20),col=c('red','darkblue') )

points(d$l[d$delta==1], y= rep(-10.51, sum(d$delta)), pch=3,col=ifelse(d$a==1, 'red', 'darkblue'))
points(d$l[d$delta==0], d$y[d$delta==0], pch=20,col=ifelse(d$a==1, 'red', 'darkblue'))


plot(d$l[d$delta==0], d$y[d$delta==0], xlim=c(-2,2), pch=20,col=ifelse(d$a==1, 'red', 'darkblue'),
     xlab='L', ylab='Y', ylim=c(-10,4), main='Posterior Draws of Regression E[Y|A=a,L=l], a=0,1 \n(MNAR)')
points(d$l[d$delta==1], y= rep(-10.51, sum(d$delta)), pch=3,col=ifelse(d$a==1, 'red', 'darkblue'))

for(j in 1:50){
  lines(grid, post1_mnar[100+j,], col=alpha('pink',alpha=.5), pch=20)
}
lines(grid, apply(post1_mnar,2,median), col='red', lwd=2) 

for(j in 1:50){
  lines(grid, post0_mnar[100+j,], col=alpha('lightblue',alpha=.5), pch=20)
}
lines(grid, apply(post0_mnar,2,median), col='darkblue', lwd=2) 

points(d$l[d$delta==0], d$y[d$delta==0], pch=20,col=ifelse(d$a==1, 'red', 'darkblue'))

dev.off()