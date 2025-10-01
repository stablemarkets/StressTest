library(rstan)
library(latex2exp)

setwd("~/google_drive/Research/epi_sensitivity_analysis/manuscript_code/misclassification_example/")

###--- Simulate Some Data ---###
set.seed(10)

## set actual misclassification rates
xi1 = .90 ## 80% sensitivity
xi2 = .1 ## implies 1-xi2=80% specificity

n = 500
l = rbinom(n,1, .5 ) ## simulate confounder
a = rbinom(n, 1, plogis( 0 + 1.1*l ) ) ## simulate true treatment
a_tilde = rbinom(n, 1, a*xi1 + xi2*(1-a) ) ## simulate error-prone treatment
y = rbinom(n, 1, plogis(1 - 1*l + 1.1*a) ) ## simulate outcome

## compute true ATE as a rough sanity-check/benchmark
true_ATE = .5*(plogis(1 - 1*1 + 1.1)-plogis(1 - 1*1 + 1.1*0)) + 
           .5*(plogis(1 - 1*0 + 1.1)-plogis(1 - 1*0 + 1.1*0))


### --- Load Stan model for exposure misclassification --- ###
misclass_mod = stan_model("exposure_misclassification.stan") 


## --- run under different sensitivity value parameters
xi1_vals = c(.99, .95, .90, .85, .80)
xi2_vals = c(.001, .01, .1, .15, .20)
n1vals = length(xi1_vals)
n2vals = length(xi2_vals)

resmat = matrix(nrow=n1vals*n2vals, ncol=5)

counter = 1
for(u in 1:n1vals){
  for(v in 1:n2vals){
    
    stan_data = list(n=n, y= y, a_tilde = a_tilde, l=l, xi1=xi1_vals[u], xi2=xi2_vals[v] )
    
    res = sampling(misclass_mod, 
                   seed=31, 
                   data = stan_data, 
                   chains = 1, iter = 2000, warmup = 1000)
    
    resmat[counter,1:2] = c(xi1_vals[u], xi2_vals[v])
    resmat[counter,3:5] = summary(res, pars=c('ATE') )$summary[,c("mean","2.5%","97.5%")]
    
    counter = counter+1
  }
}


dmat = data.frame(resmat)

png('misclass.png', width = 500, height = 500)
par(mar = c(5, 6, 5, 2))  
xi2v = .1

plot(dmat$X1[dmat$X2==xi2v], dmat$X3[dmat$X2==xi2v], pch=20, type='o', ylim=c(.1,.6), xlim=c(.77, 1.03), 
     xlab=latex2exp("sensitivity, $\\xi_1$"),
     ylab=latex2exp("Posterior Point & Intervals for ATE, $\\Psi$"), 
     main = 'Posterior ATE Estimates Under Different \nExposure Misclassification Settings',
     cex.lab=1.5, cex.axis=1.5, lwd=1.5, cex=1.5, col='orange')
segments( dmat$X1[dmat$X2==xi2v], dmat$X4[dmat$X2==xi2v], 
          dmat$X1[dmat$X2==xi2v],dmat$X5[dmat$X2==xi2v], 
          lwd=1.5, cex=1.5, col='orange' )

xi2v = .15
points(dmat$X1[dmat$X2==xi2v]-.005, dmat$X3[dmat$X2==xi2v],
       pch=20, type='o', col='purple', cex=1.5)
segments( dmat$X1[dmat$X2==xi2v]-.005, dmat$X4[dmat$X2==xi2v], 
          dmat$X1[dmat$X2==xi2v]-.005,dmat$X5[dmat$X2==xi2v], 
          lwd=1.5, cex=1.5, col='purple' )

xi2v = .2
points(dmat$X1[dmat$X2==xi2v]-0.01, dmat$X3[dmat$X2==xi2v], 
       pch=20, type='o', col='black', cex=1.5)
segments( dmat$X1[dmat$X2==xi2v]-0.01, dmat$X4[dmat$X2==xi2v], 
          dmat$X1[dmat$X2==xi2v]-0.01,dmat$X5[dmat$X2==xi2v], 
          lwd=1.5, cex=1.5, col='black' )

xi2v = .01
points(dmat$X1[dmat$X2==xi2v]+.005, dmat$X3[dmat$X2==xi2v],
       pch=20, type='o', col='blue', cex=1.5)
segments( dmat$X1[dmat$X2==xi2v]+.005, dmat$X4[dmat$X2==xi2v], 
          dmat$X1[dmat$X2==xi2v]+.005,dmat$X5[dmat$X2==xi2v], 
          lwd=1.5, cex=1.5, col='blue' )


xi2v = .001
points(dmat$X1[dmat$X2==xi2v]+.01, dmat$X3[dmat$X2==xi2v], 
       pch=20, type='o', col='green', cex=1.5)
segments( dmat$X1[dmat$X2==xi2v]+.01, dmat$X4[dmat$X2==xi2v], 
          dmat$X1[dmat$X2==xi2v]+.01,dmat$X5[dmat$X2==xi2v], 
          lwd=1.5, cex=1.5, col='green' )

abline(h=dmat$X3[dmat$X2==.001][1], lty=2, lwd=1.5, col='red')


legend("topright",
       legend = latex2exp::TeX(c(
         "$\\xi_2 = 0.20$",
         "$\\xi_2 = 0.15$",
         "$\\xi_2 = 0.10$",
         "$\\xi_2 = 0.01$",
         "$\\xi_2 = 0.001$"
       )),
       col = c("black", "purple", "orange", "blue", "green"),
       pch = 20, pt.cex = 1.5,
       lty = 1, lwd = 1.5,
       bty = "n",   # no legend box
       cex = 1.2)   # legend text size
dev.off()
