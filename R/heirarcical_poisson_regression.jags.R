#Default input values for heirarcical overdispersed poisson 
#regression:
count.data <- list ("K", "y", "z", "degree")

count.parameters <- c ("offset","beta","sig.beta","lambda")

count.inits <- function (){
  list(beta=rnorm(degree), offset=rnorm(1), sig.beta=runif(1))
}

#--------------------------------------------------
#JAGS or WinBUGS heirarcical poisson 
#regression to be used with Efron's 2-groups 
#empirical Bayes methodology for local false discovery
#rates
#--------------------------------------------------
model.func.hpr.BUGS<-function() {
  
  #construct a temp vec to sum. Avoids use of inprod in the next loop.
  for(i in 1:K) {  
    for (j in 1:degree) {
      zb.tmp[i,j] <- beta[j]*z[i,j] 
    } 
  }
  
  for(i in 1:K) {
    y[i] ~ dpois(lambda[i])
    #log(lambda[i]) <- offset + inprod(beta[],z[i,])
    #log(lambda[i]) <- offset + sum(zb.tmp[i,])
    log(lambda[i]) <- offset + sum(zb.tmp[i,])
  }
  offset ~ dnorm(0,tau.beta)
  
  
  for(j in 1:degree) {
    beta[j] ~ dnorm(0,tau.beta)
  }  
  tau.beta <- pow(sig.beta,-2)
  sig.beta ~ dunif(0,100)
}