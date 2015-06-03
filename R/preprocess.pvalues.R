preprocess.pvalues <- function(p.values, num.bins, degree, sampler) {
  
  z.values <- qnorm(p.values)
  z.info <- preprocess.z(z.values, num.bins)
  
  y <- z.info$y
  K <- z.info$K
  x <- z.info$x  
  z <- ns(x, df = degree)[,]
  
  #Global variables needed for JAGS runs
  if(sampler == "jags") {
    #Default input values for heirarcical overdispersed poisson regression:
    count.data <- list ("K", "y", "z", "degree")
    
    count.parameters <- c ("offset","beta","sig.beta","lambda","eps")
    
    count.inits <- function (){
      list(beta=rnorm(degree), offset=rnorm(1), sig.beta=runif(1), eps=rnorm(K))
    }    
  }
  
 if(sampler == "stan") {
   
 } 
  
}