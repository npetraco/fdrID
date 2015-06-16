#--------------------------------------------
#' @title sampler.prep
#' @description XXXX
#' 
#' @details XXXX
#'
#' @param XXXX
#' 
#' @return XXXX
#' 
#' @references XXXX
#'
#' @examples
#' XXXX
#--------------------------------------------
sampler.prep <- function(p.values, num.bins=120, degree=7, interceptQ=FALSE, overdispersionQ=FALSE, sampler=c("jags","stan")) {
  
  z.values <- qnorm(p.values)
  z.info <- preprocess.z(z.values, num.bins)
  
  y <- z.info$counts
  K <- z.info$num.bins
  x <- z.info$bin.midpoints
  z <- ns(x, df = degree)[,]
  
  #Global variables needed for JAGS runs
  if(sampler == "jags") {
    
    #Common input values for heirarcical poisson regression:
    count.data <- list (K=K, y=y, z=z, degree=degree)
    
    if(interceptQ==TRUE & overdispersionQ==TRUE){ # Model with intercept and overdispersion terms
      
      #Get model file path:
      path.to.model.file <- system.file("jags","hierarchical.overdispersed.poission.regression.with.intercept.bug", package="fdrID")
      
      #Parameters for the model:
      count.parameters <- c ("offset","beta","sig.beta","lambda","eps")
      
      #Initalization function:
      count.inits <- function (){
        list(beta=rnorm(degree), offset=rnorm(1), sig.beta=runif(1), eps=rnorm(K))
      }    
      
    } else if(interceptQ==FALSE & overdispersionQ==TRUE) { # Model with overdispersion terms
      
      #Get model file path:
      path.to.model.file <- system.file("jags","hierarchical.overdispersed.poission.regression.bug", package="fdrID")
      
      #Parameters for the model:
      count.parameters <- c ("beta","sig.beta","lambda","eps")
      
      #Initalization function:
      count.inits <- function (){
        list(beta=rnorm(degree), sig.beta=runif(1), eps=rnorm(K))
      }    
      
    } else if(interceptQ==TRUE & overdispersionQ==FALSE) { # Model with intercept term
      
      #Get model file path:
      path.to.model.file <- system.file("jags","hierarchical.poission.regression.with.intercept.bug", package="fdrID")
      
      #Parameters for the model:
      count.parameters <- c ("offset","beta","sig.beta","lambda")
      
      #Initalization function:
      count.inits <- function (){ 
        list(beta=rnorm(degree), offset=rnorm(1), sig.beta=runif(1))
      }    
      
    } else if(interceptQ==FALSE & overdispersionQ==FALSE) { # Model with NO intercept OR overdispersion terms
      
      #Get model file path:
      path.to.model.file <- system.file("jags","hierarchical.poission.regression.bug", package="fdrID")
      
      #Parameters for the model:
      count.parameters <- c ("beta","sig.beta","lambda")
      
      #Initalization function:
      count.inits <- function (){
        list(beta=rnorm(degree), sig.beta=runif(1))
      }    
    
    } else {
      
      stop("Model Not Implemented!")
      
    }
    
  }
  
  sampler.info <- list(count.data, count.parameters, count.inits, path.to.model.file, z.values, x)
  names(sampler.info) <- c("Data", "Model.Parameters", "Initialization.Function", "BUG.Model.File.Path", "z.Values", "Bin.Midpoints")
  return(sampler.info)
  
#  if(sampler == "stan") {
#    
#    inits <- function (){
#      list(beta=rnorm(degree), offset=rnorm(1), sigma_sq_beta=runif(1))
#    }
#    
#  } 
  
}