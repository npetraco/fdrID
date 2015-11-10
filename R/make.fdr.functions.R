#--------------------------------------------
#' @title make.fdr.functions
#' @description XXXX
#' 
#' @details XXXX
#'
#' @param XXXX
#' 
#' @return XXXX
#'
#' @examples
#' XXXX
#--------------------------------------------
make.fdr.functions <- function(z.values, bin.midpoints, pct0=0.25, posterior.f.values, credibility.level=0.95, interval.type="equal.tail", plotQ=FALSE) {
  
  processed.fdr.info <- process.f(z.values, bin.midpoints, pct0=pct0, posterior.f.values)
  
  #Get the posterior point fdr(z) estimates 
  posterior.fdr.draws <- processed.fdr.info[["fdr"]]
  
  K <- ncol(posterior.fdr.draws)
  print(K)
  p0 <- processed.fdr.info[["p0"]]
  delta0 <- processed.fdr.info[["delta0"]]
  sigma0 <- processed.fdr.info[["sigma0"]]
  
  posterior.fdr.means <- colMeans(posterior.fdr.draws)
  
  posterior.fdr.medians<-apply(posterior.fdr.draws,2,median)
  
  if(interval.type=="equal.tail") {
    
    signif <- (1-(credibility.level))/2 #Two sided interval equal tail interval
    cred.ints.fdrs <- t(sapply(1:K,function(xx){quantile(posterior.fdr.draws[,xx],c(signif,(1-signif)))}))
    print(cred.ints.fdrs)    
    
  } else if (interval.type=="hpd") {
    
    cred.ints.fdrs <- HPDinterval(as.mcmc(posterior.fdr.draws, prob=credibility.level))    
    print(cred.ints.fdrs)
    
  } else {
    stop("Interval Type Not Implemented!")
  }
  
  #Local false discovery rate interpolation functions:
  fdr.means.func <- splinefun(bin.midpoints, posterior.fdr.means)
  fdr.medians.func <- splinefun(bin.midpoints, posterior.fdr.medians)
  
  upper.fdr.ci.func <- splinefun(bin.midpoints, cred.ints.fdrs[,2])
  lower.fdr.ci.func <- splinefun(bin.midpoints, cred.ints.fdrs[,1])
  
  #Local true discovery rate interpolation functions:
  tdr.means.func <- splinefun(bin.midpoints, 1 - posterior.fdr.means)
  tdr.medians.func <- splinefun(bin.midpoints, 1 - posterior.fdr.medians)
  
  upper.tdr.ci.func <- splinefun(bin.midpoints, 1 - cred.ints.fdrs[,2])
  lower.tdr.ci.func <- splinefun(bin.midpoints, 1 - cred.ints.fdrs[,1])
  
  if(plotQ==TRUE) {
    plot(bin.midpoints, cred.ints.fdrs[,2]*100,col="red",typ="l",lwd=2,ylim=c(0,100.1), xlab="z",ylab="fdr(%)")
    lines(bin.midpoints, posterior.fdr.means*100,ylim=c(0,100.1))
    lines(bin.midpoints, cred.ints.fdrs[,1]*100,col="blue",typ="l",lwd=2,ylim=c(0,100.1))
  }

  fdr.info <- list(
    fdr.means.func,
    fdr.medians.func,
    upper.fdr.ci.func,
    lower.fdr.ci.func,
    tdr.means.func,
    tdr.medians.func,
    upper.tdr.ci.func,
    lower.tdr.ci.func,
    p0,
    delta0,
    sigma0)
  
  names(fdr.info) <- c("fdr.means.func","fdr.medians.func",
                       "upper.fdr.ci.func","lower.fdr.ci.func",
                       "tdr.means.func","tdr.medians.func",
                       "upper.tdr.ci.func","lower.tdr.ci.func",
                       "p0","delta0","sigma0")
  
  return(fdr.info)
  
  
}