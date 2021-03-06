#--------------------------------------------
#' @title process.f
#' @description Computes estimates of the posterior null probabilities using f(x), the denominator of Bayes'
#' Theorem in Efron's two-groups methodology.
#' 
#' @details #Core function of the package. This algorithm take points making up the denominator
#' of Bayes theorem in in Efron's methodology f(x), computes
#' the likelihood function f0(x) = p(x|null) and (empirical)
#' estimate of the prior p0. One or many estimates of f(x) can be contained (as columns) in lam.mat.
#'
#' The algorithim is largely taken from locfdr code. Follows basically
#' Efron 2006 Size, Power and False Discovery Rates
#'
#' NOTE: Rows of lam.mat are f's generated by a Bayesian Poisson
#' regression onto the hisogram of counts generated by z-scores.
#'
#' @param zscores Efron's z-scores derived from null and non-null score distributions
#' @param xmids Midpoints of z-score histogram cells
#' @param pct0 Taken from locfdr package. Proportion of the z-score distribution used in fitting the null density f0(z) by central matching. If a 2-vector, e.g. pct0=c(0.25,0.60), the range [pct0[1], pct0[2]] is used. If a scalar, [pct0, 1-pct0] is used.
#' @param lam.mat Point estimate(s) of f(x) (lambda(x) in the Poission regression). Each column of lam.mat is an extimate of f(x). In a Bayesian regression context, each column is a draw of f(x) from the joint posterior distribution.
#' 
#' @return A list consisting of:
#' @return delta0, null mean estimates
#' @return sig0, null standard dev estimates
#' @return p0, null prior estimates
#' @return f(x) estimates.
#'
#' 
#' @references locfdr package by Efron, Turnbull and Narasimhan. 
#' Archived at http://cran.r-project.org/src/contrib/Archive/locfdr/
#'
#' @examples
#' XXXX
#--------------------------------------------
process.f<-function(zscores,xmids,pct0,lam.mat){
  
  num.sims<-nrow(lam.mat)
  KK<-ncol(lam.mat)
  
  #From locfdr
  pctup <- 1 - pct0
  pctlo <- pct0
  lo0 <- quantile(zscores, pctlo)
  hi0 <- quantile(zscores, pctup)
  nx <- length(xmids)
  i0 <- (1.:nx)[xmids > lo0 & xmids < hi0]
  x0 <- xmids[i0]
  
  delta0.vec<-rep(NA,num.sims)
  sig0.vec<-rep(NA,num.sims)
  p0.vec<-rep(NA,num.sims)
  #p00.vec<-rep(NA,num.sims)
  fdr.mat<-array(NA,dim(lam.mat))
  
  #*********Blast this over all requested processes????
  
  for(i in 1:num.sims) {
    #From locfdr:
    f<-lam.mat[i,]
    l <- log(f)  
    y0 <- l[i0]
    #Shift axis so central peak is at z=0. Cf. Efron 2006
    imax <- seq(l)[l == max(l)][1]
    xmax <- xmids[imax]
    X00 <- cbind(x0 - xmax, (x0 - xmax)^2)
    
    #Parabolic fit to central peak of f (i.e. use central matching alg described in Efron 2006 and locfdr):
    lr <- lm(y0 ~ X00)
    co <- lr$coef
    X0 <- cbind(1, x - xmax, (x - xmax)^2)
    delta0 <-  - co[2.]/(2. * co[3.]) + xmax #Same as Efron 2006
    sig0 <- 1./sqrt(-2. * co[3.])            #Same as Efron 2006
    delta0.vec[i] <- delta0
    sig0.vec[i] <- sig0
        
    l0 <- as.vector(X0 %*% co)
    f0p <- exp(l0)
    p0 <- sum(f0p)/sum(f)
    p0.vec[i] <- p0 
    
    fdr.sim <- pmin((f0p)/f, 1)
    fdr.mat[i,] <- fdr.sim
    
  }
  
  sim.proc.info<-list(delta0.vec, sig0.vec, p0.vec,fdr.mat)
  names(sim.proc.info) <- c("delta0","sigma0","p0","fdr")
  
  return(sim.proc.info)
  
}