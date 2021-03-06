#--------------------------------------------
#' @title null.sn.fit
#' 
#' @description Skewed normal fit for null or transformed null similarity scores.
#' 
#' @details The skewed normal distribution is a flexible three parameter distribution which often gives 
#' a good fit to the null log similarity scores when these scores are on a scale of 0 to 1. Other than 
#' standardization, this function assumes that if the user wants the scores to be transformed, they've 
#' transformed them. This routine calls the \code{selm} function from the sn package.
#' 
#' @param score.null.vec  Vector of null (non-match) similarity scores or transformed similarity scores
#' @param standardizeQ    Whether or not to standardize the null scores
#' @param plotQ           Diagnostic plots?
#' 
#' @return  list with the fitted parameters, fit info and chi-square goodness of fit test results
#'
#' @references XXXX
#'
#' @examples
#' XXXX
#--------------------------------------------
null.sn.fit <- function(score.null.vec, standardizeQ=FALSE, plotQ=FALSE) {
  
  lgs <- score.null.vec
  
  if(standardizeQ==TRUE){
    lgs <- (lgs - mean(lgs))/sd(lgs)
  }
  
  lgs <- sort(lgs)
  
  sn.fit <- selm(lgs~1)
  dp.vec <- coef(sn.fit, param.type = "DP")
  xi.est <- dp.vec[1]
  om.est <- dp.vec[2]
  alp.est <- dp.vec[3]
  
  dens <- dsn(lgs, xi=xi.est, omega=om.est, alpha=alp.est, tau=0)
  
  #Compute AIC and BIC for comparison to other fits
  llk <- sum(log(dens))
  N<-length(lgs)
  k<-3
  aic <- -2* llk + (2*N*k/(N-k-1))
  bic <- -2* llk + k*log(N)
  
  #This is needed for both plots and fit diagnostics:
  lgs.hist.info <- hist(lgs,plot=F)
  
  if(plotQ==TRUE){
    
    print("Rendering diagnostic plots...")
    
    split.screen( figs = c( 1, 2 ) )
    
    screen(1)
    
    ylim.max <- max(dens,lgs.hist.info$density)
    xlim.max <- max(lgs,lgs.hist.info$breaks)
    xlim.min <- min(lgs,lgs.hist.info$breaks)
    
    if(standardizeQ==TRUE){
      fittitle <- "SN fit to standardized log(null)"
    } else {
      fittitle <- "SN fit to log(null)"
    }
    
    plot(lgs, dens, xlim=c(xlim.min,xlim.max), ylim=c(0,ylim.max), typ="l", col="blue", lwd=3, xlab="", ylab="")
    par(new=T)
    hist(lgs, probability=T, xlim=c(xlim.min,xlim.max), ylim=c(0,ylim.max), xlab="STD log(KNM)", main=fittitle)
    
    #Q-Q Plot:
    #t-axis:
    tmax<-1000
    tax<-seq(1,tmax,1)/(tmax+1)
    cemp<-ecdf(lgs)(lgs)
    #Empirical quantile function (inverse CDF):
    qemp<-splinefun(cemp,lgs)
    #Empirical quantiles:
    Zt<-qemp(tax)
    #Quantiles from fit:
    Zt.hat <- qsn(tax, xi=xi.est, omega=om.est, alpha=alp.est, tau=0)
    
    if(standardizeQ==TRUE){
      QQtitle <- "Q-Q plot for SN fit to standardized log(null) hist"
    } else {
      QQtitle <- "Q-Q plot for SN fit to log(null) hist"
    }
    
    #Q-Q plot:
    screen(2)
    plot(Zt,Zt.hat, xlab="empirical quantiles", ylab="fit quantiles",main=QQtitle)
    abline(0,1)
    
    close.screen( all = TRUE )
    
  }
  
  freq.obs <- lgs.hist.info$counts
  freq.expec <- rep(-1,length(lgs.hist.info$mids))
  fit.probs <- rep(-1,length(lgs.hist.info$mids))
  
  print("Computing fit interquantile probabilities...")
  for(i in 1:(length(lgs.hist.info$breaks)-1)){
    
    upi <- lgs.hist.info$breaks[i+1]
    loi <- lgs.hist.info$breaks[i]
    
    fit.probs[i] <- psn(upi, xi=xi.est, omega=om.est, alpha=alp.est, tau=0) - psn(loi, xi=xi.est, omega=om.est, alpha=alp.est, tau=0) 
    freq.expec[i] <- fit.probs[i] * length(score.null.vec)
  }
  
  plt <- psn(lgs.hist.info$breaks[1], xi=xi.est, omega=om.est, alpha=alp.est, tau=0)
  prt <- 1-sum(c(plt,fit.probs))
  fit.probs <- c(plt,fit.probs,prt)
  
  fit.info <- cbind(fit.probs, c(plt*length(score.null.vec), freq.expec, prt*length(score.null.vec)), c(0,freq.obs,0))
  colnames(fit.info) <- c("interquant.probs", "interquant.exp.cts", "interquant.obs.cts")
  chisq.results <- chisq.test(c(0,freq.obs,0), p = fit.probs)
  
  fit.params <- c(xi.est, om.est, alp.est)
  names(fit.params) <- c("xi.hat", "omega.hat", "alpha.hat")
  
  info.list <- list(fit.params, fit.info, chisq.results, sn.fit, aic, bic)
  names(info.list) <- c("parameters", "fit.info", "chi.square.test", "fit.obj", "AIC", "BIC")
  
  return(info.list)
  
}