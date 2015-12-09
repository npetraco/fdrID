#--------------------------------------------
#' @title null.gevd.fit
#' 
#' @description Generalized extreme value fit for null similarity scores.
#' 
#' @details The generalized extreme value distribution is a flexible three parameter distribution which often gives 
#' a good fit to the null log similarity scores when these scores are on a scale of 0 to 1. Other than 
#' standardization, this function assumes that if the user wants the scores to be transformed, they've 
#' transformed them. This routine calls the \code{fgev} function from the evd package.
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
null.gev.fit <- function(score.null.vec, standardizeQ=FALSE, plotQ=FALSE) {
  
  #Take the log:
  lgs <- score.null.vec
  
  if(standardizeQ==TRUE){
    lgs <- (lgs - mean(lgs))/sd(lgs)
  }
  
  lgs <- sort(lgs)
  
  evFit.lgs <- fgev(lgs)
  
  loc.est <- evFit.lgs$param[1]
  scale.est <- evFit.lgs$param[2]
  shape.est <- evFit.lgs$param[3]
  
  dens <- dgev(lgs, loc=loc.est, scale=scale.est, shape=shape.est)
  
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
      fittitle <- "GEV fit to standardized log(null)"
    } else {
      fittitle <- "GEV fit to log(null)"
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
    Zt.hat<-qgev(tax, loc=loc.est, scale=scale.est, shape=shape.est)
    
    if(standardizeQ==TRUE){
      QQtitle <- "Q-Q plot for GEV fit to standardized log(null) hist"
    } else {
      QQtitle <- "Q-Q plot for GEV fit to log(null) hist"
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
    
    fit.probs[i] <- pgev(upi, loc=loc.est, scale=scale.est, shape=shape.est) - pgev(loi, loc=loc.est, scale=scale.est, shape=shape.est)
    freq.expec[i] <- fit.probs[i] * length(score.null.vec)
    
    #print(paste("From:", loi, "to", upi, "prob:", fit.probs[i], "expec cnt:", freq.expec[i]))
    #print(paste("Mid:", upi - loi, "prob:", fit.probs[i], "expec cnt:", freq.expec[i]))
  }
  
  plt <- pgev(lgs.hist.info$breaks[1], loc=loc.est, scale=scale.est, shape=shape.est)
  
  prt <- 1-sum(c(plt,fit.probs))
  fit.probs <- c(plt,fit.probs,prt)
  
  #print(lgs.hist.info$mids)
  
  fit.info <- cbind(c(-Inf, lgs.hist.info$mids, Inf), fit.probs, c(plt*length(score.null.vec), freq.expec, prt*length(score.null.vec)), c(0,freq.obs,0))
  colnames(fit.info) <- c("interquant.mid","interquant.probs", "interquant.exp.cts", "interquant.obs.cts")
  rownames(fit.info) <- NULL
  chisq.results <- chisq.test(c(0,freq.obs,0), p = fit.probs)
  
  fit.params <- c(loc.est, scale.est, shape.est)
  names(fit.params) <- c("loc.hat", "scale.hat", "shape.hat")
  
  info.list <- list(fit.params, fit.info, chisq.results, evFit.lgs, aic, bic)
  names(info.list) <- c("parameters", "fit.info", "chi.square.test", "fit.obj", "AIC", "BIC")
  
  return(info.list)
  
}