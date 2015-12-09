#--------------------------------------------
#' @title null.nig.fit
#' 
#' @description Normal inverse gaussian fit for null similarity scores.
#' 
#' @details The normal inverse gaussian distribution is a flexible four parameter distribution which often gives 
#' a good fit to the null log similarity scores when these scores are on a scale of 0 to 1. Other than 
#' standardization, this function assumes that if the user wants the scores to be transformed, they've 
#' transformed them. This routine calls the \code{fit.NIGuv} function from the ghyp package.
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
null.nig.fit <- function(score.null.vec, max.iter = 1000, standardizeQ=FALSE, plotQ=FALSE) {
  
  #Take the log:
  lgs <- score.null.vec
  
  if(standardizeQ==TRUE){
    lgs <- (lgs - mean(lgs))/sd(lgs)
  }
  
  lgs <- sort(lgs)
  
  #nigFit.lgs <- fBasics::nigFit(lgs, alpha = alpha.init, beta = beta.init, delta = delta.init, mu = mu.init, doplot = FALSE, method="mle")
  nigFit.lgs <- fit.NIGuv(lgs, control=list(maxit=max.iter))
  if(nigFit.lgs@converged == TRUE) {
    print("Converged.")
  } else {
    stop("NIG fit did not converge!!!! Turn up the number of iterations.")
  }
  
  
  alpb.est <- nigFit.lgs@alpha.bar
  mu.est <- nigFit.lgs@mu
  sig.est <- nigFit.lgs@sigma
  gam.est <- nigFit.lgs@gamma
  lam.est <- nigFit.lgs@lambda  
  #alp.est <- nigFit.lgs@fit$estimate[1]
  #bet.est <- nigFit.lgs@fit$estimate[2]
  #del.est <- nigFit.lgs@fit$estimate[3]
  #mu.est <- nigFit.lgs@fit$estimate[4]
  
  #This is needed for both plots and fit diagnostics:
  lgs.hist.info <- hist(lgs,plot=F)
  
  #dens <- fBasics::dnig(lgs, alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est)
  dens <- dghyp(lgs, object = nigFit.lgs)
  
  #Compute AIC and BIC for comparison to other fits
  #llk <- sum(log(dens))
  llk <- nigFit.lgs@llh
  N<-length(lgs)
  k<-4
  #aic <- -2* llk + (2*N*k/(N-k-1))
  aic <- AIC(nigFit.lgs)
  bic <- -2* llk + k*log(N)
  
  
  if(plotQ==TRUE){
    
    print("Rendering diagnostic plots...")
    
    split.screen( figs = c( 1, 2 ) )
    
    screen(1)
    
    ylim.max <- max(dens,lgs.hist.info$density)
    xlim.max <- max(lgs,lgs.hist.info$breaks)
    xlim.min <- min(lgs,lgs.hist.info$breaks)
    
    if(standardizeQ==TRUE){
      fittitle <- "NIG fit to standardized log(null)"
    } else {
      fittitle <- "NIG fit to log(null)"
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
    #Zt.hat <- fBasics::qnig(tax,alpha=alp.est,beta=bet.est,delta=del.est,mu=mu.est)
    Zt.hat <- qghyp(tax, object = nigFit.lgs)
    
    if(standardizeQ==TRUE){
      QQtitle <- "Q-Q plot for NIG fit to standardized log(null) hist"
    } else {
      QQtitle <- "Q-Q plot for NIG fit to log(null) hist"
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
    
    #fit.probs[i] <- fBasics::pnig(upi, alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est) - fBasics::pnig(loi, alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est)
    fit.probs[i] <- pghyp(upi, object = nigFit.lgs) - pghyp(loi, object = nigFit.lgs)
    freq.expec[i] <- fit.probs[i] * length(score.null.vec)
  }
  
  #plt <- fBasics::pnig(lgs.hist.info$breaks[1], alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est)
  plt <- pghyp(lgs.hist.info$breaks[1], object = nigFit.lgs)
  prt <- 1-sum(c(plt,fit.probs))
  fit.probs <- c(plt,fit.probs,prt)
  
  
  fit.info <- cbind(fit.probs, c(plt*length(score.null.vec), freq.expec, prt*length(score.null.vec)), c(0,freq.obs,0))
  colnames(fit.info) <- c("interquant.probs", "interquant.exp.cts", "interquant.obs.cts")
  chisq.results <- chisq.test(c(0,freq.obs,0), p = fit.probs)
  
  fit.params <- c(alpb.est, mu.est, sig.est, gam.est, lam.est)
  names(fit.params) <- c("alpha.bar.hat", "mu.hat", "sigma.hat", "gamma.hat", "lambda.hat")
  
  info.list <- list(fit.params, fit.info, chisq.results, nigFit.lgs, aic, bic)
  names(info.list) <- c("parameters", "fit.info", "chi.square.test", "fit.obj", "AIC", "BIC")
  
  return(info.list)
  
}