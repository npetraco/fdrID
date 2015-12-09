#--------------------------------------------
#' @title check.ps.and.zs
#' @description Run diagnostic tests on p-values and z-values for both null (non-match) and non-null (match) scores.
#' 
#' @details Check to see if: a. Null z-values are = +/- inf, Non-Null z-values are = -/+. If things are operating 
#' according to plan (cf. Efron 2006) the Null z-values should be between +/- 3 with ~ 99% probability 
#' (i.e. z_Null ~ N(0,1)  approximately). Extreme Null z-values are pathological. Check the KNM scores and/or
#' density fit to the KNM scores for weird-ness. It may be due to overly thin tails.
#'
#' @param null.p.values    Optional. The p-values from the null (non-match, KM) scores.
#' @param nonnull.p.values Optional. The p-values from the non-null (match, KNM) scores.
#' @param plotQ            Diagnostic plots?
#' 
#' @return A list with indices of problem p/z-values and a list containing statistical test results for 
#' uniformity(p-values)/normality(z-values). 
#'
#' @examples
#' XXXX
#--------------------------------------------
check.ps.and.zs <- function(null.p.values = NULL, nonnull.p.values = NULL, printQ=FALSE, plotQ=FALSE) {
  
  #Check p-values for major problems:
  if(!is.null(null.p.values)) {
    too.big.pos.null.zs <- which(qnorm(null.p.values)==Inf)
    too.big.neg.null.zs <- which(qnorm(null.p.values)==-Inf)
  } else {
    too.big.pos.null.zs <- NULL
    too.big.neg.null.zs <- NULL
  }
  
  if(!is.null(nonnull.p.values)) {
    too.big.pos.nonnull.zs <- which(qnorm(nonnull.p.values)==Inf)
    too.big.neg.nonnull.zs <- which(qnorm(nonnull.p.values)==-Inf)
  } else {
    too.big.pos.nonnull.zs <- NULL
    too.big.neg.nonnull.zs <- NULL
  }
    
  if(printQ == TRUE) {
    
    print("NOTE: If there are any problem null p-values, they are being dropped before further diagnostic tests...")
    
    if(length(too.big.pos.null.zs)> 0){
      print(paste("These", length(too.big.pos.null.zs),"null p-values lead to +INFINITE z-values:"))
      print(too.big.pos.null.zs)
    }
    if(length(too.big.neg.null.zs)> 0){
      print(paste("These", length(too.big.neg.null.zs),"null p-values lead to -INFINITE z-values:"))
      print(too.big.neg.null.zs)
    }
    if(length(too.big.pos.nonnull.zs)> 0){
      print(paste("These", length(too.big.pos.nonnull.zs),"non-null p-values lead to +INFINITE z-values:"))
      print(too.big.pos.nonnull.zs)
    }
    if(length(too.big.neg.nonnull.zs)> 0){
      print(paste("These", length(too.big.neg.nonnull.zs),"null p-values lead to -INFINITE z-values:"))
      print(too.big.neg.nonnull.zs)
    }
    
  }
    
  #Drop any problem null p-values before moving on with testing:
  if(length(c(too.big.pos.null.zs, too.big.neg.null.zs))> 0) {
    cleaned.null.p.values <- null.p.values[-c(too.big.pos.null.zs, too.big.neg.null.zs)]
  } else {
    cleaned.null.p.values <- null.p.values
  }
    
  
  #Check to see if there is evidence that the null p-vals are not 0,1 uniform
  options(warn=-1) #Shut off annoying warning about ties in ks test
  ks.results <-  ks.test(cleaned.null.p.values,"punif",0,1)
  options(warn=0)
  ks.pval <- ks.results$p.value
  #print(ks.results)
  
  #Check to see if there is evidence that the null z-vals are not Gaussian
  cleaned.null.z.values <- qnorm(cleaned.null.p.values)
  
  shapiro.results <- shapiro.test(sample(cleaned.null.z.values,5000))
  shapiro.pval <- shapiro.results$p.value
  #print(shapiro.results) 
  
  #From fBasics:
  jb.results <- jbTest(cleaned.null.z.values)
  #jb.pval <- p.value(jb.results)
  #print(jb.results)
  
  dago.results <- dagoTest(cleaned.null.z.values)
  #print(dago.results)
  #dago.pval <- dago.results$p.value
  
  #From nortest:
  ad.results <- ad.test(cleaned.null.z.values)
  ad.pval <- ad.results$p.value
  
  cvm.results <- cvm.test(cleaned.null.z.values)
  cvm.pval <- cvm.results$p.value
  
  lillie.results <- lillie.test(cleaned.null.z.values)
  lillie.pval <- lillie.results$p.value
  
  options(warn=-1) #Turn off chisq warning that it may not work. 
  chisq.results <- pearson.test(cleaned.null.z.values)
  options(warn=0)
  chisq.pval <- chisq.results$p.value
  
  
  sf.results <- sf.test(sample(cleaned.null.z.values,5000))
  sf.pval <- sf.results$p.value
  
  if(printQ == TRUE) {
    
    print("Null test p-values:")
    print("-------------------")
    print(paste("K-S (uniform):       ",ks.pval))
    print(paste("Shapiro (normality): ",shapiro.pval))
    #print(paste("Adjusted J-B:",jb.pval))
    #print(paste("D'Agostino:  ",dago.pval))
    print(paste("A-D (normality):     ",ad.pval))
    print(paste("C-vM (normality):    ",cvm.pval))
    print(paste("Lillie (normality):  ",lillie.pval))
    print(paste("Chi-Sq (normality):  ",chisq.pval))
    print(paste("S-F (normality):     ",sf.pval))
    
  }
  
  if(plotQ==TRUE) {
    split.screen( figs = c( 2, 2 ) )
    
    screen(1)
    hist(cleaned.null.p.values, xlab="p", main="Null p-values",col=rgb(1,0,0,1/2))
    
    screen(2)
    hist(c(cleaned.null.p.values, nonnull.p.values), probability=TRUE, xlim=c(0,1), main="All p-values",xlab="p")
    
    screen(3)
    hist(cleaned.null.z.values, xlab="z", main="Null z-values",col=rgb(1,0,0,1/2))
    
    screen(4)
    qqnorm(cleaned.null.z.values, main="Q-Q Null z-values")
    qqline(cleaned.null.z.values)
    
    close.screen( all = TRUE )
    
  }
  
  #Later update the is an S4 class:
  
  diagnostic.info <- list(
    too.big.pos.null.zs, 
    too.big.neg.null.zs, 
    too.big.pos.nonnull.zs, 
    too.big.neg.nonnull.zs,
    list(ks.results,
        shapiro.results,
        jb.results,
        dago.results,
        ad.results,
        cvm.results,
        lillie.results,
        chisq.results,
        sf.results)
       )
  
  names(diagnostic.info) <- c("null.posinf.idxs",
                             "null.neginf.idxs",
                             "nonnull.posinf.idxs",
                             "nonnull.neginf.idxs",
                             "test.results"
                             )
  
  return(diagnostic.info)
}