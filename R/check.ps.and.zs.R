#--------------------------------------------
#' @title check.ps.and.zs
#' @description Run diagnostic tests on p-values and z-values for both null (non-match) and non-null (match) scores.
#' 
#' @details Make a nice plot of overlapped histograms for each set of scores. NOTE: the y-axis is density not frequency.
#'
#' @param null.p.values    The p-values from the null (non-match) scores.
#' @param nonnull.p.values The p-values from the non-null (match) scores
#' @param plotQ            Diagnostic plots?
#' 
#' @examples
#' XXXX
#--------------------------------------------
check.ps.and.zs <- function(null.p.values, nonnull.p.values, plotQ=FALSE) {
  
  #Check p-values for major problems:
  too.big.pos.null.zs <- which(qnorm(null.p.values)==Inf)
  too.big.neg.null.zs <- which(qnorm(null.p.values)==-Inf)
  too.big.pos.nonnull.zs <- which(qnorm(nonnull.p.values)==Inf)
  too.big.neg.nonnull.zs <- which(qnorm(nonnull.p.values)==-Inf)
  
  if(length(too.big.pos.null.zs)> 0){
    print(paste("These", length(too.big.pos.null.zs),"null p-values lead to +INFINITE z-values:"))
    print(too.big.pos.null.zs)
  }
  if(length(too.big.neg.null.zs)> 0){
    print(paste("These", length(too.neg.pos.null.zs),"null p-values lead to -INFINITE z-values:"))
    print(too.big.neg.null.zs)
  }
  if(length(too.big.pos.nonnull.zs)> 0){
    print(paste("These", length(too.big.pos.nonnull.zs),"non-null p-values lead to +INFINITE z-values:"))
    print(too.big.pos.nonnull.zs)
  }
  if(length(too.big.neg.nonnull.zs)> 0){
    print(paste("These", length(too.neg.pos.nonnull.zs),"null p-values lead to -INFINITE z-values:"))
    print(too.big.neg.nonnull.zs)
  }
  
  #Drop any problem null p-values before moving on with testing:
  print("If there are any problem null p-values, they are being dropped before further diagnostic tests...")
  if(length(c(too.big.pos.null.zs, too.big.neg.null.zs))> 0) {
    cleaned.null.p.values <- null.p.values[-c(too.big.pos.null.zs, too.big.neg.null.zs)]
  } else {
    cleaned.null.p.values <- null.p.values
  }
    
  
  #Check to see if there is evidence that the null p-vals are not 0,1 uniform
  ks.results <-  ks.test(cleaned.null.p.values,"punif",0,1)
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
  
  chisq.results <- pearson.test(cleaned.null.z.values)
  chisq.pval <- chisq.results$p.value
  
  sf.results <- sf.test(sample(cleaned.null.z.values,5000))
  sf.pval <- sf.results$p.value
  
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
  
}