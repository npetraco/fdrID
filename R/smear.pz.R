#--------------------------------------------
#' @title smear.extreme.nonnull.pvalues
#' @description XXXX
#' 
#' infinite z-values
#'
#' @param XXXX
#' 
#' @return XXXX
#'
#' @examples
#' XXXX
#--------------------------------------------
smear.extreme.nonnull.pvalues <- function(non.null.pvalues, upper.set.pvalue=1e-8, lower.set.pvalue=3.5e-16, plotQ=FALSE) {
  
  #z-value corresponding to upper part of interval to be smeared out
  ztol<-qnorm(upper.set.pvalue)
  
  non.null.zvalues <- qnorm(non.null.pvalues)
  
  #indices of non-null p-values less than upper.set.pvalue. The p-values at these indices should be smeared out over a narrow negative z-interval:
  pv.not.acceptable.idxs <- which(non.null.pvalues < upper.set.pvalue)
  
  smear.info.mat<-cbind(pv.not.acceptable.idxs, non.null.pvalues[pv.not.acceptable.idxs], non.null.zvalues[pv.not.acceptable.idxs])
  colnames(smear.info.mat)<-c("pval idx", "pval to be replaced", "corresp zval")
  #print(smear.info.mat)
  
  #Look at the p and z values that are "ACCEPTABLE"
  pv.accepable.idxs <- 1:length(non.null.pvalues)
  pv.accepable.idxs <- pv.accepable.idxs[-pv.not.acceptable.idxs]
  no.smear.info.mat <- cbind(pv.accepable.idxs, non.null.pvalues[-pv.not.acceptable.idxs], non.null.zvalues[-pv.not.acceptable.idxs])
  colnames(no.smear.info.mat)<-c("pval idx", "pval to be kept", "corresp zval")
  #print(no.smear.info.mat)
  
  
  #This is the upper bound on the smeared distribution:
  minimal.accepable.actual.pval<-min(non.null.pvalues[-pv.not.acceptable.idxs])
  print(paste("Smallest acceptable actual p-value:", minimal.accepable.actual.pval))
  
  #This better be TRUE:
  if(!(minimal.accepable.actual.pval>upper.set.pvalue)){
    stop("Whoa there! The minimal accetable actual p-value is LARGER than the upper set p-value for the interval to be smeared.")
  }
  
  #uniform spread (smear) the very small pvalues over this interval:
  smear.interval.p <- c(minimal.accepable.actual.pval,lower.set.pvalue)
  smear.interval.z <- qnorm(smear.interval.p)
  print(paste("The unacceptable p-values will be uniformly spread across this interval: [", smear.interval.p[2], ",", smear.interval.p[1], "]", sep=""))
  print(paste("This corresponds to the z interval:                                      [", smear.interval.z[2], ",", smear.interval.z[1], "]", sep=""))
  
  
  #Candidate replacement "small" pvalues:
  replacement.small.non.null.pvals<-runif(length(pv.not.acceptable.idxs),max=minimal.accepable.actual.pval, min=lower.set.pvalue)
  
  #So this is what the z-values for the km would look like:
  candidate.replacement.non.null.zvals<-qnorm(replacement.small.non.null.pvals)
  candidate.non.null.zvals<-c(candidate.replacement.non.null.zvals, non.null.zvalues[-pv.not.acceptable.idxs])
  
  if(plotQ==TRUE) {
    split.screen( figs = c( 1, 3 ) )
    
    screen(1)
    hist(replacement.small.non.null.pvals, xlab="p", main="Replacement non-null p-values")
    
    screen(2)
    hist(candidate.replacement.non.null.zvals, xlab="z", main="Replacement non-null z-values")
    
    screen(3)
    hist(candidate.non.null.zvals, xlab="z", main="All non-null z-values")
    
    close.screen( all = TRUE )    
  }
  
  #Do the replacement:
  non.null.pvalues.mod <- non.null.pvalues
  non.null.pvalues.mod[pv.not.acceptable.idxs] <- replacement.small.non.null.pvals
  
  non.null.pvalue.replacement.info <- list(smear.info.mat, no.smear.info.mat, non.null.pvalues.mod)
  
  names(non.null.pvalue.replacement.info) <- c("Replacement p-value info", "Non-replacement p-value info", "Modified p-values")
  
  return(non.null.pvalue.replacement.info)
  
}