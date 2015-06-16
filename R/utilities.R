#--------------------------------------------
#' @title split.scores
#' @description Split a set of scores into null and non-null sets
#' 
#' @details Split a set of scores into null (non-match) and non-null (match) sets. NOTE: One non-null score is assumed per row in \code{score.mat}
#'
#' @param score.mat A matrix of null (non-match) and non-null (match) scores. 
#' @param lbls      Correct ID labels.
#' 
#' @return a list of null and non-null scores
#' 
#' @examples
#' XXXX
#--------------------------------------------
split.scores<-function(score.mat,lbls)
{
  #extract vector of known match (non-null) scores:
  score.km<-score.mat[cbind(1:nrow(score.mat),as.numeric(lbls))]
  nonnull.vec<-as.numeric(score.km)
  print(paste("Number of non-null scores: ", length(nonnull.vec)))
  
  #build up vector of known non-match (null) scores:
  score.knm<-NULL
  for(j in 1:nrow(score.mat))
  {
    tmp.prob.vec<-score.mat[j,]
    tmp.prob.vec<-tmp.prob.vec[-as.numeric(lbls[j])]  #here, there are k-1 knm scores per observation. Take mean/median/max, something else  here??
    score.knm<-c(score.knm,tmp.prob.vec)
  }
  null.vec<-as.numeric(score.knm)
  print(paste("Number of null scores: ", length(null.vec)))
  
  return(list(null.vec,nonnull.vec))
  
}


#--------------------------------------------
#' @title scores.histograms
#' @description Make a nice plot of overlapped histograms for each set of scores.
#' 
#' @details Make a nice plot of overlapped histograms for each set of scores. NOTE: the y-axis is density not frequency.
#'
#' @param null.vec Non-matching scores
#' @param nonnull.vec Matching scores
#' @param ylim.max Max on the y-axis. Helps beautify the overlayed histograms of null/non-null scores
#' 
#' @examples
#' XXXX
#--------------------------------------------
scores.histograms<-function(null.vec, nonnull.vec, xlim.min=min(c(null.vec, nonnull.vec)), xlim.max=max(c(null.vec, nonnull.vec)), ylim.max=20) {
  
  hist(null.vec, probability=TRUE, xlim=c(xlim.min,xlim.max), ylim=c(0,ylim.max), col=rgb(1,0,0,1/4), main="",xlab="")
  par(new=TRUE)
  hist(nonnull.vec, probability=TRUE, xlim=c(xlim.min,xlim.max), ylim=c(0,ylim.max), col=rgb(0,1,0,1/4), main="Null(KNM, red) and Non-null(KM, green)", xlab="Score")
  
}

#--------------------------------------------
#' @title hist.mode
#' @description Handy mode function.
#' 
#' @details A handy mode function for a vector of data. Useful for poking around a troublesome non-null distribution 
#' of scires
#'
#' @param x data 
#' 
#' @examples
#' XXXX
#--------------------------------------------
hist.mode <- function(scores, bre, plotQ=FALSE) {
  
  #Taken from locfdr:
  lo <- min(scores)
  up <- max(scores)
  
  zzz <- pmax(pmin(scores, up), lo)
  
  breaks <- seq(lo, up, length = bre)
  hist.info <- hist(zzz, breaks = breaks, plot = F)
  
  bin.cts <- hist.info$counts
  bin.mids <- hist.info$mids
  
  #Mode of the histogram
  mde <- bin.mids[which(bin.cts==max(bin.cts))]
    
  if(plotQ==T) {
    hist(zzz, breaks=breaks)
  }
  
  freqs <- cbind(hist.info$mids, hist.info$counts)
  colnames(freqs) <- c("bin.mids", "counts")
  
  return(list(freqs, mde)) 
}