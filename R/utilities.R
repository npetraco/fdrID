#--------------------------------------------
#' @title split,scores
#' @description Split a set of scores into null and non-null sets
#' 
#' @details Split a set of scores into null (non-match) and non-null (match) sets.
#'
#' @param score.mat 
#' @param lbls
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
#' @details Make a nice plot of overlapped histograms for each set of scores.
#'
#' @param null.vec Non-matching scores
#' @param nonnull.vec Matching scores
#' @param ylim.max Max on the y-axis. Helps beautify the overlayed histograms of null/non-null scores
#' 
#' @examples
#' XXXX
#--------------------------------------------
scores.histograms<-function(null.vec, nonnull.vec, ylim.max) {
  
  hist(null.vec,probability=TRUE,xlim=c(0,1),ylim=c(0,ylim.max),col=rgb(1,0,0,1/4),main="",xlab="")
  par(new=TRUE)
  hist(nonnull.vec,probability=TRUE,xlim=c(0,1),ylim=c(0,ylim.max),col=rgb(0,1,0,1/4),main="Null(KNM, red) and Non-null(KM, green)",xlab="Platt-score")
  
}