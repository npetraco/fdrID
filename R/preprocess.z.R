#-------------------------------------------------------
#Internal function.
#Count useful things from the histogram of joint z-scores
#Largely taken from locfdr code by Efron, Turnbull and Narasimhan
#-------------------------------------------------------
preprocess.z<-function(zscores,bre) {
  
  #Taken from locfdr:
  lo <- min(zscores)
  up <- max(zscores)
  
  zzz <- pmax(pmin(zscores, up), lo)
  breaks <- seq(lo, up, length = bre)
  zh <- hist(zzz, breaks = breaks, plot = F)
  x <- (breaks[-1] + breaks[ - length(breaks)])/2
  y <- zh$counts
  K <- length(y)
  N <- length(zscores)
  
  zinfo<-list(x,y,K,N)
  names(zinfo)<-c("bin midpoints","counts","num bins", "num scores")
  
  return(zinfo)
  
}