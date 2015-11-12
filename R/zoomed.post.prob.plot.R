#--------------------------------------------
#' @title zoomed.post.prob.plot
#' @description XXXX
#' 
#' @details XXXX
#'
#' @param XXXX
#' 
#' @return XXXX
#'
#' @examples
#' XXXX
#--------------------------------------------
zoomed.post.prob.plot <- function(zscores, zbounds=c(NA,NA), point.est.func, upper.est.func, lower.est.func, prob.scale="percent", xlab=NULL, ylab=NULL, main=NULL) {
  
  #put in a check for NAs in zbounds
  
  if(is.na(zbounds[1]) & !is.na(zbounds[2])) {
    
    zidxs <- which(zscores < zbounds[2])
    
  } else if(!is.na(zbounds[1]) & is.na(zbounds[2])) {
    
    zidxs <- which(zscores > zbounds[1])
    
  } else if(is.na(zbounds[1]) & is.na(zbounds[2])) {
    
    zidxs <- 1:length(zscores)
    
  }  else if(!is.na(zbounds[1]) & !is.na(zbounds[2])) {
    
    zidxs <- which( (zscores >= zbounds[1]) & (zscores <= zbounds[2]) )
    
  }
  
  zs <- sort(zscores[zidxs])
  pe <- point.est.func(zs)
  ue <- upper.est.func(zs)
  le <- lower.est.func(zs)
  
  if(prob.scale=="percent") {
    pe <- 100 * pe
    ue <- 100 * ue
    le <- 100 * le
  } else {
    
  }
  
  ymax <- max(ue)
  ymin <- min(le)
  
  if(is.null(xlab)) {
    xlab.txt <- "z-score"
  } else {
    xlab.txt <- xlab
  }
  
  if(is.null(ylab)) {
    ylab.txt <- "posterior prob."
  } else {
    ylab.txt <- ylab
  }
  
  if(is.null(main)) {
    main.txt <- ""
  } else {
    main.txt <- main
  }
  
  plot(zs, pe, typ="l", ylim=c(ymin,ymax), xlab = xlab.txt, ylab=ylab.txt, main=main.txt)
  lines(zs, ue, col="red")
  lines(zs, le, col="blue")

}