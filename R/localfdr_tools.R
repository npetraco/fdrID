#--------------------------------------------
#' @title bs.null
#' @description Estimate svm Platt-score null and non-null distributions.
#' 
#' @details Estimate svm Platt-score null and non-null distributions via a group-wise bootstrap.
#' Adapted from Storey and Tibshirani permutation method in PNAS.
#'
#' @param dat         A n by p data frame or matrix of variables. One row per pattern to classify
#' @param lbls        Data ID labels
#' @param nbs         Number of bootstrap iterations.
#' @param svmtyp      Support vector machine type. See e1071 package documentation.
#' @param kern        Kernel type. See e1071 package documentation.
#' @param pparams     Kernel parameters. See e1071 package documentation.
#' @param printQ      Plot boolean
#' @param ylim.max    Plot parameter (max of the y-axis).
#' 
#' @return list of null and non-null scores
#' 
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
bs.null<-function(dat,lbls,nbs,svmtyp="C-classification",kern="linear",pparams=0.1,printQ=FALSE,ylim.max=NULL)
{
#parse lbls into groups:
lbls.idxs<-lapply(1:nlevels(lbls),function(x){which(lbls==levels(lbls)[x])})

null.vec<-NULL
nonnull.vec<-NULL
for(i in 1:nbs)
  {
   #Grab a bootstrap sample, bootstrapping out of each GROUP invididually:
   bsidx<-unlist(as.vector(sapply(1:length(lbls.idxs),function(x){sample(lbls.idxs[[x]],length(lbls.idxs[[x]]), replace=TRUE )})))
#   print(bsidx)
   bsdat<-dat[bsidx,]

   rownames(bsdat)<-NULL #Shuts off an annoying warning
   bslbl<-lbls[bsidx]      
   
   #Fit SVM to the bootstrap sample
   svm.model<-svm(bsdat, bslbl, scale=FALSE, type=svmtyp, kernel=kern, cost=pparams, fitted=TRUE, probability=TRUE)

   #Get Platt scores from the whole data set
   pred<-predict(svm.model,dat,probability=TRUE)
   iter.prob.mat<-attr(pred, "probabilities")[,]

   #build up vector of known match (non-null) scores:
   iter.score.km<-iter.prob.mat[cbind(1:nrow(iter.prob.mat),as.numeric(lbls))]
   nonnull.vec<-c(nonnull.vec,as.numeric(iter.score.km)) #here, there is only one km score per observation

   #build up vector of known non-match (null) scores:
   iter.score.knm<-NULL
   for(j in 1:nrow(iter.prob.mat))
    {
     tmp.prob.vec<-iter.prob.mat[j,]
     tmp.prob.vec<-tmp.prob.vec[-as.numeric(lbls[j])]  #here, there are k-1 knm scores per observation. Take mean/median/max, something else  here??
     iter.score.knm<-c(iter.score.knm,tmp.prob.vec)
    }
   iter.score.knm<-as.numeric(iter.score.knm)
   null.vec<-c(null.vec,iter.score.knm)
   
   print(paste("B.S. Iter: ",i))
  }

#print(warnings())

if(printQ==TRUE)
 {
#   hist(null.vec)
#   par(new=TRUE)
#   hist(nonnull.vec)
#   hist(c(null.vec,nonnull.vec)) 
  hist(null.vec,probability=TRUE,xlim=c(0,1),ylim=c(0,ylim.max),col=rgb(1,0,0,1/4),main="",xlab="")
  par(new=TRUE)
  hist(nonnull.vec,probability=TRUE,xlim=c(0,1),ylim=c(0,ylim.max),col=rgb(0,1,0,1/4),main="Null(KNM, red) and Non-null(KM, green)",xlab="Platt-score")
  #hist(c(null.vec,nonnull.vec))
 }

print(paste("Maximum null score    : ",max(null.vec)))
print(paste("Minimum non-null score: ",min(nonnull.vec)))

return(list(null.vec,nonnull.vec))

}


#------------------------------------------------------------------------------------
#Variant of above. Estimate svm Platt-score null and non-null distributions via a group-wise bootstrap
#
#**MAJOR CHANGES FROM bs.null above:
#**This time "flip a coin" to see which KNM score we keep for an observation
#**ALSO, only keep scores from observations NOT in the BS sample.  
#
#Philosophy is that this decreases depencence between scores and also makes the null distribution of scores a bit more 
#conservative. Will need to run more BS iterations however.
#
#Adapted from Storey and Tibshirani permutation method in PNAS
#------------------------------------------------------------------------------------
bs.null2<-function(dat,lbls,nbs,svmtyp="C-classification",kern="linear",pparams=0.1,printQ=FALSE,ylim.max=NULL)
{
  #parse lbls into groups:
  lbls.idxs<-lapply(1:nlevels(lbls),function(x){which(lbls==levels(lbls)[x])})
  
  #Initialize score arrays this time instead of iteratively building them up. Make make things rin faster
#   null.vec<-rep(0.0,nbs*nrow(dat))    #To keep track of KM (non-null scores)
#   nonnull.vec<-rep(0.0,nbs*nrow(dat)) #Save only 1 random KNM score per obs-iteration
  null.vec<-NULL
  nonnull.vec<-NULL
  for(i in 1:nbs)
  {
    #Grab a bootstrap sample, bootstrapping out of each GROUP invididually:
    bsidx<-unlist(as.vector(sapply(1:length(lbls.idxs),function(x){sample(lbls.idxs[[x]],length(lbls.idxs[[x]]), replace=TRUE )})))
    #   print(bsidx)
    bsdat<-dat[bsidx,]
    
    rownames(bsdat)<-NULL #Shuts off an annoying warning
    bslbl<-lbls[bsidx]      
    
    #Fit SVM to the bootstrap sample
    svm.model<-svm(bsdat, bslbl, scale=FALSE, type=svmtyp, kernel=kern, cost=pparams, fitted=TRUE, probability=TRUE)
    
    #Get Platt scores from the whole data set
    pred<-predict(svm.model,dat,probability=TRUE)
    mod.bsidx<-unique(sort(bsidx))
    iter.prob.mat<-attr(pred, "probabilities")[-mod.bsidx,] #Keep only scores not contained in the BS sample
    #rownames(iter.prob.mat)<-NULL
    #print(paste("Iter:",i))
    #print(iter.prob.mat)
    
    #build up vector of known match (non-null) scores:
    iter.score.km<-iter.prob.mat[ cbind(1:nrow(iter.prob.mat),as.numeric(lbls[-mod.bsidx])) ]
    #print(length(iter.score.km))
    nonnull.vec<-c(nonnull.vec,as.numeric(iter.score.km)) #here, there is only one km score per observation
    
    #build up vector of known non-match (null) scores:
    iter.score.knm<-NULL
    for(j in 1:nrow(iter.prob.mat))
    {
      tmp.prob.vec<-iter.prob.mat[j,]
      tmp.prob.vec<-sample(tmp.prob.vec[-as.numeric(lbls[-mod.bsidx][j])],1)  #Just randomly grab one of the KNM scores
      iter.score.knm<-c(iter.score.knm,tmp.prob.vec)
    }
    iter.score.knm<-as.numeric(iter.score.knm)
    null.vec<-c(null.vec,iter.score.knm)
    
    print(paste("B.S. Iter: ",i))
   }
  
  #print(warnings())
  
  if(printQ==TRUE)
  {
    #   hist(null.vec)
    #   par(new=TRUE)
    #   hist(nonnull.vec)
    #   hist(c(null.vec,nonnull.vec)) 
    hist(null.vec,probability=TRUE,xlim=c(0,1),ylim=c(0,ylim.max),col=rgb(1,0,0,1/4),main="",xlab="")
    par(new=TRUE)
    hist(nonnull.vec,probability=TRUE,xlim=c(0,1),ylim=c(0,ylim.max),col=rgb(0,1,0,1/4),main="Null(KNM, red) and Non-null(KM, green)",xlab="Platt-score")
    #hist(c(null.vec,nonnull.vec))
  }
  
  print(paste("Maximum null score    : ",max(null.vec)))
  print(paste("Minimum non-null score: ",min(nonnull.vec)))
  
  return(list(null.vec,nonnull.vec))
  
}


#--------------------------------------------------------
#Split a set of scores into null and non-null sets
#--------------------------------------------------------
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


#--------------------------------------------------------
#Plot histograms of output scores
#--------------------------------------------------------
scores.histograms<-function(null.vec, nonnull.vec, ylim.max) {
  
  hist(null.vec,probability=TRUE,xlim=c(0,1),ylim=c(0,ylim.max),col=rgb(1,0,0,1/4),main="",xlab="")
  par(new=TRUE)
  hist(nonnull.vec,probability=TRUE,xlim=c(0,1),ylim=c(0,ylim.max),col=rgb(0,1,0,1/4),main="Null(KNM, red) and Non-null(KM, green)",xlab="Platt-score")
  
}