library(fdrID)

dat.km <- read.table("/Users/npetraco/latex/papers/posterior_error/comparison_examples/Jacqueline/KMs_from_QAs_use_me_jacqueline.txt")
dat.knm <- read.table("/Users/npetraco/latex/papers/posterior_error/comparison_examples/Jacqueline/KNMs_use_me_jacqueline.txt")

hist(dat.knm[,3])
hist(dat.km[,3])
hist(c(dat.km[,3],dat.knm[,3]),bre=120)


hist(log(dat.knm[,3]))
hist(log(dat.km[,3]))
hist(c( log(dat.knm[,3]), log(dat.km[,3]) ))


mn.null.vec <- mean(log(dat.knm[,3]))
sd.null.vec <- sd(log(dat.knm[,3]))
mn.null.vec
sd.null.vec

null.vec <- log(dat.knm[,3])
null.vec.std <- (null.vec - mn.null.vec)/sd.null.vec
hist(null.vec.std)

nonnull.vec <- log(dat.km[,3])
hist(nonnull.vec)

nonnull.vec.std <- (nonnull.vec - mn.null.vec)/sd.null.vec
hist(nonnull.vec.std)

scores.histograms(null.vec, nonnull.vec, main = "Null-KNM(red) and Non-null-KM(green)",xlab="log(MPOC Score)")
scores.histograms(null.vec.std, nonnull.vec.std, xlim.min = -3, xlim.max = 8, main = "Null-KNM(red) and Non-null-KM(green)",xlab="standardized log(MPOC Score)")
scores.histograms(dat.knm[,3], dat.km[,3], main = "Null-KNM(red) and Non-null-KM(green)",xlab="MPOC Score") #Untransformed scores
max(null.vec.std)
min(nonnull.vec.std)

#Fit null (log(KNM_scores))
sn.null.fit.info <- null.sn.fit(null.vec, standardizeQ=T, plotQ=T)
sn.null.fit.info

gev.null.fit.info <- null.gev.fit(null.vec, standardizeQ=T, plotQ=T)
gev.null.fit.info

nig.null.fit.info <- null.nig.fit(null.vec, standardizeQ=T, plotQ=T)
nig.null.fit.info

lg.null.fit.info <- null.lg.fit(null.vec, standardizeQ=T, plotQ=T)
lg.null.fit.info


#preliminary p-values check
sn.null.fit.info$parameters
p.null <- 1 - psn(null.vec.std, xi=sn.null.fit.info[[1]][1], omega=sn.null.fit.info[[1]][2], alpha=sn.null.fit.info[[1]][3])
hist(p.null)
p.nonnull <- 1 - psn(nonnull.vec.std, xi=sn.null.fit.info[[1]][1], omega=sn.null.fit.info[[1]][2], alpha=sn.null.fit.info[[1]][3])
hist(p.nonnull, xlab="p-value", main="Non-null(KM) p-values")
hist(c(p.null,p.nonnull))

hist(qnorm(c(p.null)))
hist(qnorm(c(p.nonnull)))
hist(qnorm(c(p.null,p.nonnull)))

length(nonnull.vec)

length(null.vec)
n.sims <- 19900 #Approx num KNM 
null.vec.val.std <- rsn(n.sims, xi=sn.null.fit.info[[1]][1], omega=sn.null.fit.info[[1]][2], alpha=sn.null.fit.info[[1]][3])
hist(null.vec.val.std, main="IID random sample of modeled Null scores", xlab="standardized log(MPOC Score)")
p.null.val <- 1 - psn(null.vec.val.std, xi=sn.null.fit.info[[1]][1], omega=sn.null.fit.info[[1]][2], alpha=sn.null.fit.info[[1]][3])
hist(p.null.val)
hist(qnorm(p.null.val), xlab="z-value", main=expression(paste(Phi^-1, "(p_KNM)")))
hist(qnorm(c(p.null.val,p.nonnull)), xlab="z-value", main="f(z)")


zsc.prelim <- qnorm(c(p.null.val,p.nonnull))
p.vals <- c(p.null.val,p.nonnull)

#First Do Efron lfdr fit to examine fit to f(z)
library(locfdr)
fdr.model<-locfdr(qnorm(p.vals), bre = 120, df = 8, pct = 0, pct0 = 1/4, nulltype = 1, type =0, plot = 1, main = " ", sw = 0)

minfo <- sampler.prep(p.vals, num.bins=120, degree=8, interceptQ=T, overdispersionQ=F, sampler="jags")
jsim <- jags(data=minfo$Data, inits=minfo$Initialization.Function, minfo$Model.Parameters, n.iter=10000, n.chains=4, model.file=minfo$BUG.Model.File.Path)
jsim


posterior.f <- jsim$BUGSoutput$sims.list$lambda
x <- minfo$Bin.Midpoints
z <- minfo$z.Values

lfdr.info <- make.fdr.functions(z, x, pct0=0.25, posterior.f, credibility.level=0.95, interval.type="hpd", plotQ=T)
efdrf <- lfdr.info$fdr.means.func
med.fdrf <- lfdr.info$fdr.means.func
uci.fdrf <- lfdr.info$upper.fdr.ci.func
lci.fdrf <- lfdr.info$lower.fdr.ci.func
etdrf <- lfdr.info$tdr.means.func
med.tdrf <- lfdr.info$tdr.means.func
uci.tdrf <- lfdr.info$upper.tdr.ci.func
lci.tdrf <- lfdr.info$lower.tdr.ci.func

points(seq(min(z),max(x), length.out = 119), posterior.f[1,]/max(posterior.f[1,])*100)

p0 <- lfdr.info$p0
d0 <- lfdr.info$delta0
s0 <- lfdr.info$sigma0

hist(p0)
mean(p0)

z.nonnull <- qnorm(p.nonnull)
z.null <- qnorm(p.null)
z.unknowns <- z.nonnull

zoomed.post.prob.plot(c(z.null,z.nonnull), zbounds=c(NA, max(z.unknowns)), 
                      point.est.func=efdrf, 
                      upper.est.func=uci.fdrf, 
                      lower.est.func=lci.fdrf, prob.scale="percent", xlab=NULL, ylab="Pr(not source | z)   (%)", main=NULL)

min(z.unknowns)
min(c(z.null,z.nonnull))
max(z.unknowns)
max(c(z.null,z.nonnull))

points(z.unknowns[which(100*efdrf(z.unknowns) <=5)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) <=5)],col="green")
points(z.unknowns[which(100*efdrf(z.unknowns) > 5 & 100*efdrf(z.unknowns) <= 50)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) > 5 & 100*efdrf(z.unknowns) <= 50)],col="yellow")
points(z.unknowns[which(100*efdrf(z.unknowns) > 50)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) > 50)],col="red")

#----------------------------------------------------------------------------------
zoomed.post.prob.plot(c(z.null,z.nonnull), zbounds=c(NA, max(z.unknowns)), prbounds=c(1, NA),
                      point.est.func=etdrf, 
                      upper.est.func=uci.tdrf, 
                      lower.est.func=lci.tdrf, prob.scale="percent", xlab="MPOC z-score", ylab="", main="MPOC KM Scores Projected on the Model")
title(ylab=expression(paste(hat(Pr), "(match | z)   (%)")),mgp=c(2.5,1,0))

points(z.unknowns[which(100*etdrf(z.unknowns) >=95)], 100*etdrf(z.unknowns)[which(100*etdrf(z.unknowns) >=95)],col="green",pch=16)
points(z.unknowns[which(100*etdrf(z.unknowns) < 95 & 100*etdrf(z.unknowns) >= 50)], 100*etdrf(z.unknowns)[which(100*etdrf(z.unknowns) < 95 & 100*etdrf(z.unknowns) >= 50)],col="yellow",pch=16)
points(z.unknowns[which(100*etdrf(z.unknowns) < 50)], 100*etdrf(z.unknowns)[which(100*etdrf(z.unknowns) < 50)],col="red",pch=16)

cbind(
z.unknowns[which(100*etdrf(z.unknowns) >=95)], 
100*uci.tdrf(z.unknowns)[which(100*etdrf(z.unknowns) >=95)],
100*etdrf(z.unknowns)[which(100*etdrf(z.unknowns) >=95)],
100*lci.tdrf(z.unknowns)[which(100*etdrf(z.unknowns) >=95)]
)

ord.idx <- order(dat.km[,3])

res.prelim <- cbind(
  dat.km[ord.idx,3],
  z.unknowns[ord.idx], 
  round(100*uci.tdrf(z.unknowns)[ord.idx],4),
  round(100*etdrf(z.unknowns)[ord.idx],4),
  round(100*lci.tdrf(z.unknowns)[ord.idx],4)
)

colnames(res.prelim) <- c("MPOC_KM","z_KM", "Lower 95% CI", "Pr(match|z)", "Upper 95% CI")
res.prelim
plot(res.prelim[,1], res.prelim[,4], xlab="MPOC score", ylab=expression(paste(hat(Pr), "(match|z)")) , mgp = c(2.5, 1, 0), typ="l", main="MPOC Scores vs. Non-null (match) Posterior Probability Est.")
points(res.prelim[which(res.prelim[,4] >=95),1], res.prelim[which(res.prelim[,4] >=95), 4],col="green",pch=16)
points(res.prelim[which(res.prelim[,4] < 95 & res.prelim[,4] >= 50), 1], res.prelim[which(res.prelim[,4] < 95 & res.prelim[,4] >= 50), 4],col="yellow",pch=16)
points(res.prelim[which(res.prelim[,4] < 50),1], res.prelim[which(res.prelim[,4] < 50),4],col="red",pch=16)


#-----------------------------------------------------

zoomed.post.prob.plot(c(z.null,z.nonnull), zbounds=c(NA, NA), 
                      point.est.func=etdrf, 
                      upper.est.func=uci.tdrf, 
                      lower.est.func=lci.tdrf, prob.scale="percent", xlab=NULL, ylab="Pr(match | z)   (%)", main="")
points(seq(min(z),max(x), length.out = 119), posterior.f[1,]/max(posterior.f[1,])*100)
