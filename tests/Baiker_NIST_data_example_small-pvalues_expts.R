library(svmSCORES)
library(fdrID)
library(caret)

Ztr <- as.matrix(read.csv("/Users/npetraco/latex/papers/posterior_error/comparison_examples/NIST_Martin_screwdriver/NFIToolmark_Biasotti_Murdock_features_15-30_Ztr.csv",header=T)[,2:199])
Zval <- as.matrix(read.csv("/Users/npetraco/latex/papers/posterior_error/comparison_examples/NIST_Martin_screwdriver/NFIToolmark_Biasotti_Murdock_features_15-30_Zval.csv",header=T)[,2:199])
Zte <- as.matrix(read.csv("/Users/npetraco/latex/papers/posterior_error/comparison_examples/NIST_Martin_screwdriver/NFIToolmark_Biasotti_Murdock_features_15-30_Zte.csv",header=T)[,2:199])

label.mat <- read.csv("/Users/npetraco/latex/papers/posterior_error/comparison_examples/NIST_Martin_screwdriver/NFIToolmark_Biasotti_Murdock_features_15-30_SET-LABELS.csv",header=T)[,2:4]
lbltr <- as.factor(label.mat[,1])
lblval <- as.factor(label.mat[,2])
lblte <- as.factor(label.mat[,3])

#Get a set of bootstrapped null Platt scores:
pscs <- bootstrap.platt.scores.parallel(
  num.processes=8, 
  dat.mat=Ztr, 
  lbls=lbltr, 
  nbs=2000, 
  svmtyp="C-classification", 
  kern="linear", 
  pparams=0.1,
  timerQ=T)

#Before looking at straight empirical determination of p-values, try this first:
zscs.info <- 
  zscore.fit2(pscs, Ztr, Zval, lbltr, lblval, 
              pvalue.method = "integral",
              distribution="nig",
              num.processes=8, 
              standardizeQ=T, 
              num.bs.iter=2000, 
              C.param = 0.1, 
              printQ=F, plotQ=T)


xlms <- c(
floor(min(
c(min(zscs.info$score.null.validation),
min(zscs.info$score.null.training),
min(zscs.info$boostrapped.null.platt.scores),
min(zscs.info$score.nonnull.validation))
)),
ceiling(max(
c(max(zscs.info$score.null.validation),
max(zscs.info$score.null.training),
max(zscs.info$boostrapped.null.platt.scores),
max(zscs.info$score.nonnull.validation))
)))
xlms

xx <- seq(from=xlms[1], to=xlms[2],by=0.1)
pmv <- zscs.info$fit.info.and.diagnostics$parameters
pmv
plot(xx,fBasics::dnig(xx, alpha=pmv[1], beta=pmv[2], delta=pmv[3], mu=pmv[4]),typ="l",lwd=3, ylab="")
par(new=T)
hist(zscs.info$score.nonnull.validation, xlim=xlms, axes=F, ylab="", xlab="", main="", col="green")
par(new=T)
hist(zscs.info$score.null.validation, xlim=xlms, axes=F, ylab="", xlab="", main="")
par(new=T)
hist(zscs.info$score.null.training, xlim=xlms, axes=F, ylab="", xlab="", main="")
par(new=T)
hist(zscs.info$boostrapped.null.platt.scores, xlim=xlms, axes=F, ylab="", xlab="", main="")

dat <- cbind(xx,fBasics::dnig(xx, alpha=pmv[1], beta=pmv[2], delta=pmv[3], mu=pmv[4]),
             1-fBasics::pnig(xx, alpha=pmv[1], beta=pmv[2], delta=pmv[3], mu=pmv[4]))
plot(dat[117:131,1],1e9*dat[117:131,2], main="Null NIG Right Tail")
plot(dat[117:131,1],log(1e9*dat[117:131,2]))

alp <- pmv[1]
bet <- pmv[2]
delt <- pmv[3]
mu <- pmv[4]


xarg <- dat[,1]
fpdf <- (alp*delt*besselK(alp*sqrt(delt^2 + (xarg - mu)^2), 1, expon.scaled = FALSE))/(pi * sqrt(delt^2 + (xarg - mu)^2)) * exp(delt*sqrt(alp^2-bet^2) + bet*(xarg-mu))
cbind(dat[,2],fpdf)

fit <- lm(log(1e9*dat[117:131,2]) ~ dat[117:131,1] + I(dat[117:131,1]^2))
summary(fit)
#abline(fit)
qqnorm(residuals(fit))
qqline(residuals(fit))


fBasics::dnig(lgs, alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est)

zscs.info["Were.nonnull.zvalues.smeared?"]

z.nonnull.pot <- smear.extreme.nonnull.zvalues(
  zscs.info$nonnull.pvalues,
  upper.set.zvalue = (-6.5), 
  mu.factor = (-0.5), 
  p.factor=0.95, plotQ=T)[[3]]
#z.nonnull <- z.nonnull.pot


#Checks:
names(zscs.info)
zscs.info$fit.info.and.diagnostics
z.null <- zscs.info$null.zvalues
z.nonnull <- zscs.info$nonnull.zvalues.smeared
#p.nonnull <- zscs.info$nonnull.pvalues
p.nonnull <- pnorm(z.nonnull)

cbind(p.nonnull,z.nonnull)
check.ps.and.zs(null.p.values = pnorm(z.null), nonnull.p.values = pnorm(z.nonnull), printQ = T, plotQ = T)
hist(z.nonnull, main="Validation z Non-Null")
hist(z.null, main="Validation z Null")
hist(c(z.null,z.nonnull), main="All z")

hist(zscs.info$nonnull.pvalues)

zscs.info["Were.nonnull.zvalues.smeared?"]

#Fit lfdr models:
p.vals <- c(pnorm(z.null), pnorm(z.nonnull))

#First Do Efron lfdr fit to examine fit to f(z)
library(locfdr)
fdr.model<-locfdr(qnorm(p.vals), bre = 120, df = 7, pct = 0, pct0 = 1/4, nulltype = 1, type =0, plot = 1, main = " ", sw = 0)

minfo <- sampler.prep(p.vals, num.bins=120, degree=7, interceptQ=T, overdispersionQ=F, sampler="jags")
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
mean(d0)
mean(s0)

#-------------------------
# Test input Zte and output Platts
mln <- mean(log(pscs[,1]))
sln <- sd(log(pscs[,1]))
tinfo3 <- posterior.probs.for.unks2(training.dmat=Ztr, training.labels=lbltr, C.param=0.1, test.dmat = Zte,
                                    standardizeQ=zscs.info$standardization.flag, 
                                    mean.log.null.score=mln, sd.log.null.score=sln,
                                    pvalue.method = "integral",
                                    null.vec.training = zscs.info$score.null.training,
                                    distribution.name = zscs.info$fit.distribution.name,
                                    distribution.fit.info = zscs.info$fit.info.and.diagnostics,
                                    pp.point.est.func=efdrf,
                                    pp.uci.est.func=uci.fdrf,
                                    pp.lci.est.func=lci.fdrf)
tinfo3
plot(tinfo3[,4],typ="h")

z.unknowns <- tinfo3[,2]
z.unknowns
zoomed.post.prob.plot(c(z.null,z.nonnull), zbounds=c(NA, max(z.unknowns)), 
                      point.est.func=efdrf, 
                      upper.est.func=uci.fdrf, 
                      lower.est.func=lci.fdrf, prob.scale="percent", xlab=NULL, ylab="fdr (%)", main=NULL)
min(z.unknowns)
min(c(z.null,z.nonnull))
max(z.unknowns)
max(c(z.null,z.nonnull))

points(z.unknowns[which(100*efdrf(z.unknowns) <=5)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) <=5)],col="green")
points(z.unknowns[which(100*efdrf(z.unknowns) > 5 & 100*efdrf(z.unknowns) <= 50)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) > 5 & 100*efdrf(z.unknowns) <= 50)],col="yellow",pch=16)
points(z.unknowns[which(100*efdrf(z.unknowns) > 50)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) > 50)],col="red")
#Got these one wrong
wrong.idxs <- which((tinfo3[,1] == lblte)==FALSE)
data.frame(lblte[wrong.idxs],tinfo3[wrong.idxs,])
points(z.unknowns[wrong.idxs], 100*efdrf(z.unknowns)[wrong.idxs],col="black", pch=16)

#Error rate est from test set:
length(wrong.idxs)/nrow(Zte) * 100
