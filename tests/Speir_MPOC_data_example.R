library(fdrID)

dat.km <- read.table("/Users/npetraco/latex/papers/posterior_error/comparison_examples/KMs_from_QAs_use_me_jacqueline.txt")
dat.knm <- read.table("/Users/npetraco/latex/papers/posterior_error/comparison_examples/KNMs_use_me_jacqueline.txt")

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

scores.histograms(null.vec, nonnull.vec)
scores.histograms(null.vec.std, nonnull.vec.std, xlim.min = -3, xlim.max = 8)
max(null.vec.std)
min(nonnull.vec.std)

#Fit null (log(KNM_scores))
sn.null.fit.info <- null.sn.fit(null.vec, standardizeQ=T, plotQ=T)
sn.null.fit.info

gev.null.fit.info <- null.gev.fit(null.vec, standardizeQ=T, plotQ=T)
gev.null.fit.info

library(ghyp)
nig.null.fit.info <- null.nig.fit(null.vec, standardizeQ=T, plotQ=T)
nig.null.fit.info

lg.null.fit.info <- null.lg.fit(null.vec, standardizeQ=T, plotQ=T)
lg.null.fit.info


#p-values
sn.null.fit.info$parameters
p.null <- 1 - psn(null.vec.std, xi=sn.null.fit.info[[1]][1], omega=sn.null.fit.info[[1]][2], alpha=sn.null.fit.info[[1]][3])
hist(p.null)
p.nonnull <- 1 - psn(nonnull.vec.std, xi=sn.null.fit.info[[1]][1], omega=sn.null.fit.info[[1]][2], alpha=sn.null.fit.info[[1]][3])
hist(p.nonnull)
hist(c(p.null,p.nonnull))

hist(qnorm(c(p.null)))
hist(qnorm(c(p.nonnull)))
hist(qnorm(c(p.null,p.nonnull)))
