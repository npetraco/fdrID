\name{fdrID-package}
\alias{fdrID-package}
\alias{fdrID}
\docType{package}
\title{
Computes local false discovery rate probabilities using JAGS, Stan and code from the locfdr package
}
\description{
Using ideas from Ann. Stat. 35(4):1351-1377 (2007), and some handy code from the former locfdr package, this package computes local false discovery rate (i.e. posterior error probability) estimates for classification/discrimination tasks using a Bayesian Poission regression. Specifically, the Bayesian Poission regression is used to estimate the denominator of Bayes Theorem within Efron's two-groups empirical Bayes methodology. The regression is a heirarical model and carried out with either JAGS (mcmc-jags.sourceforge.net) or Stan (mc-stan.org), via their R interfaces.
}
\details{
\tabular{ll}{
Package: \tab fdrID\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-05-25\cr
License: \tab GPL >=2 \cr
}
In Ann. Stat. 35(4):1351-1377 (2007), Efron describes fitting methods for his two-groups empirical Bayes' methodology. This methodology is used to compute local false discovery rates, also known as posterior error probabilities (Kall 2008). The heart of the method requires that a histogram of z-score counts, representing the denominator in Bayes' theorem, is fit via poission regression. Efron uses a maximum likelihood fit. This package uses a Bayesian hierarchical fit instead to ease the inclusion of overdispersion terms in the event that the histogram counts are correlated. The user supplys vectors of (null/non-null) p-values or z-scores. In the event that the user needs to generate such values for a classification/discrimination task, they can use the (forthcomming) svm-zscore package.
}
\author{
Maintainer: Nick Petraco <npetraco@gmail.com>
}
\references{
Efron B. Size Power and False Discovery Rates. Ann Stat 35(4):1351-1377 (2007)

Efron B, Turnbull B, Narasimhan B. locfdr. http://cran.r-project.org/src/contrib/Archive/locfdr/

Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)

Kall L, Storey JD, MacCross MJ, Noble WS. Posterior Error Probabilities and False Discovery Rates: Two Sides of the Same Coin. J Proteome Res 7:40-44 (2008)
}
\keyword{ package }