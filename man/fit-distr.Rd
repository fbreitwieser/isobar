\name{fit distributions}
\alias{fitCauchy}
\alias{fitTlsd}
\alias{fitNorm}
\alias{fitWeightedNorm}
\alias{fitNormalCauchyMixture}
\alias{fitGaussianMixture}
\title{Fit weighted and unweighted Cauchy and Normal distributions}
\description{
  Functions to fit the probability density functions
    on ratio distribution.
}
\usage{
fitCauchy(x)
fitNorm(x, portion = 0.75)
fitWeightedNorm(x, weights)
fitNormalCauchyMixture(x)
fitGaussianMixture(x, n = 500)
fitTlsd(x)
}
\arguments{
  \item{x}{Ratios}
  \item{weights}{Weights}
  \item{portion}{Central portion of data to take for computation}
  \item{n}{number of sampling steps}
}
\value{
  \code{\link[distr]{Cauchy}},\code{\link[distr]{Norm}}
}
\author{
  Florian P Breitwieser, Jacques Colinge.
}

\seealso{
  \link{proteinRatios}
}
\examples{
library(distr)
data(ibspiked_set1)
data(noise.model.hcd)
# calculate protein ratios of Trypsin and CERU_HUMAN. Note: this is only
#  for illustration purposes. For estimation of sample variability, data
#  from all protein should be used
pr <- proteinRatios(ibspiked_set1,noise.model=noise.model.hcd,
                    cl=as.character(c(1,1,2,2)),combn.method="intraclass",protein=c("136429","P00450"))

# fit a Cauchy distribution
ratiodistr <- fitCauchy(pr$lratio)
plot(ratiodistr)
}
