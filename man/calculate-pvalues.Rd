\name{calculate-pvalues}
\alias{calculate.sample.pvalue}
\alias{calculate.ratio.pvalue}
\alias{calculate.mult.sample.pvalue}
\alias{adjust.ratio.pvalue}
\title{
  Calculate and Adjust Ratio and Sample p-values.
}

\description{
  Functions for calculating and adjusting ratios and sample p-values.
  Usually, these are called by proteinRatios or peptideRatios.
}

\usage{
  calculate.ratio.pvalue(lratio, variance, ratiodistr = NULL)
  calculate.sample.pvalue(lratio, ratiodistr)

  calculate.mult.sample.pvalue(lratio, ratiodistr, strict.pval, 
                               lower.tail, n.possible.val, n.observed.val)

  adjust.ratio.pvalue(quant.tbl, p.adjust, sign.level, globally = FALSE)
}

\arguments{
  \item{lratio}{log 10 protein or peptide ratios.}
  \item{ratiodistr}{Fitted ratio distribution/}
  \item{variance}{Variance of lratios.}

  \item{strict.pval}{If FALSE, missing ratios are ignored. If TRUE, missing ratios are penalized by giving them a sample.pval of 0.5.}
  \item{lower.tail}{lower.tail of distribution?}
  \item{n.possible.val}{Number of possible ratios.}
  \item{n.observed.val}{Number of observed ratios.}

  \item{quant.tbl}{Quantification table (from proteinRatios or peptideRatios).}
  \item{p.adjust}{p-value adjustment method (see ?p.adjust).}
  \item{sign.level}{Ratio significance level.}
  \item{globally}{Whether the p-values should be adjusted over all conditions,
                  or individually in each condition.}

}
\author{Florian P. Breitwieser}

\seealso{
 \code{\link{proteinRatios}},\code{\link{peptideRatios}}
}

\examples{
  lratio <- c(-1,-1,seq(from=-1,to=1,by=.25),1,1)
  variance <- c(0,1,rep(0.1,9),0,1)
  ratiodistr.precise <- new("Norm",mean=0,sd=.25)
  ratiodistr.wide <- new("Norm",mean=0,sd=.5)

  # ratio p-value is impacted only by the variance
  # sample p-value captures whether the ratio distribution is narrow ('precise')
  #  or wide
  data.frame(lratio, variance,
             ratio.pvalue=calculate.ratio.pvalue(lratio, variance),
             sample.pvalue.precise=calculate.sample.pvalue(lratio,ratiodistr.precise),
             sample.pvalue.wide=calculate.sample.pvalue(lratio,ratiodistr.wide))

}
