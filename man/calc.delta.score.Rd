\name{calc.delta.score}
\alias{calc.delta.score}
\alias{calc.pep.delta.score}
\alias{filterSpectraDeltaScore}

\title{Calculate Delta Score from Ion Score}
\usage{
calc.delta.score(my.data)
filterSpectraDeltaScore(my.data, min.delta.score=10, do.remove=FALSE)
}

\description{
  Calculates delta score from raw search engine score by substracting the best matching hit with the second best matching.
  data needs to have not only the best hit per spectrum, but multiple, to be able to calculate the delta score.
  filterSpectraDeltaScore calls \code{calc.delta.score} and filters spectra below a minum delta score.}

\arguments{
  \item{my.data}{IBSpectra data frame.}
  \item{min.delta.score}{Minimum delta score.}
  \item{do.remove}{If TRUE, spectra below the min.prob threshold are not just set as 'use.for.quant=FALSE' but removed.}
}
\value{
  Returns data with additional column 'delta.score'.
}
\author{
  Florian P. Breitwieser
}
