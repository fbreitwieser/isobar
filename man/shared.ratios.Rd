\name{shared.ratios}
\alias{shared.ratios}
\title{
  Shared ratio calculation
}
\description{
  Calculate ratios of reporter proteins and subset proteins with shared peptides.
}
\usage{
shared.ratios(ibspectra, noise.model, channel1 , channel2 , protein = reporterProteins(proteinGroup(ibspectra)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ibspectra}{IBspectra object.}
  \item{noise.model}{NoiseModel object.}
  \item{channel1}{channel1 to compare.}
  \item{channel2}{channel2 to compare.}
  \item{protein}{proteins for which the calculation should be made.}
  \item{...}{Additional arguments passed to estimteRatio.}
}
\value{
  data.frame
}
\author{
  Florian P.\ Breitwieser
}

\seealso{
  \code{\link{shared.ratios.sign}}
}
