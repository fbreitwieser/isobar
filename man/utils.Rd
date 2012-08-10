\name{Isobar util functions}

\alias{paste0}
\alias{\%inrange\%}

\title{
  Isobar util functions
}

\description{
  Utility functions.
  paste0 as a shorthand to paste(...,sep="") in versions of R pre 2.14.
}

\usage{
paste0(..., sep = "")
a \%inrange\% b
}

\arguments{
  \item{\dots}{Arguments to paste.}
  \item{sep}{Separator.}
  \item{a}{values.}
  \item{b}{range.}
}

\author{
  Florian P Breitwieser
}

\examples{
  1:10 %inrange% c(2,5)
}

