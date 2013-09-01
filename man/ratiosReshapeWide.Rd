\name{ratiosReshapeWide}
\alias{ratiosReshapeWide}

\title{Reshape output of proteinRatios into wide format}
\description{
  Reshape output of proteinRatios into wide format
}

\usage{
ratiosReshapeWide(quant.tbl, vs.class = NULL, sep = ".", cmbn = NULL,
                  short.names = FALSE)
}
\arguments{
  \item{quant.tbl}{Output of proteinRatios or peptideRatios.}
  \item{vs.class}{Only return ratios where class1 is vs.class}
  \item{sep}{Separator for column names in the reshape.}
  \item{cmbn}{Not functional.}
  \item{short.names}{If vs.class is set and short.names=TRUE, then the comparision 
    name will be i.e. 'class2' instead of 'class2/class1'.}
}
\author{Florian P. Breitwieser}
