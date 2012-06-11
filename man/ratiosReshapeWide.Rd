\name{ratiosReshapeWide}
\alias{ratiosReshapeWide}

\title{Reshape output of proteinRatios into wide format}
\description{
  Reshape output of proteinRatios into wide format
}

\usage{
ratiosReshapeWide(quant.tbl, grouped.cols = TRUE, 
                  vs.class = NULL, sep = ".")
}
\arguments{
  \item{quant.tbl}{Output of proteinRatios or peptideRatios.}
  \item{grouped.cols}{Whether the columns should be grouped next to each other.}
  \item{vs.class}{Only return ratios where class1 is vs.class}
  \item{sep}{Separator for column names in the reshape.}
}
\author{Florian P. Breitwieser}
