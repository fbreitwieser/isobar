\name{writeIBSpectra}
\alias{writeIBSpectra}
\title{
  Write IBSpectra file as CSV in a format readable by readIBSpectra.
}
\description{
  Write IBSpectra file using write.table with defaults in a format readable by readIBSpectra.
}
\usage{
writeIBSpectra(ibspectra, file, sep = "\t", row.names = FALSE, ...)
}
\arguments{
  \item{ibspectra}{IBSpectra object}
  \item{file}{file name.}
  \item{sep}{field separator string.}
  \item{row.names}{indicates whether row.names should be written.}
  \item{\dots}{further arguments to \code{\link{write.table}}}
}
\author{
  Florian P Breitwieser
}
