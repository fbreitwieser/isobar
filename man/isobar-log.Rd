\name{IBSpectra.log}
\alias{do.log}
\alias{do.log,IBSpectra,character-method}
\alias{get.log}
\alias{get.log,IBSpectra,character-method}
\alias{is.logged}
\alias{is.logged,IBSpectra,character-method}
\title{
    Log functions for IBSpectra objects
}
\description{
    The slot log of IBSpectra objects contains a matrix with
    two columns which contain a timestamp and message. Rownames
    relate to the item logged. 

    For now, correctIsotopeImpurities and normalize use and write
    to the log.
}
\usage{
do.log(x, name, msg)

get.log(x, name)

is.logged(x, name)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{IBSpectra object}
  \item{name}{Name of property to be logged (translates to row name).}
  \item{msg}{Message to be logged for name.}
}
\details{
    A warning message will be displayed if a already logged property
    is logged again.
}
\value{
    IBSpectra object with updated log.
}
\author{Florian P Breitwieser}

\seealso{IBSpectra-class}
\examples{
    data(ibspiked_set1)
    ib <- normalize(correctIsotopeImpurities(ibspiked_set1))
    ib@log
}
