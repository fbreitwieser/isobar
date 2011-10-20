\name{number.ranges}
\alias{number.ranges}
\title{
    Helper function to transform number lists to ranges
}
\description{
    1,2,3,4,5,8,9,10 -> 1-5,8-10
}

\usage{
number.ranges(numbers)
}
\arguments{
  \item{numbers}{numeric}
}
\value{character}
\examples{
    number.ranges(c(1,2,3,9,3,10,8,11))
}
\author{Florian P Breitwieser}
