\name{Tlsd-class}
\Rdversion{1.1}
\docType{class}
\alias{Tlsd-class}
\alias{Tlsd}

\title{Class \code{"Tlsd"}}
\description{
  Location scale family T distribution, based on the original T function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Tlsd", df, location, scale)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{gaps}:}{Object of class \code{"OptionalMatrix"} ~~ }
    \item{\code{img}:}{Object of class \code{"rSpace"} ~~ }
    \item{\code{param}:}{Object of class \code{"OptionalParameter"} ~~ }
    \item{\code{r}:}{Object of class \code{"function"} ~~ }
    \item{\code{d}:}{Object of class \code{"OptionalFunction"} ~~ }
    \item{\code{p}:}{Object of class \code{"OptionalFunction"} ~~ }
    \item{\code{q}:}{Object of class \code{"OptionalFunction"} ~~ }
    \item{\code{.withSim}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.withArith}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.logExact}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.lowerExact}:}{Object of class \code{"logical"} ~~ }
    \item{\code{Symmetry}:}{Object of class \code{"DistributionSymmetry"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{AbscontDistribution}"}, directly.
Class \code{"\linkS4class{UnivariateDistribution}"}, by class "AbscontDistribution", distance 2.
Class \code{"\linkS4class{AcDcLcDistribution}"}, by class "AbscontDistribution", distance 2.
Class \code{"\linkS4class{Distribution}"}, by class "AbscontDistribution", distance 3.
Class \code{"\linkS4class{UnivDistrListOrDistribution}"}, by class "AbscontDistribution", distance 3.
}
\section{Methods}{
No methods defined with class "Tlsd" in the signature.
}

\author{
  Florian P. Breitwieser, based on original T distribution class.
}

\examples{
showClass("Tlsd")
}
\keyword{classes}
