\name{TlsParameter-class}
\Rdversion{1.1}
\docType{class}
\alias{TlsParameter-class}

\title{Class \code{"TlsParameter"}}
\description{The parameter of a location scale t distribution, used by Tlsd-class
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TlsParameter", ...)}.
Usually an object of this class is not needed on its own, it is generated automatically when an object of the class
\code{Tlsd} is instantiated. 
}
\section{Slots}{
  \describe{
    \item{\code{df}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{location}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{scale}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Parameter}"}, directly.
Class \code{"\linkS4class{OptionalParameter}"}, by class "Parameter", distance 2.
}
\section{Methods}{
No methods defined with class "TlsParameter" in the signature.
}
\author{
  Florian P. Breitwieser, based on original TParameter class.
}
\seealso{
\code{\linkS4class{Tlsd}}
}
\examples{
showClass("TlsParameter")
}
\keyword{classes}
