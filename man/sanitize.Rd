\name{sanitize}

\alias{sanitize}

\title{
        Helper function for LaTeX export
}
\description{
        Sanitizes strings for LaTeX
}

\usage{
    sanitize(str, dash = TRUE)
}
\arguments{
  \item{str}{character string to be escaped}
  \item{dash}{shoud a dash ('-') should be escaped to a '\\nobreakdash-'?}
}
\value{
    escaped character
}
\author{iQuantitator,Florian P Breitwieser
}
\examples{
    sanitize("\\textbf{123-123}")
}
