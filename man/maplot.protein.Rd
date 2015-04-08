\name{maplot.protein}
\alias{maplot.protein}
\title{
    Ratio intensity plot for individual proteins
}
\description{
    Plots ratio-versus-intensity for a selected protein against a reference channel.
}
\usage{
maplot.protein(x, relative.to, protein, noise.model = NULL, channels = NULL, 
               xlim = NULL, ylim = NULL, identify = FALSE, add = FALSE, 
               pchs = NULL, log="xy", legend.pos = "topright", names = NULL, 
               legend.cex = 0.8, cols = pchs, ltys = 1, main = protein, 
               xlab = NULL, ylab = NULL, type="ma", show.lm = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{IBSpectra object}
  \item{relative.to}{
      a character vector specifying reporter tag names. Either of 
      length 1 or same length as channels.
  }
  \item{protein}{
      Protein group identifier.
  }
  \item{noise.model}{
      NoiseModel object.
  }
  \item{channels}{
       Reporter tag names.
  }
  \item{xlim}{
      See par.
  }
  \item{ylim}{
      See par.
  }
  \item{identify}{
      boolean. If \code{true}, \code{\link{identify}} is called with peptide labels.
  }
  \item{add}{
%%     ~~Describe \code{add} here~~
  }
  \item{pchs}{
    a vector of the same length as \code{channels}. See \code{pch} in
    \code{\link{plot.default}}.
  }
  \item{log}{
    a character string which contains \code{x} if the x axis is to
    be logarithmic, \code{y} if the y axis is to be logarithmic and
    \code{xy} or \code{yx} if both axes are to be logarithmic.
  }
  \item{legend.pos}{
    see \code{pos} in \code{\link{legend}}.
  }
  \item{names}{
    a character string of the same length as \code{channels}, legend text.
  }
  \item{legend.cex}{
    see \code{cex} in \code{\link{legend}}.
  }
  \item{cols}{
    a vector of the same length as \code{channels}. See \code{col} in
    \code{\link{plot.default}}.
  }
  \item{ltys}{
    a vector of the same length as \code{channels}. See \code{lty} in
    \code{\link{plot.default}}.
  }
  \item{main}{
    a main title for the plot  
  }
  \item{xlab}{
    a label for the x axis, defaults to a description of \code{x}.
  }
  \item{ylab}{
    a label for the y axis, defaults to a description of \code{y}.
  }
  \item{type}{type of plot}
  \item{\dots}{
    passed to \code{\link{plot}}.
  }
  \item{show.lm}{show LM}
}
\author{Florian P. Breitwieser}
