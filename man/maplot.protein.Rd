\name{maplot.protein}
\alias{maplot.protein}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    MAplot for individual proteins
}
\description{
    Plots ratio-versus-intensity for a selected protein against a reference channel.
}
\usage{
maplot.protein(x, relative.to, protein, noise.model = NULL, channels = NULL, ylim = c(0.01, 100), identify = FALSE, add = FALSE, pchs = NULL, legend.pos = "topright", names = NULL, legend.cex = 0.8, cols = pchs, ltys = NULL, main = protein, xlab = "average intensity", ylab = "ratio", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{IBSpectra object}
  \item{relative.to}{
      A character factor specifying reporter tag names. Either of 
      length 1 or same length as channels.
}
  \item{protein}{
      Protein group identifier.
}
  \item{noise.model}{
      NoiseModel object.
}
  \item{channels}{
       Reproter tag names.
}
  \item{ylim}{
      See par.
}
  \item{identify}{
      Allows identifiy to be called.
}
  \item{add}{
%%     ~~Describe \code{add} here~~
}
  \item{pchs}{
%%     ~~Describe \code{pchs} here~~
}
  \item{legend.pos}{
%%     ~~Describe \code{legend.pos} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
  \item{legend.cex}{
%%     ~~Describe \code{legend.cex} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
  \item{ltys}{
%%     ~~Describe \code{ltys} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\author{Florian P. Breitwieser}
