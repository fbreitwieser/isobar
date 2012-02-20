\name{isobar-preprocessing}

\alias{isobar-preprocessing}
\alias{correctIsotopeImpurities}
\alias{correctIsotopeImpurities,IBSpectra-method}
\alias{normalize}
%\alias{normalize,IBSpectra-method}
\alias{subtractAdditiveNoise}
\alias{subtractAdditiveNoise,IBSpectra-method}
\alias{exclude}
\alias{exclude,IBSpectra,character-method}

\title{IBSpectra preprocessing}

\description{
  Preprocessing is a necessary step prior to analysis of data.
  In a sequential order, it is often neccassary to correct isotope
  impurities, to normalize, and subtract additive noise.
}

\section{Isotope impurity correction}{
  \describe{
    \item{\code{correctIsotopeImpurities(x)}:}{
      Returns impurity corrected IBSpectra object by solving a linear
      system of equations.
      See also \link{isotopeImpurities}.
      \describe{
        \item{\code{}}{}
    }
  }}
}

\section{Normalization}{
  \describe{
    \item{\code{normalize(x,f=median,target="intensity",exclude.protein=NULL,
	use.protein=NULL,f.doapply=TRUE,log=TRUE,channels=NULL,na.rm=FALSE)}:}{
      Normalizes the intensities for multiplicative errors. Those changes
      are most likely produced by pipetting errors, and different
      hybridization efficencies, but can also be due to biological
      reasons. 
      By default, tag intensities are multiplied by a factor so that the
      median intensity is equal across tags.
      \describe{
        \item{\code{f}:}{\code{f} is applied to each column, unless f.doapply is FALSE. 
          Then \code{f} is supposed to compute column-wise statistics of the matrix of 
          intensities. E.g. \code{colSums} and \code{colMeans}.}
	    \item{\code{target}:}{One of "intensity" and "ratio".}
      	\item{\code{exclude.proteins}}{
          Spectra of peptides which might come from these proteins are excluded. Use for example 
          for contaminants and proteins depleted in the experiment.}
        \item{\code{use.protein}:}{
          If specified, only spectra coming from this protein are used. Use when a protein 
          is spiked-in as normalization control.}
        \item{\code{f.isglobal}:}{If true, \code{f} is applied on each column. 
            If false, \code{f} is supposed to compute column-wise statistics of the matrix 
            of intensities. E.g. \code{colSums} and \code{colMeans}.} 
        \item{\code{log}:}{Used when target=ratio.}
      }
    }
  }
}

\section{Substract additive noise}{
  \describe{
    \item{\code{subtractAdditiveNoise(x,method="quantile",shared=TRUE,prob=0.01)}:}{
      
      \describe{
        \item{\code{method}}{'quantile' method is supported for now. It
	  take's the prob (0.01) quantile to estimate the noise level.
	  This value is subtracted from all intensities, and all
	  remaining intensities have to be at least that value.
	}
	\item{\code{prob}}{See 'method'.}
	\item{\code{shared}}{
	  If channels are assumed similar in intensity and hence a
	  shared noise level is reasonable.
	  If not, then one level per channel is necessary.
	}

      }
    }
  }
}


\section{Exclusion of proteins}{
  \describe{
    \item{\code{exclude(x,proteins.to.exclude)}:}{
      Removes spectra which are assigned to proteins in
      protein.to.exclude from the object. This can be useful to remove
      contaminants. It create a new grouping based on the data which
      is left.
      \describe{
        \item{\code{proteins.to.exclude}}{Proteins to exclude.}
    }
  }}
}


\author{Florian P. Breitwieser, Jacques Colinge}

\seealso{
  \link{ProteinGroup},
  \link{IBSpectra},
  \link{isobar-analysis},
  \link{isobar-plots}
}

\examples{
data(ibspiked_set1)
maplot(ibspiked_set1,main="IBSpiked, not normalized")
maplot(normalize(ibspiked_set1),main="IBSpiked, normalized")

}
