\name{isobar-plots}

\alias{isobar-plots}
\alias{reporterMassPrecision}
\alias{reporterMassPrecision,IBSpectra,missing-method}
\alias{reporterMassPrecision,IBSpectra,logical-method}
\alias{reporterIntensityPlot}
\alias{reporterIntensityPlot-methods}
\alias{reporterIntensityPlot,IBSpectra-method}
\alias{raplot}
\alias{raplot,IBSpectra-method}
\alias{plotRatio}
\alias{plotRatio,IBSpectra,character,character,character-method}
\alias{maplot}
\alias{maplot,IBSpectra,missing,missing-method}
\alias{maplot,IBSpectra,character,character-method}
\alias{maplot,missing,numeric,numeric-method}
\alias{maplot2}
\alias{maplot2,ANY,character,character-method}
\alias{maplot2,list,character,character-method}
\alias{protGgdata}
\alias{protGgdata,ANY,character,character-method}

\title{IBSpectra plots}

\description{
  Various plots are implement to assure data quality, and accompany
  preprocessing  and analysis.
}

\section{reporterMassPrecision}{
  \describe{
    \item{\code{reporterMassPrecision(x)}:}{
      Calculates and displays the deviation from the 'true' tag mass - 
      as specified in the IBSpectra object - of each channel.
  }}
}

\section{reporterIntensityPlot}{
  \describe{
    \item{\code{reporterIntensityPlot(x)}:}{
      Displays boxplots of intensity of channels before and after
      normalization - useful to check the result of normalization.
  }}
}


\section{raplot}{
  \describe{
    \item{\code{raplot(x,...)}:}{
      Ratio-Absolute intensity plot - will be deprecated by maplot
      \describe{
        \item{\code{x}}{IBSpectra object}
      	\item{\code{...}}{Parameters to plot function.}
    }
  }}
}

\section{plotRatio}{
  \describe{
    \item{\code{plotRatio(x,channel1,channel2,protein,...)}:}{
      Plots abundances of one protein
      
      \describe{
        \item{\code{x}}{IBSpectra object}
        \item{\code{channel1}}{}
        \item{\code{channel2}}{}
        \item{\code{protein}}{}
      	\item{\code{...}}{Parameters to plot function.}
      }
    }
  }
}


\section{maplot}{
  \describe{
    \item{\code{maplot(x,channel1,channel2,...)}:}{
      Creates a ratio-versus-intensity plot.
      \describe{
        \item{\code{x}}{IBSpectra object.}
    }
  }}
}

\section{maplot2}{
  \describe{
    \item{\code{maplot2()}:}{
  }}
}


\author{Florian P. Breitwieser, Jacques Colinge}

\seealso{
  \link{IBSpectra},
  \link{isobar-preprocessing}
  \link{isobar-analysis}
}

\examples{
data(ibspiked_set1)
maplot(ibspiked_set1,main="IBSpiked, not normalized")
maplot(normalize(ibspiked_set1),main="IBSpiked, normalized")

}
