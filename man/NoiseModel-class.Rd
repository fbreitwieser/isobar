\name{NoiseModel-class}
\docType{class}

% Class:
\alias{class:NoiseModel}
\alias{NoiseModel-class}
\alias{ExponentialNoiseModel-class}
\alias{ExponentialNoANoiseModel-class}
\alias{InverseNoiseModel-class}
\alias{InverseNoANoiseModel-class}
\alias{GeneralNoiseModel-class}

% Constructor:
\alias{initialize,NoiseModel-method}
\alias{NoiseModel}
\alias{NoiseModel,IBSpectra-method}

% Accessor methods:
\alias{variance}
\alias{variance,NoiseModel,numeric,numeric-method}
\alias{variance,NoiseModel,numeric,missing-method}

\alias{stddev}
\alias{stddev,NoiseModel-method}

\alias{noiseFunction}
\alias{noiseFunction,NoiseModel-method}

\alias{parameter}
\alias{parameter<-}
\alias{parameter,NoiseModel-method}
\alias{parameter<-,NoiseModel-method}

\alias{lowIntensity}
\alias{lowIntensity<-}
\alias{lowIntensity,NoiseModel-method}
\alias{lowIntensity<-,NoiseModel-method}

\alias{naRegion}
\alias{naRegion<-}
\alias{naRegion,NoiseModel-method}
\alias{naRegion<-,NoiseModel-method}

% show method:
\alias{show,NoiseModel-method}

% plot method:
\alias{plot.NoiseModel}

\title{NoiseModel objects}

\description{
  A NoiseModel represent the technical variation which is dependent on
  signal intensity.
}


\section{Constructor}{
  \describe{
    \item{\code{new(type,ibspectra,reporterTagNames=NULL,one.to.one=TRUE,min.spectra=10,plot=FALSE,
                    pool=FALSE)}:}{
    Creates a new NoiseModel object based on ibspectra object.
      \describe{
        \item{\code{type}:}{A non-virtual class deriving from NoiseModel:
          \code{ExponentialNoiseModel}, \code{ExponentialNoANoiseModel},
          \code{InverseNoiseModel}, \code{InverseNoANoiseModel}}
        \item{\code{reporterTagNames}:}{When NULL, all channels from ibspectra are taken
          (i.e. \code{sampleNames(ibspectra)}). Otherwise, specify
	  subset of names, or a matrix which defines the desireed combination of channels (nrow=2).}
	\item{\code{one.to.one}:}{Set to false to learn noise model one a non
	  one-to-one dataset}
	\item{\code{min.spectra}:}{When one.to.one=FALSE, only take proteins
	  with min.spectra to learn noise model.}
	\item{\code{plot}:}{Set to true to plot data the noise model is learnt on.}
    \item{\code{pool}:}{If false, a NoiseModel is estimated on each combination
         of channels indivdually, and then the parameters are averaged. If true,
         the ratios of all channels are pooled and then a NoiseModel is estimated. }
      }
    }
  }
}

\section{Accessor methods}{
  \describe{
    \item{\code{noiseFunction}:}{Gets the noise function.}
    \item{\code{parameter}:}{Gets and sets the parameters for the noise function.}
    \item{\code{variance}:}{Gets the variance for data points based on
      the noise function and parameters.}
    \item{\code{stddev}:}{Convenience function, \code{sqrt(variance(...))}.}
    \item{\code{lowIntensity}:}{Gets and sets the low intensity slot, denoting the noise region.}
    \item{\code{naRegion}:}{Gets and sets the na.region slot.}
  }
}

\examples{

data(ibspiked_set1)

ceru.proteins <- protein.g(proteinGroup(ibspiked_set1),"CERU")

# normalize
ibspiked_set1 <- normalize(correctIsotopeImpurities(ibspiked_set1))

# remove spiked proteins
ibspiked_set1.noceru <- exclude(ibspiked_set1,ceru.proteins)
ibspiked_set1.justceru <- subsetIBSpectra(ibspiked_set1,protein=ceru.proteins,direction="include")

# learn noise models
nm.i <- new("InverseNoiseModel",ibspiked_set1.noceru)
nm.e <- new("ExponentialNoiseModel",ibspiked_set1.noceru)

#learn on non-one.to.one data: not normalized, with spiked proteins
nm.n <- new("ExponentialNoiseModel",ibspiked_set1.justceru,one.to.one=FALSE)

maplot(ibspiked_set1,noise.model=c(nm.e,nm.i,nm.n),ylim=c(0.1,10))

}
