\name{IBSpectra-class}
\docType{class}

% Class
\alias{class:IBSpectra}
\alias{IBSpectra-class}
\alias{iTRAQSpectra-class}
\alias{iTRAQSpectra}
\alias{iTRAQ4plexSpectra-class}
\alias{iTRAQ4plexSpectra}
\alias{iTRAQ8plexSpectra-class}
\alias{iTRAQ8plexSpectra}
\alias{TMTSpectra-class}
\alias{TMTSpectra}
\alias{TMT2plexSpectra-class}
\alias{TMT2plexSpectra}
\alias{TMT6plexSpectra-class}
\alias{TMT6plexSpectra}
\alias{IBSpectra}

% Constructor and readIBSpectra function:
\alias{initialize,IBSpectra-method}
\alias{VARMETADATA}

% Coercion methods:
\alias{coerce,IBSpectra,data.frame-method}
\alias{as.data.frame,IBSpectra-method}

% Accessor-like methods:
\alias{proteinGroup}
\alias{proteinGroup,IBSpectra-method}
\alias{proteinGroup<-}
\alias{proteinGroup<-,IBSpectra-method}
\alias{isotopeImpurities}
\alias{isotopeImpurities,IBSpectra-method}
\alias{isotopeImpurities<-}
\alias{isotopeImpurities<-,IBSpectra-method}
\alias{reporterIntensities}
\alias{reporterIntensities,IBSpectra-method}
\alias{reporterIntensities<-}
\alias{reporterIntensities<-,IBSpectra-method}
\alias{reporterMasses}
\alias{reporterMasses,IBSpectra-method}
\alias{reporterMasses<-}
\alias{reporterMasses<-,IBSpectra-method}
\alias{reporterNames}
\alias{reporterNames,IBSpectra-method}
\alias{reporterData}
\alias{reporterData,IBSpectra-method}
\alias{reporterData<-}
\alias{reporterData<-,IBSpectra-method}
\alias{classLabels}
\alias{classLabels,IBSpectra-method}
\alias{classLabels<-}
\alias{classLabels<-,IBSpectra-method}

% Util methods:
\alias{subsetIBSpectra}
\alias{spectrumTitles}
\alias{spectrumTitles,IBSpectra-method}
\alias{spectrumSel}
\alias{spectrumSel,IBSpectra,missing,missing-method}
\alias{spectrumSel,IBSpectra,character,missing-method}
\alias{spectrumSel,IBSpectra,missing,character-method}
\alias{IBSpectraTypes}
\alias{writeData}
\alias{writeData,IBSpectra-method}

% show method:
\alias{show,IBSpectra-method}

\title{IBSpectra objects}

\description{
  This class is a representation for Mass Spectrometry Proteomics spectrum
  level data labeled using Isobaric tags (iTRAQ, TMT). 
  
  It contains per-spectrum meassurements of the reporter tag intensity and 
  m/z in assayData, and protein grouping in proteinGroup.
}

\section{Slots}{
  It extends \code{\link[Biobase:class.eSet]{eSet}}.
  Introduced in \code{eSet}:
   \describe{
      \item{\code{assayData}:}{Contains matrices 'ions' and 'mass',
	i.e. reporter intensities and masses, respectively. Columns are
	reporter tags, rows are spectra. Class:\code{\link{AssayData-class}}}
      \item{\code{phenoData}:}{Contains experimenter-supplied variables
	describing phenotypes behind reporter tags.
	Class: \code{\link{AnnotatedDataFrame-class}}}
      \item{\code{featureData}:}{Contains variables describing spectra
	with their retention time, identified peptide, etc. Class:
	\code{\link{AnnotatedDataFrame-class}}} 
      \item{\code{experimentData}:}{Contains details of experimental
OB	methods. Class: \code{\link{MIAME-class}}}
      \item{\code{annotation}:}{Label associated with the annotation
	package used in the experiment. Class: \code{character}}
      \item{\code{protocolData}:}{Contains equipment-generated
      variables describing reporter tags. Class: \code{\link{AnnotatedDataFrame-class}}}
  }
  
  Introduced in \code{IBSpectra}:
  \describe{
    \item{\code{proteinGroup}:}{A \code{\link{ProteinGroup-class}}
      object describing and grouping the identifications.}
    \item{\code{reporterNames}:}{A character vector denoting the
      reporter tag labels.}
    \item{\code{reporterMasses}:}{The 'true' m/z of the reporter tags
      in the MS/MS spectrum, used to isolate m/z-intensity pairs from
      peaklist.}
    \item{\code{isotopeImpurities}:}{Manufacturer supplied isotope
      impurities, need to be set per batch and used for correction.
    }
  }
}


\section{Constructor}{
  See \code{\link{readIBSpectra}} for creation based on peaklist
  (e.g. MGF format) and identification files (Mascot and Phenyx output).
  \describe{
    \item{\code{new(type,data)}:}{
      Creates a IBSpectra object.
      \describe{
        \item{\code{type}}{Denotes the type of IBSpectra, either
	  'iTRAQ4plexSpectra','iTRAQ8plexSpectra','TMT2plexSpectra' or
	  'TMT6plexSpectra'. Call IBSpectraTypes() to see a list of
	  the implemented types.
	}
        \item{\code{data}}{A 'data.frame' in a ibspectra-csv format.}
      }
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a IBSpectra object.

  \describe{
    \item{\code{as.data.frame(x, row.names = NULL, optional = FALSE)}:}{
      Creates an extended data.frame with all identification and
      quantitation information from x.
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a IBSpectra object.

  \describe{
    \item{\code{subsetIBSpectra(x, protein=NULL, direction="exclude",specificity=c("unspecific","group-specific","reporter-specific"))}:}{ 
      Get a 'subset' of IBSpectra: include or exclude specific proteins. }
    \item{\code{proteinGroup(x)}:}{Gets and sets the proteinGroup.}
    \item{\code{isotopeImpurities(x)}:}{Gets and sets the isotope
    impurities of the isobaric tags as defined by the manufacturers
    per batch.}
  
    \item{\code{reporterData(x,element="ions",...)}:}{
      Gets and sets the element ('ions' or 'mass') for each tag and spectrum.
      '...' is handed down to spectrumSel, so it is possible to
      select for peptides or proteins.
    }
    
    \item{\code{reporterIntensities(x,...)}:}{
      Convenience function, calls reporterData(...,element="ions")
    }
    \item{\code{reporterMasses(x,...)}:}{
      Convenience function, calls reporterData(...,element="mass")
    }

    \item{\code{spectrumSel(x,peptide,protein,specificity="reporter-specific")}:}{
      Gets a boolean vector selecting the corresponding spectra:
      If peptide is given, all spectra assigned to this peptide.
      If protein is given, all spectra assigned to peptides of this
      protein with specificity 'specificity'.
      See also \code{\link[isobar]{ProteinGroup}}.
    }
    \item{\code{spectrumTitles(x,...)}:}{
      Gets the spectrum titles. '...' is passed down to spectrumSel.
    }
    \item{\code{classLabels(x)}:}{
      Gets and sets the class labels in phenoData. Used for
      summarization, see also \code{\link{estimateRatio}} and
      \code{\link[Biobase]{phenoData}}.
    }
  }
}

\author{Florian P. Breitwieser}

\seealso{
  \link{ProteinGroup},
  \link{isobar-preprocessing},
  \link{isobar-analysis},
  \link{isobar-plots}
}

\examples{
data(ibspiked_set1)
ibspiked_set1
}
