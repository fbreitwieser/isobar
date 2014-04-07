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
\alias{TMT6plexSpectra2-class}
\alias{TMT6plexSpectra}
\alias{TMT6plexSpectra2}
\alias{TMT10plexSpectra-class}
\alias{TMT10plexSpectra}
\alias{IBSpectra}

% Constructor and readIBSpectra function:
\alias{initialize,IBSpectra-method}
\alias{VARMETADATA}

% Coercion methods:
\alias{coerce,IBSpectra,data.frame-method}
\alias{as.data.frame,IBSpectra-method}
\alias{as.data.frame.IBSpectra}

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
\alias{reporterTagNames}
\alias{reporterTagNames,IBSpectra-method}
\alias{reporterTagMasses}
\alias{reporterTagMasses,IBSpectra-method}
\alias{reporterData}
\alias{reporterData,IBSpectra-method}
\alias{reporterData<-}
\alias{reporterData<-,IBSpectra-method}
\alias{classLabels}
\alias{classLabels,IBSpectra-method}
\alias{classLabels<-}
\alias{classLabels<-,IBSpectra-method}

% Util methods:

\alias{spectrumTitles}
\alias{spectrumTitles,IBSpectra-method}
\alias{spectrumSel}
\alias{spectrumSel,IBSpectra,missing,missing-method}
\alias{spectrumSel,IBSpectra,matrix,missing-method}
\alias{spectrumSel,IBSpectra,character,missing-method}
\alias{spectrumSel,IBSpectra,missing,character-method}
\alias{spectrumSel,IBSpectra,data.frame,missing-method}
\alias{IBSpectraTypes}
\alias{writeData}
\alias{writeData,IBSpectra-method}
\alias{ibSpectra.as.concise.data.frame}

% show method:
\alias{show,IBSpectra-method}

\title{IBSpectra Class for Isobarically Tagged Quantitative MS
  Proteomics Data}

\description{ This class represents a quantitative MS proteomics
  experiment labeled using Isobaric tags (iTRAQ, TMT). \code{IBSpectra}
  is a abstract class which is implemented in the
  \code{\link{IBSpectraTypes}} classes \code{iTRAQ4plexSpectra},
  \code{iTRAQ8plexSpectra}, \code{TMT2plexSpectra}, \code{TMT6plexSpectra} and \code{TMT10plexSpectra}.
    
  It contains per-spectrum meassurements of the reporter tag intensity and 
  m/z in assayData, and protein grouping in proteinGroup.
}

\section{Objects from the Class}{
  \code{IBSpectra} objects are typically created using the
  \code{\link{readIBSpectra}} method or by calls of the form
  \code{new("iTRAQ4plexSpectra",data=NULL,data.ions=NULL,...).}
}

\section{Slots}{
  \code{IBSpectra} extends \code{\link[Biobase:class.eSet]{eSet}} which
  is a container for high-throughput assays and experimental metadata.
  Slots introduced in \code{eSet} (for more details on slots and methods
  refer to \code{\link{eSet}} help):
   \describe{
     \item{\code{assayData}:}{Contains matrices 'ions' and 'mass
       storing reporter tag intensities and m/z values for each tag and
       spectrum. Can be accessed by \code{\link{reporterIntensities}}
       and \code{\link{reporterMasses}}.  Class: \code{\link{AssayData}}}
      \item{\code{phenoData}:}{Contains experimenter-supplied variables
	describing phenotypes behind reporter tags.
	Class: \code{\link{AnnotatedDataFrame-class}}}
      \item{\code{featureData}:}{Describes the spectra's retention time,
	charge, peptide sequence, etc and can be accessed by
	\code{\link{fData}}. Class: \code{\link{AnnotatedDataFrame}}} 
      \item{\code{experimentData}:}{Contains details of experimental
	methods. Class: \code{\link{MIAME}}}
      \item{\code{annotation}:}{UNUSED. Label associated with the annotation
	package used in the experiment. Class: \code{character}}
      \item{\code{protocolData}:}{UNUSED. Contains equipment-generated
      variables describing reporter tags. Class: \code{\link{AnnotatedDataFrame}}}
      \item{\code{log}:}{\code{character matrix} logging isotope impurity
        correction, normalization, etc.}
   }
  
  Slots introduced in \code{IBSpectra}:
  \describe{
    \item{\code{proteinGroup}:}{A \code{\link{ProteinGroup}}
      object describing peptide and protein identifications
      grouped by shared peptides.}
    \item{\code{reporterTagNames}:}{A character vector denoting the
      reporter tag labels.}
    \item{\code{reporterMasses}:}{The 'true' m/z of the reporter tags
      in the MS/MS spectrum, used to isolate m/z-intensity pairs from
      peaklist.}
    \item{\code{isotopeImpurities}:}{Manufacturer supplied isotope
      impurities, need to be set per batch and used for correction
      by \code{\link{correctIsotopeImpurities}}.
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
	  'iTRAQ4plexSpectra','iTRAQ8plexSpectra','TMT2plexSpectra', 'TMT6plexSpectra'
	  or 'TMT10plexSpectra'. Call IBSpectraTypes() to see a list of
	  the implemented types.
	}
        \item{\code{data}}{A 'data.frame' in a ibspectra-csv format.}
      }
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a IBSpectra object.
  \code{IBSpectra} object can be coerced to 
  \describe{
    \item{\code{as(x, "data.frame")}:}{ Creates a data.frame containing
      all identification and quantitation information. Peptide matching
      to multiple proteins produce multiple lines. }
    \item{\code{ibSpectra.as.concise.data.frame(x)}:}{ Creates a data.frame
      containing all identification and
      quantitation information. Proteins are concatenated - so the
      resulting \code{data.frame} has one line per spectrum.}
    \item{\code{as(x, "MSnSet")}:}{
      Coerces to a \code{\link[MSnbase]{MSnSet}} object (package
      \code{\link[MSnbase]{MSnbase}}). }
    \item{\code{as(msnset,"IBSpectra")}:}{
      Coerces a \code{\link[MSnbase]{MSnSet}} to \code{IBSpectra} object.
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a IBSpectra object.

  \describe{
    \item{\code{proteinGroup(x)}:}{Gets and sets the ProteinGroup.}
    \item{\code{isotopeImpurities(x)}:}{Gets and sets the isotope
    impurities of the isobaric tags as defined by the manufacturers
    per batch.}  
    \item{\code{reporterData(x,element="ions",na.rm=FALSE,na.rm.f='any',...)}:}{
      Gets and sets the element ('ions' or 'mass') for each tag and spectrum.
      '...' is handed down to spectrumSel, so it is possible to
      select for peptides or proteins.
      If na.rm is TRUE, than spectra missing quantitative information in 'any' or 'all'
      channels (parameter na.rm.f) are removed.
    }
    
    \item{\code{reporterIntensities(x,...)}:}{
      Convenience function, calls reporterData(...,element="ions")
    }
    \item{\code{reporterMasses(x,...)}:}{
      Convenience function, calls reporterData(...,element="mass")
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

\section{Methods}{
  In the following code snippets, \code{x} is a IBSpectra object.

  \describe{
    \item{\code{subsetIBSpectra(x, protein=NULL, peptide=NULL, direction="exclude",specificity)}:}{
      Get a 'subset' of IBSpectra: include or exclude proteins
      or peptides. When selection  is based on proteins, it can be
      defined to exclude only peptides which are specific to the protein
      ('reporter-specific'), specific to the group ('group-specific') or
      which are shared with other proteins ('unspecific'). See \code{\link{subsetIBSpectra}}.}
    \item{\code{spectrumSel(x,peptide,protein,specificity="reporter-specific")}:}{
      Gets a boolean vector selecting the corresponding spectra:
      If peptide is given, all spectra assigned to this peptide.
      If protein is given, all spectra assigned to peptides of this
      protein with specificity 'specificity'.
      See also \code{\link[isobar]{ProteinGroup}}.
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
head(reporterIntensities(ibspiked_set1))
head(reporterMasses(ibspiked_set1))
proteinGroup(ibspiked_set1)
isotopeImpurities(ibspiked_set1)

# create new object
set.seed(123)
data <- data.frame(spectrum=letters,
                   peptide=sample(c("pepA","pepB","pepC"),26,TRUE),
                   start.pos=1,
                   modif=sample(c("::X:::",":Y::::","::Z:::"),26,TRUE),
                   accession=c("protein1","protein2"))
data.ions <- matrix(rnorm(26*2,1000,50),
                    ncol=2,dimnames=list(letters,NULL))
data.mass <- matrix(rep(c(126.1,127.1),26),
                    ncol=2,byrow=TRUE,dimnames=list(letters,NULL))
ib <- new("TMT2plexSpectra",data,data.ions,data.mass)
ib
reporterIntensities(ib)
isotopeImpurities(ib) <- matrix(c(0.8,0.1,0.2,0.9),nrow=2)
reporterIntensities(correctIsotopeImpurities(ib))
}
