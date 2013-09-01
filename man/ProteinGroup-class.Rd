\name{ProteinGroup-class}
\docType{class}

% Class
\alias{class:ProteinGroup}
\alias{ProteinGroup-class}
\alias{ProteinGroup}

% Constructor-like function:
\alias{ProteinGroup}
\alias{ProteinGroup,data.frame,missing-method}
\alias{ProteinGroup,data.frame,NULL-method}
\alias{ProteinGroup,data.frame,ProteinGroup-method}
\alias{readProteinGroup}
\alias{readProteinGroup2}

% Coercion methods:
\alias{coerce,data.frame,ProteinGroup-method}
\alias{as.data.frame,ProteinGroup-method}
\alias{as.data.frame.ProteinGroup}
\alias{proteinGroup.as.concise.data.frame}

% Accessor-like methods:
\alias{peptideSpecificity}
\alias{peptideSpecificity,ProteinGroup-method}
\alias{peptideNProtein}
\alias{peptideNProtein,ProteinGroup-method}
\alias{proteinGroupTable}
\alias{proteinGroupTable,ProteinGroup-method}
\alias{peptideSpecificity}
\alias{peptideSpecificity,ProteinGroup-method}
\alias{peptideInfo}
\alias{peptideInfo-methods}
\alias{peptideInfo,ProteinGroup-method}

\alias{indistinguishableProteins}
\alias{indistinguishableProteins,ProteinGroup-method}
\alias{indistinguishableProteins,ProteinGroup,ANY,ANY-method}
\alias{indistinguishableProteins,ProteinGroup,character,missing-method}
\alias{indistinguishableProteins,ProteinGroup,missing,character-method}
\alias{indistinguishableProteins,ProteinGroup,missing,missing-method}
\alias{spectrumToPeptide}
\alias{spectrumToPeptide,ProteinGroup-method}
\alias{get.pep.group}
\alias{peptides}
\alias{peptides,ProteinGroup,missing-method}
\alias{peptides,ProteinGroup,character-method}
\alias{reporterProteins}
\alias{reporterProteins,ProteinGroup-method}
\alias{protein.ac}
\alias{protein.g}
\alias{protein.ac,ProteinGroup,character-method}
\alias{protein.ac,ProteinGroup,missing-method}
\alias{protein.g,ProteinGroup,character-method}
\alias{protein.g,ProteinGroup,character,character-method}

% show and method:
\alias{show,ProteinGroup-method}
\alias{summary.ProteinGroup}


\title{ProteinGroup objects}

\description{
  The ProteinGroup class is a container for identified peptides and 
	  proteins, and groups them to distinguish proteins with specific 
	  peptides.
}

\usage{
ProteinGroup(from,template=NULL,proteinInfo=data.frame())

protein.ac(x, protein.g)
protein.g(x, pattern, variables=c("AC","name"), ...)
}

\arguments{
  \item{from}{\code{data.frame} object to create a ProteinGroup from.
  See Details from column specifications}
  \item{template}{'template' ProteinGroup object for grouping.}
  \item{x}{ProteinGroup object}
  \item{protein}{character string}
  \item{proteinInfo}{data.frame for proteinInfo slot}
  \item{protein.g}{character string, denoting a 'protein group'.}
  \item{pattern}{character string, see \code{\link{grep}} for details.}
  \item{variables}{\code{AC} maps a protein accession code to a protein group.
  \code{name} maps using protein information from \code{proteinInfo}.}
  \item{\dots}{Passed on to \code{\link{grep}}.}
}

\details{
  The ProteinGroup class stores spectrum to peptide to protein mapping.
  
  The proteins are grouped by their evidence, i. e. peptides:
  \itemize{
    \item Peptides with changes only from Leucin to Isoleucin are
      considered the same, as they cannot be distinguished by MS.
    \item Proteins which are detected with the same peptides are
      grouped together to a 'indistinguishable protein'- normally
      these are splice variants.
    \item Proteins with specific peptides are 'reporters'.
    \item Proteins with no specific peptides are grouped under these 'reporters.
  }

  This information is stored in six slots:
  \describe{
    \item{\code{spectra.n.peptides}}{a named 'character' vector, names
  being spectrum identifier and values are peptides.}
    \item{\code{peptide.n.proteins}}{a 'data.frame' containing the
  number of proteins the peptides could derive from.}
    \item{\code{peptide.n.protein}}{a character 'matrix' linking
  peptides to proteins.}
    \item{\code{indistinguishable.proteins}}{a 'matrix' contain.}
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{ProteinGroup(tbl.prot.pep,template=NULL)}:
      Creates a ProteinGroup object.
      \describe{
        \item{\code{tbl.prot.pep}}{A 'data.frame' with three columns:
	  1. Protein, 2. Peptide, 3. Spectrum.}
        \item{\code{template}}{Optional ProteinGroup object the grouping
  is based upon.}
      }
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a ProteinGroup object.

  \describe{
    \item{\code{as(from, "ProteinGroup")}:}{
      Creates a ProteinGroup object from a data.frame.
    }
    \item{\code{as.data.frame(x, row.names = NULL, optional = FALSE)}:}{
      Creates a data.frame with columns \code{protein} (character),
      \code{peptide} (character), \code{spectrum}.}
    \item{\code{as.concise.data.frame(from)}:}{
      Creates a 'concise' data.frame with one spectrum per row, and protein ACs combined
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a ProteinGroup object.

  \describe{
    \item{\code{spectrumToPeptide(x)}:}{ Gets spectrum to peptide assignment.}
    \item{\code{peptideInfo(x)}:}{ Peptide information such as protein start position.}
    \item{\code{peptideSpecificity(x)}:}{ Gets a 'data.frame'
    containing the peptide specificity: they can be reporter-specific,
    group-specific, or non-specific.}
    \item{\code{peptideNProtein(x)}:}{ Gets peptide to protein assignment.}
    \item{\code{indistinguishableProteins(x)}:}{
      Gets the proteins which cannot be distinguished based on peptide evidence.
    }
    \item{\code{proteinGroupTable}:}{ Gets the protein grouping, listing
    reporters and group members. }
    \item{\code{peptides(x,protein=NULL,specificity=c("reporter-specific",
	"group-specific","unspecific"),columns="peptide",set=union)}:}{
      Gets all peptides detected, or just those for a protein with the
      defined specificity. \code{columns} might define multiple columns of
      \code{peptideSpecificity(x)}. set=union returns the union of
      peptides of all proteins defined, set=intersect returns the intersection.
      
    }
  }
}

\author{Florian P. Breitwieser}

\seealso{
  \link[isobar]{IBSpectra}
}

\examples{

tbl <- data.frame(spectrum=1:14,peptide=c(rep(letters[1:3],4),"a","x"),
                  modif=":",start.pos=1,
                  protein=c(rep(c("A","B"),each=6),"C","D"))
pg <- ProteinGroup(tbl)
pg
proteinGroupTable(pg)

data(ibspiked_set1)
pg <- proteinGroup(ibspiked_set1)
ceru.proteins <- protein.g(pg,"CERU")

## all ceru peptides
peptides(pg,ceru.proteins)

## peptides shared by all ceru proteins
peptides(pg,ceru.proteins, set=intersect) 

}
