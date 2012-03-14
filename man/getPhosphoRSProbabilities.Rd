\name{getPhosphoRSProbabilities}
\alias{getPhosphoRSProbabilities}
\alias{readPhosphoRSOutput}
\alias{writePhosphoRSInput}

\title{
  Get phosphoRS probabilities
}
\description{
  Get phosphorylation site localization probabilities by calling PhosphoRS
  and parsing its output. \code{getPhosphoRSProbabilities} generates a XML
  input file for PhosphoRS calling \code{writePhosphoRSInput}, then executes
  phosphoRS.jar with java, and parses the XML result file with 
  \code{readPhosphoRSOutput}.
}
\usage{
  getPhosphoRSProbabilities(id.file, mgf.file, massTolerance, activationType, simplify = FALSE,
                            mapping.file = NULL, mapping = c(peaklist = "even", id = "odd"), 
                            phosphors.jar = system.file("phosphors", "phosphoRS.jar", package = "isobar", mustWork = TRUE))

  writePhosphoRSInput(phosphoRS.infile, id.file, mgf.file, massTolerance, activationType, 
                      mapping.file = NULL, mapping = c(peaklist = "even", id = "odd"))

  readPhosphoRSOutput(phosphoRS.outfile, simplify = FALSE)
}
\arguments{
  \item{id.file}{Database search results file in \code{ibspectra.csv} or
    \code{mzIdentML} format. See \code{\link{IBSpectra}} and isobar vignette for
    information on converting Mascot dat and Phenyx pidres files
    into ibspectra format.}
  \item{mgf.file}{Peaklist file}
  \item{massTolerance}{Fragment ion mass tolerance (in Da)}
  \item{activationType}{Activation types of spectra. CID, HCD, or ETD.}
  \item{simplify}{If \code{TRUE}, returns a \code{data.frame} instead of a \code{list}.}
  \item{mapping.file}{Mapping file. See also \code{\link{readIBSpectra}}.}
  \item{mapping}{Mapping columns. }
  \item{phosphors.jar}{Location of the PhosphoRS jar file.}
  \item{phosphoRS.infile}{PhosphoRS input XML file name.}
  \item{phosphoRS.outfile}{PhosphoRS output XML file name.}
}

\details{
  PhosphoRS is described in Taus et al., 2011.
  It can be downloaded from http://cores.imp.ac.at/protein-chemistry/download/ 
  and used as Freeware. Java is required at runtime.
}
\value{
  If \code{simplify=TRUE}, a \code{data.frame} with the following columns:
    \code{spectrum, peptide, modif, PepScore, PepProb, seqpos}

  If \code{simplify=FALSE}, a list (of spectra) of lists (of peptide 
    identifications) of lists (with information about identification 
    and localization).
    spectrum -> {peptide 1, peptides 2, ...} -> {peptide}.
    First level:
      - spectrum
    Second level:
      - peptide identifications for spectrum (might be more than one)
    Third level:
      - peptide: vector with peptide sequence and modification stirng
      - site.probs: matrix with site probabilities for each phospho site
      - isoforms: peptide score and probabilities for each isoform  
}
\references{
  Taus et al., 2011
}
\author{
  Florian P Breitwieser
}

\keyword{ ~phospho }
