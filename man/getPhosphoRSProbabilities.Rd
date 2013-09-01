\name{getPhosphoRSProbabilities}
\alias{getPhosphoRSProbabilities}
\alias{readPhosphoRSOutput}
\alias{writePhosphoRSInput}
\alias{filterSpectraPhosphoRS}

\title{
  Generate input files for PhosphoRS, call it, and get modification site probabilities
}
\description{
  Get phosphorylation site localization probabilities by calling PhosphoRS
  and parsing its output. \code{getPhosphoRSProbabilities} generates a XML
  input file for PhosphoRS calling \code{writePhosphoRSInput}, then executes
  phosphoRS.jar with java, and parses the XML result file with 
  \code{readPhosphoRSOutput}.
}
\usage{
  getPhosphoRSProbabilities(id.file, mgf.file, massTolerance, activationType,
                 simplify = FALSE, mapping.file = NULL, mapping =
                 c(peaklist = "even", id = "odd"), pepmodif.sep =
                 "##.##", besthit.only = TRUE, phosphors.cmd =
                 paste("java -jar", system.file("phosphors",
                 "phosphoRS.jar", package = "isobar")), file.basename =
                 tempfile("phosphors."))


  writePhosphoRSInput(phosphoRS.infile, id.file, mgf.file, massTolerance,
                 activationType, mapping.file = NULL, mapping =
                 c(peaklist = "even", id = "odd"), pepmodif.sep =
                 "##.##", modif.masses = rbind(c("PHOS", "1",
                 "1:Phospho:Phospho:79.966331:PhosphoLoss:97.976896:STY"),
                 c("Oxidation_M", "2",
                 "2:Oxidation:Oxidation:15.994919:null:0:M"),
                 c("Cys_CAM", "3",
                 "3:Carbamidomethylation:Carbamidomethylation:57.021464:null:0:C"),
                 c("iTRAQ4plex", "4",
                 "4:iTRAQ4:iTRAQ4:144.1544:null:0:KX"), c("iTRAQ8plex",
                 "5", "5:iTRAQ8:iTRAQ8:304.308:null:0:KX"),
                 c("TMT6plex", "7",
                 "7:TMT6:TMT6:229.162932:null:0:KX"), c("TMTsixplex",
                 "6", "6:TMT6:TMT6:229.162932:null:0:KX")))


  readPhosphoRSOutput(phosphoRS.outfile, simplify = FALSE, pepmodif.sep = "##.##", besthit.only = TRUE)

  filterSpectraPhosphoRS(id.file, mgf.file, ..., min.prob = NULL, do.remove=FALSE)
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
  \item{besthit.only}{Only show best hit, simplifies result to data.frame instead of list.}
  \item{phosphors.cmd}{PhosphoRS script.}
  \item{file.basename}{Base name for creating phosphoRS input and output files.}
  \item{phosphoRS.infile}{PhosphoRS input XML file name.}
  \item{phosphoRS.outfile}{PhosphoRS output XML file name.}
  \item{pepmodif.sep}{separator of peptide and modification in XML id}
  \item{modif.masses}{masses and ID used for PhosphoRS}
  \item{min.prob}{Threshold for PhosphoRS peptide probability to consider it for quantification}
  \item{...}{Further arguments to getPhosphoRSProbabilities}
  \item{do.remove}{If TRUE, spectra below the min.prob threshold are not just set as 'use.for.quant=FALSE' but removed.}
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
