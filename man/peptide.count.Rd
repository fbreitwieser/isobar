\name{peptide.count}
\alias{peptide.count}
\alias{spectra.count}
\alias{sequence.coverage}
\title{
  Peptide counts, spectral counts and sequence coverage for ProteinGroup objects.
}
\description{
  Report the peptide count, spectral count and sequence coverage for supplied proteins.
}
\usage{
peptide.count(protein.group, protein.g = reporterProteins(protein.group), 
              specificity = c("reporter-specific", "group-specific", "unspecific"), ...)

spectra.count(protein.group, protein.g = reporterProteins(protein.group), 
              specificity = c("reporter-specific", "group-specific", "unspecific"), 
              modif = NULL, ...)

sequence.coverage(protein.group, protein.g = reporterProteins(protein.group), 
                  specificity = c("reporter-specific", "group-specific", "unspecific"), 
                  simplify = TRUE, ...)
}
\arguments{
  \item{protein.group}{ProteinGroup object.}
  \item{protein.g}{Protein group identifier.}
  \item{specificity}{Specificity of peptides.}
  \item{modif}{Only count peptides having a certain modification.}
  \item{simplify}{If simplify=TRUE, a named numeric vector is returned, 
                  with the mean sequence coverage of the ACs of each protein.g supplied.
                  Else, a list with the length of protein.g is returned having the sequence 
                  coverage for each protein AC.}
  \item{\dots}{Further arguments to \code{\link{peptides}}}
}
\author{
  Florian P Breitwieser
}

\seealso{
  \code{\link{calculate.emPAI}},
  \code{\link{calculate.dNSAF}},
  \code{\link{ProteinGroup}}
}
\examples{
  data(ibspiked_set1)
  sc <- spectra.count(proteinGroup(ibspiked_set1))
  pc <- peptide.count(proteinGroup(ibspiked_set1)) 
  plot(jitter(sc),jitter(pc),log="xy")
}
