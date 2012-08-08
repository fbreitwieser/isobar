\name{spectra.count2}
\alias{spectra.count2}
\title{
  Spectral count for peptides and proteins in ProteinGroup objects.
}
\description{
  Spectral count for peptides and proteins in ProteinGroup objects. 
  It can - other than \code{\link{spectra.count}} - quantify the spectra count on
  the level of peptides, potenitally modifed, too,
}
\usage{
spectra.count2(protein.group, value = reporterProteins(protein.group), 
               type = "protein.g", specificity = c("reporter-specific", "group-specific", "unspecific"), 
               modif = NULL, combine = FALSE, ...)
}
\arguments{
  \item{protein.group}{ProteinGroup object.}
  \item{value}{List of protein group identifiers or peptides.}
  \item{type}{Either 'protein.g' or 'peptide'.}
  \item{specificity}{Specificity of peptides.}
  \item{modif}{Only count peptides having a certain modification.}
  \item{combine}{If TRUE, only one combined result is returned.}
  \item{\dots}{Further arguments to \code{\link{peptides}}}
}
\author{Florian P Breitwieser}

\seealso{
  \code{\link{spectra.count}},
  \code{\link{ProteinGroup}}
}
\examples{
  data(ibspiked_set1)
  pg <- proteinGroup(ibspiked_set1)
  protein.gs <- protein.g(pg,"CERU")
  sc <- spectra.count2(pg,protein.gs)
  sc.ik <- spectra.count2(pg,protein.gs,modif="iTRAQ4plex_K")
  rbind(spectra.counts=sc,spectra.counts_iTRAQk=sc.ik)
}
