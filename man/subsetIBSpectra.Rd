\name{subsetIBSpectra}
\alias{subsetIBSpectra}
\title{
  Subset IBSpectra objects
}
\description{
  Returns an IBSpectra object which is a subset of the input, excluding
  or exclusively containing the peptides or proteins supplied.
}
\usage{
subsetIBSpectra(x, protein = NULL, peptide = NULL,
                direction = "exclude",
                specificity = c(REPORTERSPECIFIC, GROUPSPECIFIC, UNSPECIFIC), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{IBSpectra object.}
  \item{protein}{Protein group identifiers. Use \code{\link{protein.g}}
    to get protein group identifiers from protein database ACs.}
  \item{peptide}{Peptide sequences.}
  \item{direction}{either 'include' or 'exclude'.}
  \item{specificity}{When 'protein' is supplied: Which peptides should
    be selected? See \code{\link{specificities}}.}
  \item{\dots}{Further arguments passed to \code{\link{spectrumSel}}}
}
\author{Florian P Breitwieser}

\seealso{
  \code{\link{protein.g}}, \code{\link{spectrumSel}}, \code{\link{specificities}}
}
\examples{
data(ibspiked_set1)

# get Keratin proteins
keratin.proteins <- protein.g(proteinGroup(ibspiked_set1),"Keratin")

# exclude Keratin proteins
subsetIBSpectra(ibspiked_set1,protein=keratin.proteins,direction="exclude")
}