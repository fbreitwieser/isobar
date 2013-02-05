\name{calculate.emPAI}
\alias{calculate.emPAI}
\alias{n.observable.peptides}
\alias{observable.peptides}
\title{
  emPAI approximate abundance calculations.
}
\description{
  The Exponentially Modified Protein Abundance Index (emPAI) is a
  label free quantitative measure of protein abundance based on protein
  coverage by peptide matches. The original publication is 
  Ishihama Y, et al., Proteomics (2005).
}
\usage{
calculate.emPAI(protein.group, protein.g = reporterProteins(protein.group), normalize = FALSE,
                observed.pep = c("pep", "mod.charge.pep"), use.mw = FALSE, combine.f = mean, 
                ..., nmc = 0, report.all = FALSE)
n.observable.peptides(...)
observable.peptides(seq, nmc = 1, min.length = 6, min.mass = 600, max.mass = 4000, 
                    custom = list(code = c("B", "Z", "J", "U"),
                                  mass = c(164.554862, 278.61037, 213.12392, 150.953636)), ...)
}
\arguments{
  \item{protein.group}{ProteinGroup object. Its \code{@proteinInfo}
    slot \code{data.frame} must contain a \code{sequence} column to
    calculate the number of observable peptides per protein.
  }
  \item{protein.g}{Protein group identifiers.}
  \item{normalize}{Normalize to sum = 1?.}
  \item{observed.pep}{What counts as observed peptide?}
  \item{report.all}{TOADD}
  \item{use.mw}{Use MW to normalize for protein size}
  \item{combine.f}{How to handle proteins seen only with shared peptides?}
  \item{seq}{Protein sequence.}
  \item{nmc}{Number of missed cleavages.}
  \item{min.length}{Minimum length of peptide.}
  \item{min.mass}{Minimum mass of peptide.}
  \item{max.mass}{Maximum mass of peptide.}
  \item{custom}{User defined residue for \code{\link[OrgMassSpecR]{Digest}}.}
  \item{\dots}{Further arguments to \code{\link{observable.peptides}}/\code{\link[OrgMassSpecR]{Digest}}.}
}
\details{
  The formula is \deqn{emPAI = 10^{\frac{N <- {observed}}{N <- {observable}}} -1}{emPAI = 10^(N_observed/N_observable) -1}
  N_observed is the number of observed peptides - we use the count
  of unique peptide without consideration of charge state.
  N_observable is the number of observable peptides. 
  Sequence cleavage is done using \code{\link[OrgMassSpecR]{Digest}}.
}
\value{
  Named numeric vector of emPAI values.
}
\references{
  Ishihama Y, et al., Proteomics (2005)
}
\author{
  Florian P Breitwieser
}

\seealso{
  \code{\link[OrgMassSpecR]{Digest}},
  \code{\link{proteinInfo}},
  \code{\link{getProteinInfoFromUniprot}},
  \code{\link{calculate.dNSAF}},
  \code{\link{ProteinGroup}}
}
\examples{
data(ibspiked_set1)
protein.group <- proteinGroup(ibspiked_set1)
calculate.emPAI(protein.group,protein.g=protein.g(protein.group,"CERU"))
}
\keyword{ ~emPAI }
