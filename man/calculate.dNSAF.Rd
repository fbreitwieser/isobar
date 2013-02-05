\name{calculate.dNSAF}
\alias{calculate.dNSAF}

\title{
  dNSAF approximate abundance calculations.
}
\description{
  Distributed normalized spectral abundance factor (dNSAF) is a 
  label free quantitative measure of protein abundance based on 
  spectral counts which are corrected for peptides shared by multiple
  proteins. Original publication: Zhang Y et al., Analytical Chemistry
  (2010).
}
\usage{
calculate.dNSAF(protein.group, use.mw = FALSE, normalize = TRUE,
                combine.f = mean)
}
\arguments{
  \item{protein.group}{ProteinGroup object. Its \code{@proteinInfo}
    slot \code{data.frame} must contain a \code{length} column.}
  \item{use.mw}{Use MW to account for protein size}
  \item{normalize}{Normalize dSAF to dNSAF?}
  \item{combine.f}{How to handle proteins seen only with shared peptides?}
}
\value{
  Named numeric vector of dNSAF values.
}
\references{
  Zhang Y et al., Analytical Chemistry (2010)
}
\author{
  Florian P Breitwieser
}

\seealso{
  \code{\link{proteinInfo}},
  \code{\link{getProteinInfoFromUniprot}},
  \code{\link{calculate.emPAI}},
  \code{\link{ProteinGroup}}
}

\examples{
data(ibspiked_set1)
protein.group <- proteinGroup(ibspiked_set1)
calculate.dNSAF(protein.group)
}

\keyword{ ~dNSAF }
