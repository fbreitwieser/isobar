\name{peptide.count}
\alias{peptide.count}
\alias{spectra.count}
\title{
  Peptide and spectral counts for ProteinGroup objects.
}
\description{
  Reports the peptide and spectral count for supplied proteins.
}
\usage{
peptide.count(protein.group, protein.g = reporterProteins(protein.group), 
              specificity = c("reporter-specific", "group-specific", "unspecific"))
spectra.count(protein.group, protein.g = reporterProteins(protein.group), 
              specificity = c("reporter-specific", "group-specific", "unspecific"))
}
\arguments{
  \item{protein.group}{ProteinGroup object.}
  \item{protein.g}{Protein group identifier.}
  \item{specificity}{Specificity of peptides.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
\keyword{ ~kwd1 }
