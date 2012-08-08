\name{proteinNameAndDescription}
\alias{proteinNameAndDescription}
\alias{proteinDescription}
\alias{proteinGeneName}
\alias{proteinID}

\title{Get protein gene names and description from protein info of protein group.}
\description{
  Convenience functions to retrieve protein gene names and description 
  for a list of protein group identifiers.
}
\usage{
proteinNameAndDescription(protein.group, protein.g = reporterProteins(protein.group), collapse = FALSE)
proteinGeneName(protein.group, protein.g = reporterProteins(protein.group))
proteinDescription(protein.group, protein.g = reporterProteins(protein.group))
proteinID(protein.group, protein.g = reporterProteins(protein.group))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein.group}{ProteinGroup object.}
  \item{protein.g}{protein group identifier.}
  \item{collapse}{If TRUE, the information for all protein.gs is combined.}
}
\author{Florian P Breitwieser}

\examples{
data(ibspiked_set1)
pg <- proteinGroup(ibspiked_set1)
protein.gs <- protein.g(pg,"CERU")
protein.gs
proteinNameAndDescription(pg,protein.gs)
proteinNameAndDescription(pg,protein.gs,collapse=TRUE)
proteinGeneName(pg,protein.gs)
proteinDescription(pg,protein.gs)  
proteinID(pg,protein.gs)  
}
