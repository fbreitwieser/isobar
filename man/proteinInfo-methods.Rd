\name{proteinInfo-methods}
\alias{getProteinInfoFromBiomart}
\alias{proteinInfo}
\alias{proteinInfo-methods}
\alias{proteinInfo,ProteinGroup-method}
\alias{proteinInfo,ProteinGroup,character-method}
\alias{proteinInfo,ProteinGroup,missing-method}
\alias{proteinInfo<-}
\alias{proteinInfo<-,ProteinGroup-method}

\title{Methods for Function proteinInfo}
\description{
  proteinInfo slot in Proteingroup objects contains information about proteins. 
  \code{proteinInfo} method allows to get and set it, 
  \code{getProteinInfoFromBiomart} downloads information of contained proteins
  from Biomart.
}
\usage{

  \S4method{proteinInfo}{ProteinGroup}(x)

  \S4method{proteinInfo}{ProteinGroup,character}(x, protein.g, select="name", collapse=", ")

  getProteinInfoFromBiomart(x, database = "Uniprot")
}
\arguments{
  \item{x}{ProteinGroup object}
  \item{protein.g}{Protein group identifier. If supplied, only information for 
    these proteins is returned.}
  \item{select}{indicating columns to select. See Details.}
  \item{collapse}{passed to \code{\link{paste}} to concatenate information of
    multiple protein in one protein group}
  \item{database}{database from which the ACs stem from.}
}

\details{proteinInfo contains columns \code{accession}, \code{name}, 
  \code{gene_name}, and \code{protein_name}. \code{accession} is mapped with
  the entry AC is mapped to the entry AC in the database.
}
\seealso{
  \link{protein.g}
}
\examples{
data(ibspiked_set1)
pg <- proteinGroup(ibspiked_set1)

\dontrun{
  proteinInfo(pg) <- getProteinInfoFromBiomart(pg)
}

proteinInfo(pg,protein.g="P13635")
protein.g(pg,"CERU")
}
\keyword{methods}
