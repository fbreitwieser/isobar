\name{proteinInfo-methods}
\alias{getProteinInfoFromBiomart}
\alias{getProteinInfoFromUniprot}
\alias{getProteinInfoFromBioDb}
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
  \code{proteinInfo} method allows to get and set it.

  \code{getProteinInfoFromUniprot} downloads information of contained proteins
  from Uniprot, \code{getProteinInfoFromBiomart} from Biomart.
}
\usage{

  \S4method{proteinInfo}{ProteinGroup}(x)

  \S4method{proteinInfo}{ProteinGroup,character}(x, protein.g, select="name", collapse=", ",
                                                 simplify = TRUE, do.warn = TRUE)

  getProteinInfoFromUniprot(x, splice.by = 200)

  getProteinInfoFromBiomart(x, database = "Uniprot")

  getProteinInfoFromBioDb(x, ..., con = NULL)
}
\arguments{
  \item{x}{ProteinGroup object}
  \item{protein.g}{Protein group identifier. If supplied, only information for 
    these proteins is returned.}
  \item{select}{indicating columns to select. See Details.}
  \item{collapse}{passed to \code{\link{paste}} to concatenate information of
    multiple protein in one protein group.}
  \item{simplify}{If true, a vector or matrix is returned, with the pasted protein information.
                  If false, a list is returned.}
  \item{do.warn}{If true, report diagnostic warning messages.}
  \item{splice.by}{Chunk size for query of Uniprot database.}
  \item{database}{database from which the ACs stem from. Only Uniprot is supported for now.}
  \item{con}{database connection}
  \item{...}{arguments to build database connection.}
}

\details{proteinInfo contains columns \code{accession}, \code{name}, 
  \code{gene_name}, \code{protein_name}, and possibly \code{length}
  and \code{sequence}. \code{accession} is mapped with
  the entry AC is mapped to the entry AC in the database.
  \code{getProteinInfoFromUniprot} is the preferred methods to get the information.
  \code{getProteinInfoFromBioDb} is an example how to implement the query on a
  local database.
}
\seealso{
  \link{protein.g}
}
\examples{
data(ibspiked_set1)
pg <- proteinGroup(ibspiked_set1)

\dontrun{
  proteinInfo(pg) <- getProteinInfoFromUniprot(pg)
  proteinInfo(pg) <- getProteinInfoFromBiomart(pg)
}

proteinInfo(pg,protein.g="P13635")
protein.g(pg,"CERU")
}
\keyword{methods}
