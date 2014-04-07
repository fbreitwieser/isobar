\name{proteinInfo-methods}
\alias{getProteinInfoFromBiomart}
\alias{getProteinInfoFromUniprot}
\alias{getProteinInfoFromTheInternet}
\alias{getProteinInfoFromNextProt}
\alias{getProteinInfoFromEntrez}
\alias{getProteinInfoFromBioDb}
\alias{proteinInfoIsOnSpliceVariants}
\alias{proteinInfo}
\alias{proteinInfo-methods}
\alias{proteinInfo,ProteinGroup-method}
\alias{proteinInfo,ProteinGroup,missing,missing-method}
\alias{proteinInfo,ProteinGroup,character,missing-method}
\alias{proteinInfo,ProteinGroup,missing,character-method}
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

  \S4method{proteinInfo}{ProteinGroup,character,missing}(x, protein.g, select="name", collapse=", ",
                                                 simplify = TRUE, do.warn = TRUE)

  \S4method{proteinInfo}{ProteinGroup,missing,character}(x, protein.ac, select="name", collapse=", ",
                                                 simplify = TRUE, do.warn = TRUE)


  proteinInfoIsOnSpliceVariants(protein.info)

#  getProteinInfoFromUniprot(x, splice.by = 200, fields = c(accession = "id", name
#                 = "entry%20name", protein_name = "protein%20names",
#                 gene_name = "genes", organism = "organism", length =
#                 "length", sequence = "sequence"))


  getProteinInfoFromTheInternet(x)

  getProteinInfoFromNextProt(x)

  getProteinInfoFromBiomart(x, database = "Uniprot")

  getProteinInfoFromBioDb(x, ..., con = NULL)

  getProteinInfoFromEntrez(x, splice.by = 200)
}
\arguments{
  \item{x}{ProteinGroup object}
  \item{protein.g}{Protein group identifier. If supplied, only information for 
    these proteins is returned.}
  \item{protein.ac}{Protein ACs. If supplied, only information for 
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
  \item{fields}{mapping of CSV field names to proteinInfo field names}
  \item{...}{arguments to build database connection.}
  \item{protein.info}{protein info data.frame}
}

\details{proteinInfo contains columns \code{accession}, \code{name}, 
  \code{gene_name}, \code{protein_name}, and possibly \code{length}
  and \code{sequence}. \code{accession} is mapped with
  the entry AC is mapped to the entry AC in the database.
  \code{getProteinInfoFromUniprot} is the preferred methods to get the information.
  \code{getProteinInfoFromBioDb} is an example how to implement the query on a
  local database. Depending on the database, protein information might be
  available on protein ACs or also on the specific splice variants. This can be
  queried with the \code{proteinInfoIsOnSpliceVariants} function.
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
