\name{getPtmInfo}
\alias{getPtmInfoFromPhosphoSitePlus}
\alias{getPtmInfoFromNextprot}

\title{
  Get PTM site information for idenfied proteins from public databases.
}
\description{
  Get PTM site information for idenfied proteins from public databases.
}
\usage{
  getPtmInfoFromPhosphoSitePlus(protein.group, file.name = NULL, modif = "PHOS", 
                                psp.url = "http://www.phosphosite.org/downloads/", 
                                mapping = c(PHOS = "Phosphorylation_site_dataset.gz", 
                                            ACET = "Acetylation_site_dataset.gz", 
                                            METH = "Methylation_site_dataset.gz", 
                                            SUMO = "Sumoylation_site_dataset.gz", 
                                            UBI = "Ubiquitination_site_dataset.gz"))

  getPtmInfoFromNextprot(protein.group, nextprot.url = "http://www.nextprot.org/rest/entry/NX_XXX/ptm?format=json",
                         url.wildcard = "XXX")
}
\arguments{
  \item{protein.group}{ProteinGroup object.}
  \item{file.name}{File name to save downloaded data, defaults to the original file name (see mapping).}
  \item{modif}{Selects dataset to download (see mapping).}
  \item{psp.url}{PhosphoSitePlus main URL for datasets.}
  \item{mapping}{Names of PhosphoSitePlus modification datasets, mapped by modif name.}
  \item{nextprot.url}{URL for fetching Nextprot results. \code{url.wildcard} will be replaced by the Uniprot Protein AC.}
  \item{url.wildcard}{wildcard to replace with Uniprot protein AC in \code{nextprot.url}.}
}
\details{PhosphoSitePlus datasets are downloaded and written to the working directory with its original name (see mapping) 
  unless a file with that name exists, which is then parsed into a data.frame of suitable format.
}
\value{
  data.frame with (at least) the columns:
    isoform_ac, description, evidence, position
}
\references{
  PhosphoSitePlus: a comprehensive resource for investigating the structure and function of experimentally determined post-translational modifications in man and mouse.
  Hornbeck PV, Kornhauser JM, Tkachev S, Zhang B, Skrzypek E, Murray B, Latham V, Sullivan M.
  Nucleic Acids Res. 2012 Jan;40(Database issue):D261-70. Epub 2011 Dec 1.

  neXtProt: a knowledge platform for human proteins.
  Lane L, Argoud-Puy G, Britan A, Cusin I, Duek PD, Evalet O, Gateau A, Gaudet P, Gleizes A, Masselot A, Zwahlen C, Bairoch A.
  Nucleic Acids Res. 2012 Jan;40(Database issue):D76-83. Epub 2011 Dec 1.
}
\author{
  Florian P. Breitwieser
}
\note{
  PhosphoSitePlus is licensed under Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License and is freely available for non-commercial purpose,
  see http://www.phosphosite.org/staticDownloads.do.

  neXtProt is licensed under the Creative Commons Attribution-NoDerivs License, see: http://creativecommons.org/licenses/by-nd/3.0.

  Please read the conditions and use the data only if you agree.
}

\examples{

\dontrun{
  data(ib_phospho)
  ptm.info.np <- getPtmInfoFromNextprot(proteinGroup(ib_phospho))
  ptm.info.np <- ptm.info.np[grep("Phospho",ptm.info.np$modification.name),]
  ptm.info.psp <- getPtmInfoFromPhosphoSitePlus(proteinGroup(ib_phospho),modif="PHOS")

  str(ptm.info.np)
  str(ptm.info.psp)

  }
}
