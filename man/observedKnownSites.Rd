\name{observedKnownSites}
\alias{observedKnownSites}
\alias{modif.site.count}
\alias{modif.sites}

\title{
  Observed modification sites.
}
\description{
  Functions to display the modification sites observed for each protein isoform
  and count the number of modified residues per protein.
}
\usage{
observedKnownSites(protein.group, protein.g, ptm.info, modif, modification.name = NULL)

modif.site.count(protein.group, protein.g = reporterProteins(protein.group), modif, take = max)

modif.sites(protein.group, protein.g = reporterProteins(protein.group), modif)

}

\arguments{
  \item{protein.group}{ProteinGroupb object.}
  \item{protein.g}{protein group identifier.}
  \item{ptm.info}{ptm information data.frame, see ?getPtmInfo.}
  \item{modif}{Modification to track, e.g. 'PHOS'.}
  \item{modification.name}{Value to filter 'modification.name' column in ptm.info.}
  \item{take}{should be either max or min: When multiple isoforms are present, which value should be taken for the count?}
}
\author{
  Florian P. Breitwieser
}
\examples{
data(ib_phospho)
data(ptm.info)

# Modification sites of reporter proteins:
#  a list of protein groups, 
#     containing sub-lists of identified sites for each isoform
protein.modif.sites <- sort(modif.site.count(proteinGroup(ib_phospho),modif="PHOS"))

# Details on modification sites of proteins
#   detected with most modifications
modif.sites(proteinGroup(ib_phospho),modif="PHOS",protein.g=names(tail(protein.modif.sites)))

# How many sites are known, and how many known sites have been observed?
observedKnownSites(proteinGroup(ib_phospho),modif="PHOS",protein.g=names(tail(protein.modif.sites)),ptm.info=ptm.info,modification.name="Phospho")
}
