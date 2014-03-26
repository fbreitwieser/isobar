\name{getPeptideModifContext}
\alias{getPeptideModifContext}
\title{Get context of modification}
\description{Gets neighboring amino acids around modification which can be used to find enriched motifs.
}
\usage{
getPeptideModifContext(protein.group, modif, n.aa.up = 7, n.aa.down = 7)
}
\arguments{
  \item{protein.group}{ProteinGroup object.}
  \item{modif}{Modification of interest.}
  \item{n.aa.up}{Number of AA downstream to report.}
  \item{n.aa.down}{Number of AA upstream to report.}
}
