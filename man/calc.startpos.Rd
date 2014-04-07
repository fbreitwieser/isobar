\name{calcPeptidePosition}
\alias{calcPeptidePosition}
\title{Recalculate peptide start positions based on protein sequence
}
\description{Function to recalculate start position of peptide in protein when it is missing or wrong.
}
\usage{
calcPeptidePosition(peptide.info, protein.info, calc.il.peptide)
}
\arguments{
  \item{peptide.info}{Peptide info object of ProteinGroup. }
  \item{protein.info}{Protein info object of ProteinGroup. }
  \item{calc.il.peptide}{Should the 'real' peptide (I/L difference) be calculated?}
}
