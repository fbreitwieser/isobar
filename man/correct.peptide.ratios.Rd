\name{correct.peptide.ratios}
\alias{correct.peptide.ratios}

\title{Correct peptide ratios with protein ratios from a separate experiment.}
\description{Correct peptide ratios with protein ratios from a separate experiment.}

\usage{
correct.peptide.ratios(ibspectra, peptide.quant.tbl, protein.quant.tbl,
         protein.group.combined, adjust.variance = TRUE,
         correlation = 0, recalculate.pvalue = TRUE)
}

\arguments{
  \item{ibspectra}{IBSpectra object.}
  \item{peptide.quant.tbl}{Calculated with peptideRatios.}
  \item{protein.quant.tbl}{Calculated with proteinRatios.}
  \item{protein.group.combined}{ProteinGroup object generated on both PTM and protein data.}
  \item{adjust.variance}{Adjust variance of ratios.}
  \item{correlation}{Assumed correlation between peptide and protein ratios for variance adjustment.}
  \item{recalculate.pvalue}{Recalculate p-value after variance adjustment.}
}
\author{Florian P. Breitwieser}
