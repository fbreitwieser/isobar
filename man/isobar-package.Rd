\name{isobar-package}
\alias{isobar-package}
\alias{isobar}
\docType{package}
\title{
Analysis and quantitation of isobarically tagged MSMS proteomics data
}
\description{
\code{isobar} provides methods for preprocessing, normalization, and
report generation for the analysis of quantitative mass spectrometry
proteomics data labeled withOA isobaric tags, such as iTRAQ and TMT.
}
\details{
\tabular{ll}{
Package: \tab isobar\cr
Version: \tab 1.1.2\cr
biocViews: \tab Proteomics, MassSpectrometray, Bioinformatics, MultipleComparisons, QualityControl\cr
Depends: \tab R (>= 2.9.0), Biobase, stats, methods, ggplot2\cr
Imports: \tab distr, biomaRt\cr
Suggests: \tab MSnbase,XML\cr
LazyLoad: \tab yes\cr
License: \tab LGPL-2\cr
URL: \tab http://bioinformatics.cemm.oeaw.ac.at\cr
Collate: \tab utils.R
ProteinGroup-class.R
IBSpectra-class.R
NoiseModel-class.R
ratio-methods.R
sharedpep-methods.R
MSnSet-methods.R
zzz.R\cr
}

Index:
\preformatted{
IBSpectra-class         IBSpectra objects
NoiseModel-class        NoiseModel objects
ProteinGroup-class      ProteinGroup objects
do.log                  Log functions for IBSpectra objects
fitCauchy               Fit weighted and unweighted Cauchy and Normal
                        distributions
groupMemberPeptides     Peptide info for protein group members
human.protein.names     Info on proteins
ibspiked_set1           Isobar Data packages
isobar-analysis         IBSpectra analysis: Protein and peptide ratio
                        calculation
isobar-import           Loading data into IBSpectra objects using
                        readIBSpectra
isobar-package          Analysis and quantitation of isobaric tag
                        Proteomics data
isobar-plots            IBSpectra plots
isobar-preprocessing    IBSpectra preprocessing
isobar-reports          Isobar reports
maplot.protein          MAplot for individual proteins
number.ranges           Helper function to transform number lists to
                        ranges
proteinInfo-methods     Methods for Function proteinInfo
proteinRatios           protein and peptide ratios
sanitize                Helper function for LaTeX export
shared.ratios           Shared ratio calculation
shared.ratios.sign      Plot and get significantly shared ratios.
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{isobar} \tab Isobar Overview (source, pdf)\cr
\code{isobar-devel} \tab Isobar for developers (source, pdf)\cr
}
}
\author{
Florian P Breitwieser <fbreitwieser@cemm.oeaw.ac.at> and
Jacques Colinge <jcolinge@cemm.oeaw.ac.at>, with contributions from
Xavier Robin <xavier.robin@unige.ch>

Maintainer: Florian P Breitwieser <fbreitwieser@cemm.oeaw.ac.at>
}
\keyword{ package }
