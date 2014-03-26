\name{Protein and peptide ratio calculation and summarization}
\alias{proteinRatios}
\alias{peptideRatios}
\alias{peptideRatiosNotQuant}
\alias{combn.matrix}
\alias{combn.protein.tbl}
\alias{summarize.ratios}
\alias{weightedMean}
\alias{weightedVariance}
\alias{weightedVariance,numeric,numeric,missing-method}
\alias{weightedVariance,numeric,numeric,numeric-method}
\alias{weightedMean,numeric,numeric-method}

\title{Calculating and Summarizing Protein and Peptide Ratios}
\description{
  A set of functions to create ratios within groups and summarize them. 
  \code{proteinRatios} serves as hub and calls \code{combn.matrix}, 
  \code{combn.protein.tbl} and \code{summarize.ratios} successively. It can be 
  used to calculate intra-class and inter-class ratios, to assess ratios and
  variability within and over cases.
}
\usage{
proteinRatios(ibspectra, noise.model, reporterTagNames = NULL,
                 proteins = reporterProteins(proteinGroup(ibspectra)),
                 peptide = NULL, cl = classLabels(ibspectra),
                 combn.method = "global", combn.vs = NULL, 
                 symmetry = FALSE, summarize =
                 FALSE, summarize.method = "mult.pval", min.detect =
                 NULL, strict.sample.pval = TRUE, strict.ratio.pval =
                 TRUE, orient.div = 0, sign.level = 0.05,
                 sign.level.rat = sign.level, sign.level.sample =
                 sign.level, ratiodistr = NULL, zscore.threshold =
                 NULL, variance.function = "maxi", combine = FALSE,
                 p.adjust = NULL, reverse = FALSE, cmbn = NULL,
                 before.summarize.f = NULL, ...)

peptideRatiosNotQuant(ibspectra, ..., peptide = unique(fData(ibspectra)[!fData(ibspectra)[["use.for.quant"]], c("peptide", "modif", "site.probs")]))
peptideRatios(ibspectra, ..., peptide = peptides(proteinGroup(ibspectra), columns = c("peptide", "modif")))
              
combn.matrix(x, method = "global", cl = NULL, vs = NULL)

combn.protein.tbl(cmbn, reverse = FALSE, ...)

summarize.ratios(ratios, by.column = "ac", summarize.method = "mult.pval",
                 min.detect = NULL, n.combination = NULL, 
                 strict.sample.pval = TRUE, strict.ratio.pval = TRUE, 
                 orient.div = 0, sign.level = 0.05, sign.level.rat = 
                 sign.level, sign.level.sample = sign.level, 
                 variance.function = "maxi", ratiodistr = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ibspectra}{IBSpectra object}
  \item{x}{for combn.matrix: reporter names. See reporterTagNames. argument of proteinRatios.}
  \item{ratios}{result of combn.protein.tbl}
  \item{by.column}{Column(s) which are the identifiers. Usually 'ac', 'peptide' or c('peptide','modif') }
  \item{cmbn}{result of combn.matrix}
  \item{before.summarize.f}{Function which is called after calculating ratios before summarizing them.}
  \item{noise.model}{NoiseModel for spectra variances}
  \item{reporterTagNames}{Reporter tags to use. By default all reporterTagNames of ibspectra object.}
  \item{proteins}{proteins for which ratios are calculated - defaults to all 
     proteins with peptides specific to them.}
  \item{peptide}{peptides for which ratios are calculated.}
  \item{cl}{Class labels. See also ?classLabels. }
  \item{vs}{Class label or reporter tag name. When \code{combn.method} is "versus.class", 
            all combinations against class \code{vs} are computed, 
            when \code{combn.method} is "verus.channel", all combinations 
            against channel \code{vs}.}
  \item{combn.method}{"global", "interclass", "intra-class", "versus.class" or "versus.channel". Defines which ratios
     are computed, based on class labels cl}
  \item{method}{See combn.method}
  \item{combn.vs}{vs argument for combn, if combn.method is "versus.class" or "versus.channel".}
  \item{symmetry}{If true, reports also the inverse ratio}
  \item{summarize}{If true, ratios for each protein are summarized.}
  \item{summarize.method}{"isobar", for now.}
  \item{min.detect}{How many times must a ratio for a protein be present when
     summarizing? When NULL, defaults to the maximum number of combinations.}
  \item{strict.sample.pval}{If true, missing ratios are penalized by giving 
     them a sample.pval of 0.5.}
  \item{strict.ratio.pval}{If true, take all ratios into account. If false, only
    take ratios into account which are in the same direction as the majority of 
    ratios}
  \item{orient.div}{Number of ratios which might go in the wrong direction.}
  \item{sign.level}{Significance level}
  \item{sign.level.rat}{Significance level on ratio p-value}
  \item{sign.level.sample}{Significance level on sample p-value}
  \item{ratiodistr}{Protein ratio distribution}
  \item{variance.function}{Variance function}
  \item{zscore.threshold}{z-score threshold to apply}
  \item{\dots}{Passed to estimateRatio()}
  \item{combine}{If true, a single ratio for all proteins and peptides, resp., is calculated. See \code{\link{estimateRatio}}.}
  \item{p.adjust}{Set to one of p.adjust.methods to adjust ratio p-values for multiple comparisions. See \code{\link{p.adjust}}.}
  \item{reverse}{reverse}
  \item{n.combination}{number of combinations possible}
}
\value{
'data.frame':	11 variables:
\item{lratio }{log ratio}
\item{variance }{variance}
\item{n.spectra }{Number of spectra used for quantification}
\item{p.value.rat }{Signal p-value (NA if ratiodistr is missing)}
\item{p.value.sample }{Sample p-value (NA if ratiodistr is missing)}
\item{is.significant }{Is the ratio significant? (NA if ratiodistr is
  missing)}
\item{protein }{Protein quantified}
\item{r1 }{r1}
\item{r2 }{r2}
}
\author{
  Florian P Breitwieser, Jacques Colinge
}
\seealso{
  \link{IBSpectra},
  \link{isobar-preprocessing}
  \link{isobar-analysis}
}
\examples{

combn.matrix(114:117,method="interclass",cl=as.character(c(1,1,2,2)))
combn.matrix(114:117,method="interclass",cl=as.character(c(1,1,2,2)))
combn.matrix(114:117,method="global")

data(ibspiked_set1)
data(noise.model.hcd)

ceru.proteins <- c("P13635","Q61147")
proteinRatios(ibspiked_set1,noise.model=noise.model.hcd,proteins=ceru.proteins,cl=c("T","T","C","C"),combn.method="interclass",summarize=TRUE)

}
