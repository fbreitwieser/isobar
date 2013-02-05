\name{isobar-analysis}

\alias{isobar-analysis}

% estimateRatio methods
\alias{estimateRatio}
\alias{estimateRatio,IBSpectra,ANY,character,character,character,missing-method}
\alias{estimateRatio,IBSpectra,ANY,missing,missing,missing,character-method}
\alias{estimateRatio,IBSpectra,ANY,missing,missing,character,missing-method}
\alias{estimateRatio,IBSpectra,ANY,character,character,missing,matrix-method}
\alias{estimateRatio,IBSpectra,ANY,character,character,NULL,matrix-method}
\alias{estimateRatio,IBSpectra,ANY,character,character,character,NULL-method}
\alias{estimateRatio,IBSpectra,ANY,character,character,missing,character-method}
\alias{estimateRatio,IBSpectra,ANY,character,character,NULL,character-method}
\alias{estimateRatio,IBSpectra,ANY,character,character,NULL,data.frame-method}
\alias{estimateRatio,IBSpectra,ANY,character,character,missing,data.frame-method}
\alias{estimateRatioNumeric}
\alias{estimateRatioNumeric,numeric,numeric,missing-method}
\alias{estimateRatioNumeric,numeric,numeric,NULL-method}
\alias{estimateRatioNumeric,numeric,numeric,NoiseModel-method}

\alias{estimateRatioForProtein}
\alias{estimateRatioForPeptide}

% multUnifPValues and multUnifDensity
\alias{getMultUnifPValues}
\alias{getMultUnifDensity}

\title{IBSpectra analysis: Protein and peptide ratio calculation}
\description{Calculates the relative abundance of a peptide or protein in one tag compared to another.}
\usage{
estimateRatio(ibspectra, noise.model = NULL, channel1, channel2, protein, peptide, ...)
estimateRatioForPeptide(peptide, ibspectra, noise.model, channel1, channel2, combine = TRUE, ...)
estimateRatioForProtein(protein, ibspectra, noise.model, channel1, channel2, combine = TRUE, method = "isobar", specificity = REPORTERSPECIFIC, quant.w.grouppeptides = NULL, ...)

\S4method{estimateRatioNumeric}{numeric,numeric,missing}(channel1,channel2,summarize.f=median, ...)

\S4method{estimateRatioNumeric}{numeric,numeric,NoiseModel}(channel1,channel2,noise.model,ratiodistr=NULL,variance.function="maxi",
                                                            sign.level=0.05,sign.level.rat=sign.level,sign.level.sample=sign.level,
                                                            remove.outliers=TRUE,outliers.args=list(method = "iqr", outliers.coef = 1.5),
                                                            method="isobar",fc.threshold=1.3,
                                                            channel1.raw=NULL,channel2.raw=NULL,use.na=FALSE,preweights=NULL)

\S4method{estimateRatio}{IBSpectra,ANY,character,character,character,missing}(ibspectra,noise.model,channel1,channel2,
                                                                              protein,peptide,...)

\S4method{estimateRatio}{IBSpectra,ANY,character,character,character,NULL}(ibspectra,noise.model,channel1,channel2,
                                                                           protein,peptide=NULL,...)

\S4method{estimateRatio}{IBSpectra,ANY,character,character,missing,character}(ibspectra,noise.model,channel1,channel2,protein,peptide,...)
\S4method{estimateRatio}{IBSpectra,ANY,character,character,NULL,character}(ibspectra,noise.model,channel1,channel2,protein=NULL,peptide,...)

}
\arguments{
  \item{ibspectra}{ IBSpectra object. } 
  \item{noise.model}{ NoiseModel object.}
  \item{channel1}{ Tag channel 1. 
    Can either be a character denoting a 'reporter name' or a numeric vector
    whose value should be summarized.Ratio is calculated as channel2/channel1. }
  \item{channel2}{ Tag channel 2.
    Can either be a character denoting a 'reporter name' or a numeric vector
    whose value should be summarized. Ratio is calculated as channel2/channel1. }
  \item{protein}{ Protein(s) of interest. If present, channel1 and channel2 
    must be reporter names. Provide either proteins or peptides.  }
  \item{peptide}{ Peptide(s) of interest. If present, channel1 and channel2 
    must be reporter names. Provide either proteins or peptides. }
  \item{combine}{ If true, a single ratio is returned even for multiple
    peptides/spectra. If false, a data.frame with a row for each peptide/protein
    is returned.}
  \item{specificity}{ See \code{\link{specificities}}. }
  \item{quant.w.grouppeptides}{Proteins which should be quantified with group specific peptides.
        Normally, only reporter specific peptides are used.}
  \item{ratiodistr}{\code{distr} object of ratio distribution.}
  \item{variance.function}{Defines how the variance for ratio is calculated. 
    'ev' is the estimator variance and thus 1/sum(1/variances).
    'wsv' is the weighted sample variance.
    'maxi' method takes the maximum of the former two variances.
  }
  \item{sign.level}{Significiance level. }
  \item{sign.level.rat}{Signal p-value significiance level. }
  \item{sign.level.sample}{Sample p-value significiance level. }
  \item{remove.outliers}{Should outliers be removed? }
  \item{outliers.args}{Arguments for outlier removal, see OUTLIERS function (TODO).}
  \item{method}{ method taken for ratio computation and selection: one of
    'isobar','libra','multiq','pep','ttest' and 'compare.all'. }
  \item{fc.threshold}{When method equals fc, takes this as fold change threshold. }
  \item{summarize.f}{A method for summarizing spectrum ratios when no other 
    information is available. For example \code{median} or \code{mean}. }
  \item{channel1.raw}{When given, noise estimation is based on channel1.raw and channel2.raw. 
                      These are the intensities of the channels before normalization. }
  \item{channel2.raw}{See channel1.raw. }
  \item{use.na}{Use NA values to calculate ratio. Experimental feature - use with caution. }
  \item{preweights}{Specifies weigths for each spectrum. Experimental feature - use with caution. }
%  \item{correct.ratio}{Can be used to scale ratio - should be in log 10. }
  \item{\dots}{Passed down to \code{estimateRatioNumeric} methods.}
}
\value{
  In general, a named character vector with the following elements:
    - lratio: log ratio
    - variance
    - n.spectra: number of spectra available in the ratio calculation
    - p.value.rat: Signal p-value. NA if called w/o ratiodistr
    - p.value.sample: Sample p-value. NA if called w/o ratiodistr
    - is.significant: NA if called w/o ratiodistr

   If combine=FALSE, estimateRatio returns a data.frame, with columns as 
   described above.
}

\author{Florian P. Breitwieser, Jacques Colinge}

\seealso{
  \link{ProteinGroup},
  \link{IBSpectra},
  \link{isobar-preprocessing},
  \link{isobar-plots}
  \link{proteinRatios}
}

\examples{

  data(ibspiked_set1)
  data(noise.model.hcd)
  ceru.human <- protein.g(proteinGroup(ibspiked_set1),"CERU_HUMAN")
  ceru.rat <- protein.g(proteinGroup(ibspiked_set1),"CERU_RAT")
  ceru.mouse <- protein.g(proteinGroup(ibspiked_set1),"CERU_MOUSE")
  ceru.proteins <- c(ceru.human,ceru.rat,ceru.mouse)

## Calculate ratio based on all spectra of peptides specific 
##  to CERU_HUMAN, CERU_RAT or CERU_MOUSE. Returns a named
##  numeric vector.
10^estimateRatio(ibspiked_set1,noise.model.hcd,
                 channel1="114",channel2="115",
                 protein=ceru.proteins)['lratio']

## If argument 'combine=FALSE', estimateRatio returns a data.frame 
##  with one row per protein
10^estimateRatio(ibspiked_set1,noise.model.hcd,
                 channel1="114",channel2="115",
                 protein=ceru.proteins,combine=FALSE)[,'lratio']
## spiked material channel 115 vs 114: 
##                 CERU_HUMAN (P00450): 1
##                 CERU_RAT   (P13635): 2
##                 CERU_MOUSE (Q61147): 0.5

}
