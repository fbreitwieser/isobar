##
## Isobar properties.R file
##     for automatic report generation
## 
## It is standard R code and parsed using sys.source

#####################################################################
## General properties

## Report type: Either 'protein' or 'peptide'
# report.level="peptide"
report.level="protein"
#attr(report.level,"allowed.values") <- c("protein","peptide")

## Isobaric tagging type. Use one of the following:
# type='iTRAQ4plexSpectra'
# type='iTRAQ8plexSpectra'
# type='TMT2plexSpectra'
# type='TMT6plexSpectra'
type=NULL
#attr(type,"allowed.values") <- IBSpectraTypes()

isotope.impurities=NULL
correct.isotope.impurities=TRUE

## Name of project, by default the name of working directory
## Will be title and author of the analysis reports.
name=basename(getwd())
author=paste0("isobar R package v",packageDescription("isobar")$Version)

## specifes the IBSpectra file or object
##   - can be a data.frame (e.g. ibspectra=as.data.frame(ibspiked_set1) )
##   - if it is a character string, it is assumed to be a file
##     - if it ends on .rda, then it is assumed to be a R data object
##     - if it does not exists, then it is may generated based on 
##        the peaklist and identifications properties
ibspectra=paste(name,"ibspectra.csv",sep=".")

## When replicates or 'samples belonging together' are analyzed, a
## ProteinGroup object based on all data should be constructed
## beforehand. This then acts as a template and a subset is used.
protein.group.template=NULL

## Via database or internet connection, informations on proteins (such
## as gene names and length) can be gathered. protein.info.f defines
## the function which takes a ProteinGroup object as argument
protein.info.f=getProteinInfoFromTheInternet

## Where should cached files be saved? Will be created if it does not
## exist
# cachedir="."
cachedir="cache"
## Regenerate cache files? By default, chache files are used.
regen=FALSE

## An ibspectra object can be generated from peaklists and
## identifications.

## peaklist files for quantitation, by default all mgf file in
## directory
peaklist=list.files(pattern="*\\.mgf$")
## id files, by default all id.csv files in directory
identifications=list.files(pattern="*\\.id.csv$")
## mapping files, for data quantified and identified with different but
## correspoding spectra. For example corresponding HCD-CID files.

## masses and intensities which are outside of the 'true' tag mass 
## +/- fragment.precision/2 are discarded
fragment.precision=0.01
## filter mass outliers
fragment.outlier.prob=0.001

## Additional arguments of readIBSpectra can be set here
##  decode.titles should be set to TRUE for Mascot search results
##  as Mascot encodes the spectrum title (e.g. space -> %20)
readIBSpectra.args = list(
    mapping.file=NULL,
    decode.titles=FALSE
)

#####################################################################
## Quantification properties

normalize=TRUE
# if defined, normalize.factors will be used for normalization
normalize.factors=NULL 
normalize.channels=NULL
normalize.use.protein=NULL
normalize.exclude.protein=NULL
normalize.function=median
normalize.na.rm=FALSE

peptide.specificity=REPORTERSPECIFIC

use.na=FALSE

## the parameter noise.model can be either a NoiseModel object or a file name
data(noise.model.hcd)
noise.model=noise.model.hcd
## If it is a file name, a noise model is estimated as non one-to-one
## and saved into the file. otherwise, the noise model is loaded from
## the file
# noise.model="noise.model.rda"

## Define channels for creation of a noise model, ideally a set of 
##  channels which are technical replicates. 
noise.model.channels=NULL

## If noise.model.is.technicalreplicates is FALSE, the intensities
##   are normalized for protein means, creating artifical technical
##   replicates. For this procedure, only proteins with more than
##   noise.model.minspectra are considered.
noise.model.is.technicalreplicates=FALSE
noise.model.minspectra=50

## class labels. Must by of type character and of same length as
## number of channels I. e. 4 for iTRAQ 4plex, 6 for TMT 6plex Example
## for iTRAQ 4plex:
## Class definitions of the isobaric tag channels. 
## A character vector with the same length as channels 
##   (e.g. 4 for iTRAQ 4plex, 6 for TMT 6plex)
## Example for iTRAQ 4plex:
# class.labels=as.character(c(1,0,0,0))
# class.labels=c("Treatment","Treatment","Control","Control")
## Also names are possible - these serves as description in the report
##  and less space is used in the rows
# class.labels=c("Treatment"="T","Treatment"="T","Control"="C","Control"="C")
class.labels=NULL

## The following parameters define which ratios are calculated.

## summarize ratios with equal class labels, set to TRUE when replicates are used
summarize=FALSE

## combn.method defines which ratios are calculated - versus a channel or a class,
##   all the ratios within or across classes, or all possible combinatioins.
## When summarize=TRUE is set, use "interclass", "versus.class", or "intraclass"
# combn.method="global"
# combn.method="versus.class"
# combn.method="intraclass"
# combn.method="interclass"
combn.method="versus.channel"
vs.class=NULL

cmbn=NULL

## Arguments given to 'proteinRatios' function. See ?proteinRatios
ratios.opts = list(
    sign.level.sample=0.05,
    sign.level.rat=0.05,
    groupspecific.if.same.ac=TRUE)

quant.w.grouppeptides=c()

min.detect=NULL

preselected=c()

### Biological Variability Ratio Distribution options
## ratiodistr can be set to a file or a 'Distribution object. ' If
##  NULL, or the specified file is not existent, the biological
##  variability of ratios is estimated on the sample at hand and
##  written to cachedir/ratiodistr.rda or the specified file.
ratiodistr=NULL

## Ideally, when the biological variability is estimated for the
## sample at hand, a biological replicate is present (/ie/ same class
## defined in class labels).  Classes can also be assigned just for
## estimation of the ratio distribution, /eg/ to choose biologically
## very similar samples as pseudo replicates.
ratiodistr.class.labels=NULL

## Function for fitting. Available: fitCauchy, fitTlsd
ratiodistr.fitting.f=fitCauchy

## Use symetrical ratios - i.e. for every ratio r add a ratio -r
##   prior to fitting of a distribution
ratiodistr.symmetry=TRUE

## If defined, use z-score instead of ratio distribution
# zscore.threshold=2.5
zscore.threshold=NULL

####################################################################
## PTM properties

## PhosphoSitePlus dataset which can be used to annotate known
## modification sites. Download site:
## http://www.phosphosite.org/staticDownloads.do
phosphosite.dataset <- NULL

## Modification to track. Use 'PHOS' for phosphorylation.
# ptm <- c('ACET','METH','UBI','SUMO', 'PHOS')
ptm <- NULL

## file name of rda or data.frame with known modification sites
## gathered with ptm.info.f.  defaults to 'cachedir/ptm.info.rda'
ptm.info <- NULL

## Function to get PTM modification sites from public datasets
# ptm.info.f <- getPtmInfoFromNextprot
# ptm.info.f <- function(...)
#    getPtmInfoFromPhosphoSitePlus(...,modification="PHOS")
# ptm.info.f <- function(...)
#    getPtmInfoFromPhosphoSitePlus(...,modification=ptm)
ptm.info.f <- getPtmInfoFromNextprot

## A protein quantification data.frame (generated with
##  'proteinRatios').  The ratio and variance are used to correct the
##  observed modified peptide ratios Needs to have the experimental
##  setup as the modified peptide experiment
correct.peptide.ratios.with <- NULL
## Protein groups to use with correct.peptide.ratios()
correct.peptide.ratios.with_protein.group <- NULL

## The correlation between peptide and protein ratios defines the
## covariance

##  Var(ratio m) = Var(ratio mp) + Var(ratio p)
##                                    + 2 * Cov(ratio mp, ratio p),
##  Cov(ratio mp, ratio p) = 2 * cor * Sd(ratio mp) * Sd(ratio p),
##  with m = modifcation, mp = modified peptide, p = protein
peptide.protein.correlation <- 0

## quantification table whose columns are attached to the XLS
## quantification table
compare.to.quant <- NULL

#####################################################################
## Report properties 

write.qc.report=TRUE
write.report=TRUE
write.xls.report=TRUE

## Use name for report, ie NAME.quant.xlsx instead of
## isobar-analysis.xlsx
use.name.for.report=TRUE

## PDF Analysis report sections: Significant proteins and protein
## details
show.significant.proteins=FALSE
show.protein.details=TRUE

### QC REPORT OPTIONS ###
#qc.maplot.pairs=FALSE # plot one MA plot per tag (versus all others)
qc.maplot.pairs=TRUE # plot MA plot of each tag versus each tag

### XLS REPORT OPTIONS ###
## Spreadsheet format: Either 'xlsx' or 'xls'
# spreadsheet.format="xlsx"
spreadsheet.format="xlsx"

## XLS report format 'wide' or 'long '.

## 'wide' format outputs ratios in separate columns of the same record
## (i.e. one line per protein)
## 'long' format outputs ratios in separate records (i.e. one line per
## ratio)
# xls.report.format="wide"
xls.report.format="long"

## XLS report columns in quantification tab
##  possible values: ratio, is.significant, CI95.lower, CI95.upper, 
##                   ratio.minus.sd, ratio.plus.sd,
##                   p.value.ratio, p.value.sample, n.na1, n.na2, 
##                   log10.ratio, log10.variance,
##                   log2.ratio, log2.variance
##  only for summarize=TRUE: n.pos, n.neg
xls.report.columns <- c("ratio","is.significant","ratio.minus.sd",
                        "ratio.plus.sd","p.value.ratio","p.value.sample",
                        "log10.ratio","log10.variance")

## Perl command to be used for Excel report generation
# perl.cmd = "C:/Strawberry/perl/bin/perl.exe"
# perl.cmd = "C:/Perl/bin/perl.exe"
# perl.cmd = "perl5"
perl.cmd = "perl"

#####################################################################
## Etc

sum.intensities=FALSE

datbase="Uniprot"

scratch=list()

##
#  compile LaTeX reports into PDF files
compile=TRUE

# zip final report files into archive
zip=FALSE

# warning level (see 'warn' in ?options)
warning.level=1


###########################################################
## Novel options

shrink.mean=TRUE
use.t.stat=TRUE
