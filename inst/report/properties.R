##
## Isobar properties.conf file
##     for automatic report generation
## 
## It is standard R code and parsed using system.file

## Isobaric tagging type. Use one of the following:
# type='iTRAQ4plexSpectra'
# type='iTRAQ8plexSpectra'
# type='TMT2plexSpectra'
# type='TMT6plexSpectra'
type=NULL
isotope.impurities=NULL

## Name of project, by default the name of working directory
## Will be title and author of the analysis reports.
name=basename(getwd())
author="isobar R package"
ibspectra=paste(name,"ibspectra.csv",sep=".")

## When replicates or 'samples belonging together' are analyzed,
## a ProteinGroup object based on all data should be constructed
## beforehand. This then acts as a template and a subset is used.
protein.group.template=NULL

## Where should cached files be saved?
# cachedir="cache"
cachedir="."

## An ibspectra object can be generated from peaklists and identifications.

## peaklist files, by default all mgf file in directory
peaklist=list.files(pattern="*\\.mgf")
## id files, by default all id.csv files in directory
identifications=list.files(pattern="*\\.id.csv")
## mapping files, for data quantified and identified with different but
## correspoding spectra. For example corresponding HCD-CID files.

## masses and intensities which are outside of the 'true' tag mass 
## +/- fragment.precision/2 are discarded
fragment.precision=0.01
## filter mass outliers
fragment.outlier.prob=0.001

readIBSpectra.args = list(
    mapping.file=NULL
)

correct.isotope.impurities=TRUE

normalize=TRUE
normalize.use.protein=NULL
normalize.exclude.protein=NULL
normalize.function=median
normalize.exclude.set = list (seppro_igy14=c(
        "P02763",   #  Alpha1-Acid Glycoprotein
        "P01009-1", #  Alpha1-Antitrypsin
        "P19652",   #  Alpha1-Acid Glycoprotein
        "P01023",   #  Alpha2-Macroglobulin
        "P02768-1", #  Albumin
        "P02647",   #  HDL: Apolipoprotein A1
        "P02652",   #  HDL: Apolipoprotein A1
        "P04114",   #  LDL: Apolipoprotein B
        "P01024",   #  Complent C3
        "P02671-1", #  Fibrinogen
        "P00738",   #  Haptoglobin
        "P01876",   #  IgA 1
        "P01877",   #  IgA 2
        "P01857",   #  IgG 1
        "P01859",   #  IgG 2
        "P01860",   #  IgG 3
        "P01861",   #  IgG 4
        "P01871-1", #  IgM
        "P02787"    #  Transferrin
        ));


use.na=FALSE

## the parameter noise.model can be either a NoiseModel object or a file name
data(noise.model.hcd)
noise.model=noise.model.hcd
## If it is a file name, a noise model is estimated as non one-to-one and saved
## into the file. otherwise, the noise model is loaded from the file
# noise.model="noise.model.rda"

## Certain channels can be defined for creation of a noise model
## e.g. if the first and second channel are technical repeats
## If NULL, all channel combinations are taken into account when creating a 
## noise model.
noise.model.channels=NULL
noise.model.minspectra=50

summarize=FALSE
combn.method="interclass"
## class labels. Must by of type character and of same length as number of channels
## I. e. 4 for iTRAQ 4plex, 6 for TMT 6plex
## Example for iTRAQ 4plex:
# class.labels=as.character(c(1,0,0,0))
# class.labels=c("Treatment","Treatment","Control","Control")
class.labels=NULL
combn=NULL

## Analysis report sections: Significant proteins and protein details
show.significant.proteins=FALSE
show.protein.details=TRUE

ratios.opts = list(
    sign.level.sample=0.01,
    sign.level.rat=0.01)

quant.w.grouppeptides=c("bcrabl","bcrabl,bcrabl_p185,bcrabl_t315i")

min.detect=NULL

datbase="Uniprot"
preselected=c()

ratiodistr=NULL
ratiodistr.summarize=FALSE
ratiodistr.summarize.method="global"

write.qc.report=TRUE
write.report=TRUE
write.xls.report=TRUE

sum.intensities=FALSE
regen=FALSE

scratch=list()
