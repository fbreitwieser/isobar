report.level="protein"
compile=TRUE

## report name
name  <- basename(getwd())

## Samples to be combined:
##  Name subdirectories which contain ibspectras and quant tables
samples <- sub("./","",list.dirs(recursive=FALSE),fixed=TRUE)

## regenerate cache files
regen <- FALSE

## The samples to be combined must have been analyzed using a 
## common protein group template, such that the protein group 
## identifiers are the same
protein.group <- "protein.group.rda"

significance.method='p.value'

## Via database or internet connection, informations on proteins 
## (such as gene names and length) can be gathered. protein.info.f
## defines the function which takes a ProteinGroup object as argument
protein.info.f <- getProteinInfoFromUniprot

## Where should cached files be saved? Will be created if it does not exist
# cachedir="cache"
cachedir="."

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

## XLS report format 'wide' or 'long '.

## 'wide' format outputs ratios in separate columns of the same record
## (i.e. one line per protein)
## 'long' format outputs ratios in separate records (i.e. one line per
## ratio)
xls.report.format="wide"
#xls.report.format="long"

zscore.threshold <- 2.5


spreadsheet.format="xlsx"
