
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

## Via database or internet connection, informations on proteins 
## (such as gene names and length) can be gathered. protein.info.f
## defines the function which takes a ProteinGroup object as argument
protein.info.f <- getProteinInfoFromUniprot

## Where should cached files be saved? Will be created if it does not exist
# cachedir="cache"
cachedir="."

spreadsheet.format="xlsx"
