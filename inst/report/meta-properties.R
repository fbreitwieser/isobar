
name  <- basename(getwd())
samples <- sub("./","",list.dirs(recursive=FALSE),fixed=TRUE)
regen <- FALSE

protein.group <- "protein.group.rda"
protein.info.f <- getProteinInfoFromUniprot

cachedir="."
