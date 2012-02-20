library(gplots)
library(boot)
library(RColorBrewer)
library(isobar)

source("meta-properties.R")
load("protein.group.rda")
proteinInfo(protein.group) <- getProteinInfoFromBiomart(protein.group)

source("~/projects/quant/isobar/inst/report/meta-functions.R")

merged.table <- get.merged.table(samples)
merged.table <- subset(merged.table,r1==merged.table$r1[1])

tbl.wide <- reshape(merged.table,idvar="ac",timevar=c("r2"),direction="wide",drop="r1")
rownames(tbl.wide) <- tbl.wide$ac
all.names <- do.call(rbind,lapply(tbl.wide[,"ac"],get.names))

summarized.table <- write.summarized.table(tbl.wide,all.names,cols=unique(merged.table$r2))

ratio.matrix <- as.matrix(tbl.wide[,grep("lratio",colnames(tbl.wide))])
variance.matrix <- as.matrix(tbl.wide[,grep("var",colnames(tbl.wide))])
sel <- !apply(is.na(ratio.matrix),1,any)
ratio.matrix <- ratio.matrix[sel,]
variance.matrix <- variance.matrix[sel,]

m.median <- apply(ratio.matrix,2,median)
normalized.ratio.matrix <- ratio.matrix-matrix(m.median,nrow=nrow(ratio.matrix),ncol=ncol(ratio.matrix),byrow=T)

plot.heatmaps(ratio.matrix)
plot.pairs()
