#!/usr/bin/Rscript
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(boot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(isobar))

source("meta-properties.R")
load("protein.group.rda")
proteinInfo(protein.group) <- getProteinInfoFromBiomart(protein.group)

if (!exists("cols"))
  cols <- c()
cols <- unique(c("ac","r1","r2","lratio","variance",cols))

message("Calculating dNSAF ...")
if (file.exists("dnsaf.rda")) {
  load("dnsaf.rda")
} else {
  dnsaf <- calculate.dNSAF(protein.group)
  save(dnsaf,file="dnsaf.rda")
}
names(dnsaf) <- reporterProteins(protein.group)

message("Merging tables ...")
merged.table <- get.merged.table(samples,cols=cols)
merged.table <- subset(merged.table,r1==merged.table$r1[1])

tbl.wide <- reshape(merged.table,idvar="ac",timevar=c("r2"),direction="wide",drop="r1")
rownames(tbl.wide) <- tbl.wide$ac
all.names <- do.call(rbind,lapply(tbl.wide[,"ac"],get.names,protein.group=protein.group))
tbl.wide$dNSAF <- dnsaf[as.character(tbl.wide$ac)]

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
