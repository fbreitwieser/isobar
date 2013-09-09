#!/usr/bin/Rscript

suppressPackageStartupMessages(library(isobar))
suppressPackageStartupMessages(library(plyr))

all.files <- commandArgs(TRUE)
if (length(all.files) < 2)
  stop("usage: msgfpParser.R result.id.csv file1.msgfp.tsv file2.msgfp.tsv ...")

out.file <- all.files[1]
msgfp.files <- all.files[2:length(all.files)]

message("reading MSGF files:")
message("\t",paste(msgfp.files,collapse="\n\t"))

ib.df <- ldply(msgfp.files,isobar:::.read.msgfp.tsv)
ib.df <- ib.df[ib.df[,'msgf.specevalue'] <= 0.05 & ib.df[,'msgf.rawscore'] > 0,]
writeIBSpectra(ib.df,file=out.file,quote=FALSE)
