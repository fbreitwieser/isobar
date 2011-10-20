#!/bin/Rscript

.bnds <- function(x,bnd=NULL,min.x=-bnd,max.x=bnd) {
  max(
      min(x,max.x),
      min.x)
}

print_help_and_exit <- function() {
  cat("create_reports.R       Command line tool to create reports using the\n",
      "                       isobar R/Bioconductor package.\n",
      "\n",
      "Usage:\n",
      "  Rscripts create_reports.R -- 'type=\"iTRAQ4plexSpectra\"'\n",
      "\n",
      "  create_reports.R expects all files to be present in the working directory.\n",
      "  Options can be specified on command line or in a properties.conf file\n",
      "  in the working directory.",
      "\n",
      "Options:\n",
      "  Options defined on the command line overwrite the local properties file\n",
      "  overwrite the global properties file. For reference of possible options\n",
      "  take a look at the file ",
      system.file("report","properties.conf",package="isobar"),".\n");
  quit("no");
}

check_args <- function(args) {
  if (any(grepl("--help",args)) || any(grepl("-h",args))) {
    print_help_and_exit();
  } 
}

check_args(commandArgs(TRUE));

message("Loading package isobar ...")
suppressPackageStartupMessages(library(isobar))

source(system.file("report","report-utils.R",package="isobar"))
initialize.report("properties.conf",args=commandArgs(TRUE),env=.GlobalEnv)

pg <- proteinGroup(ibspectra)
ip <- indistinguishableProteins(pg)

# generate XLS report
if(properties.env$write.xls.report) {
  message("Writing isobar-analysis.xls")
  write.t <- function(...,sep="\t",row.names=FALSE,quote=FALSE)
    write.table(...,sep=sep,row.names=row.names,quote=quote)

  ## add 'group' to protein identification table:
  ibs <- as.data.frame(ibspectra)
        #read.table(properties.env$ibspectra,sep="\t",
        #            stringsAsFactors=F,header=TRUE)
  pg.df <- unique(data.frame(group=protein.tbl$group,
                             protein.g=protein.tbl$protein,stringsAsFactors=FALSE))
  ip.df <- data.frame(accession=names(ip),protein.g=ip,stringsAsFactors=FALSE)
  
  mm.df <- merge(pg.df,ip.df,by="protein.g",all.x=TRUE,all.y=FALSE)
  mm.df <- merge(mm.df,ibs,all.x=TRUE,all.y=FALSE,by="accession")
  mm.df <- mm.df[,c("group",colnames(ibs))]

  ## Analysis Properties:
  nn <- reporterNames(ibspectra)
  ii <- rbind(c(":centeracross:Analysis Properties",rep(":centeracross:",length(nn))),
              "",
              c(":centeracross:Isotope Impurity Correction Matrix",
                rep(":centeracross:",length(nn))),
              cbind(c("",nn),rbind(nn,isotopeImpurities(ibspectra))))

  cl <- attr(protein.tbl,"classLabels")
  if (!is.null(cl)) {
    ii <- rbind(ii,
                "",
                c(":centeracross:Class Labels",":centeracross:",rep("",length(nn)-1)))

    for (i in seq_along(reporterNames(ibspectra))) {
      ii <- rbind(ii,c(nn[i],cl[i],rep("",length(nn)-1)))
    }
  }

  ## write tables to tab seperated files files:
  write.t(xls.protein.tbl,file=paste(properties.env$cachedir,"protein_quant.csv",sep="/"))
  write.t(mm.df,file=paste(properties.env$cachedir,"protein_id.csv",sep="/"))  
  write.t(ii,file=paste(properties.env$cachedir,"properties.csv",sep="/"),col.names=FALSE)
  write.t(ibspectra@log,file=paste(properties.env$cachedir,"logged_operations.csv",sep="/"),col.names=NA,row.names=TRUE)

  ## generate perl command line:
  perl.cl <- paste(system.file("pl","tab2xls.pl",package="isobar")," isobar-analysis.xls",
                   " ':autofilter,freeze_col 3:Identifications=",properties.env$cachedir,"/protein_id.csv'",
                   " ':autofilter,freeze_col 3:Quantifications=",properties.env$cachedir,"/protein_quant.csv'",
                   " 'Analysis Properties=",properties.env$cachedir,"/properties.csv'",
                   " 'Log=",properties.env$cachedir,"/logged_operations.csv'",sep="")

  ## generate Excel report (using Spreadsheet::WriteExcel)
  message(perl.cl)
  system(perl.cl)
}

## generate Latex/Sweave report
if(properties.env$write.qc.report) {
  message("Weaving isobar-qc report")
  Sweave(system.file("report","isobar-qc.Rnw",package="isobar"))
}

if(properties.env$write.report) {
  message("Weaving isobar-analysis report")
  Sweave(system.file("report","isobar-analysis.Rnw",package="isobar"))
}
