#!/bin/Rscript

{
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
      system.file("report","properties.R",package="isobar"),".\n");
  quit("no");
}

args <- commandArgs(TRUE)

if ("--help" %in% args || "-h" %in% args) {
  print_help_and_exit();
} 

get.arg <- function(name) {
  if (name %in% args) {
    args <<- setdiff(args,name)
    return(TRUE)
  } else {
    return(FALSE)
  }
}

meta.report <- get.arg("--meta")

if (meta.report) 
  properties.file <- "meta-properties.R"
else 
  properties.file <- "properties.R"

if (length(args) > 0 && file.exists(args[length(args)])) {
  properties.file <- args[length(args)]
  args[length(args)] <- NULL
}

do.compile <- get.arg("--compile")
do.zip <- get.arg("--zip")
protein.report <- get.arg("--protein")
peptide.report <- get.arg("--peptide")

xls.report <- get.arg("--xls")
xlsx.report <- get.arg("--xlsx")
qc.report <- get.arg("--qc")
pdf.report <- get.arg("--pdf")


## TODO: parse further arguments

message("started at ",date())
message("Loading package isobar v",packageDescription("isobar")$Version," ...")
suppressPackageStartupMessages(library(isobar))

if (!exists("properties.env",inherits=FALSE)) {
  properties.env <- load.properties(properties.file,
                                    system.file("report",properties.file,package="isobar"),
                                    args=args)
}

if (xls.report || xlsx.report || qc.report || pdf.report) {
  properties.env$write.xls.report <- (xls.report || xlsx.report)
  if (properties.env$write.xls.report) {
    if (xls.report) properties.env$spreadsheet.format <- 'xls'
    if (xlsx.report) properties.env$spreadsheet.format <- 'xlsx'
  }

  properties.env$write.qc.report <- qc.report
  properties.env$write.report <- pdf.report
  do.compile <- TRUE
}


tryCatch({
  if (meta.report) create.reports.f <- create.meta.reports
  else create.reports.f <- create.reports
  create.reports.f(report.type=ifelse(peptide.report,"peptide","protein"),
                         compile=do.compile,zip=do.zip)},
         error=function(e) {
           save.image(file="isobar.fail.rda")
           stop("create.reports exited with an error - saving session to isobar.fail.rda.\n\n  Message: ",
                 as.character(e))
         }
         )
message("finished at ",date())
}
