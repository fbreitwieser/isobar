#!/usr/bin/Rscript

create.protein.group <- function(id.files,file="protein.group.rda") {
	message("Identification files used for ProteinGroup creation: \n\t",paste(id.files,collapse="\n\t"),"\n")
	suppressPackageStartupMessages(library(isobar))

	# create an IBSpectra object with all identifications
	protein.group <- readProteinGroup(id.files)
	proteinInfo(protein.group) <- getProteinInfoFromBiomart(protein.group)
	print(protein.group)
	# save ProteinGroup object
	save(protein.group,file=file)
	message("\nSaved ProteinGroup object to ",file,".\n")
}

message("## CREATE PROTEINGROUP SCRIPT ##\n")
if (length(commandArgs(TRUE)) == 0) 
  stop("Provide at least one id file as argument!\n")
  
create.protein.group(commandArgs(TRUE))

