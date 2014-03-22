### =========================================================================
### ProteinGroup objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("ProteinGroup",
    contains = "VersionedBiobase",
    representation(
                   spectrumToPeptide = "character",
                   spectrumId = "data.frame",
                   peptideSpecificity = "data.frame",
                   peptideNProtein = "matrix",
                   indistinguishableProteins = "character",
                   proteinGroupTable = "data.frame",
                   overlappingProteins = "matrix",
                   isoformToGeneProduct= "data.frame",
                   proteinInfo = "data.frame",
                   peptideInfo = "data.frame"
    ),
    prototype = prototype(
        new("VersionedBiobase",versions=c(ProteinGroup="1.0.0")),
        peptideSpecificity = data.frame(peptide=c(),specificity=c())
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.ProteinGroup.slots <- function(object) {
  peptides <- unique(object@peptideNProtein[,"peptide"])
  proteins <- unique(object@peptideNProtein[,"protein.g"])
  msg <- NULL

# TODO: check this if
#  if (!setequal(peptides,peptideSpecificity[,'peptide']) ||
#      !setequal(proteins,indistinguishableProteins[,'protein.g']) ||
#      !setequal(peptides,spectrumToPeptide[,'peptide']))
#    msg <- "slots are not of equal length"
  return(msg)
}

.valid.ProteinGroup <- function(object) {
  .valid.ProteinGroup.slots(object)
}
 
setValidity("ProteinGroup", .valid.ProteinGroup)

UNSPECIFIC="unspecific"
GROUPSPECIFIC="group-specific"
REPORTERSPECIFIC="reporter-specific"
#PROTEINSPECIFIC="protein-specific"
SPECIFICITIES=c(UNSPECIFIC,GROUPSPECIFIC,REPORTERSPECIFIC)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor and readProteinGroup.
###

setGeneric("ProteinGroup",function(from,template=NULL,proteinInfo=data.frame())
           standardGeneric("ProteinGroup"))


setMethod("ProteinGroup",signature(from="data.frame",template="ProteinGroup",proteinInfo="ANY"),
          function(from,template,proteinInfo=data.frame()) {
      # TODO: exclude groups from beeing reporters when
      #       there's no reporter-specific peptide ?
      message("Creating ProteinGroup from template ... ",appendLF=FALSE)

      from <- .factor.as.character(from)
      if ('accession' %in% colnames(from)) 
        colnames(from)[colnames(from)=='accession'] <- 'protein'
        
      required.cols <- c("spectrum","peptide","modif","protein")
      required.cols.present <- sapply(required.cols,function(x) x %in% colnames(from))
      if (!all(required.cols.present)) {
        stop("Not all required columns ['spectrum','peptide','modif','accession' or 'protein'] present:\n",
             paste(required.cols[!required.cols.present],collapse=" and ")," missing")
      }
      if (!'start.pos' %in% colnames(from)) {
        warning('start.pos not in supplied data.frame, setting to NA')
        from$start.pos <- NA
      }
      
      ## Substitute Isoleucins with Leucins (indistinguishable by Masspec)
      if (!.PEPTIDE.COLS['REALPEPTIDE'] %in% colnames(from)) 
        from <- .fix.il.peptide(from)

      peptideNProtein <- peptideNProtein(template)[peptideNProtein(template)[,"peptide"] %in% from[,'peptide'],]
      protein.group.table <-
        subset(proteinGroupTable(template),protein.g %in% peptideNProtein[,"protein.g"])
      ipt <- indistinguishableProteins(template)
      indistinguishableProteins <- ipt[ipt %in% peptideNProtein[,"protein.g"]]
      peptideSpecificity <-
        subset(peptideSpecificity(template),peptide %in% from[,"peptide"])
      spectrumId <- unique(from[,setdiff(colnames(from),c("protein","start.pos"))])
      spectrumToPeptide <- spectrumToPeptide(template)[
                     names(spectrumToPeptide(template)) %in% from[,"spectrum"]]
      
      isoforms <-template@isoformToGeneProduct[names(indistinguishableProteins),]
      peptideInfo <- subset(template@peptideInfo,peptide %in% from[,'peptide'] & modif %in% from[,"modif"])
      ## TODO: overlappingProteins are missing

      if (length(proteinInfo) == 0 && length(template@proteinInfo) > 0) {
        if (proteinInfoIsOnSpliceVariants(template@proteinInfo))
          proteinInfo <- subset(template@proteinInfo,accession %in% from$protein)
        else
          proteinInfo <- subset(template@proteinInfo,accession %in% isoforms$proteinac.wo.splicevariant)
      }

      message("done")

      return(
          new("ProteinGroup",
              peptideNProtein = peptideNProtein,
              peptideSpecificity = peptideSpecificity,
              proteinGroupTable = protein.group.table,
              indistinguishableProteins = indistinguishableProteins,
              spectrumId = spectrumId,
              spectrumToPeptide = spectrumToPeptide,
              isoformToGeneProduct = isoforms,
              proteinInfo = proteinInfo,
              peptideInfo = peptideInfo)
      )
    }
)

readProteinGroup <- function(id.file,...,identifications.format=NULL,sep="\t",
                         decode.titles=TRUE,trim.titles=FALSE) {
  pp <- .read.idfile(id.file,identifications.format=identifications.format,
                     sep=sep,decode.titles=decode.titles,trim.titles=trim.titles)
  return(ProteinGroup(from=pp,...))
}

readProteinGroup2 <- function(id.file,...,identifications.format=NULL,
                              sep="\t",decode.titles=TRUE,trim.titles=FALSE) {
  identifications <- .read.idfile(id.file,sep=sep,
                                  identifications.format=identifications.format,
                                  decode.titles=decode.titles,trim.titles=trim.titles)

  ## Check that obligatory columns are present
  identifications <- .check.columns(identifications)

  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  ## Substitute Isoleucins with Leucins (indistinguishable by Masspec)
  if (!.PEPTIDE.COLS['REALPEPTIDE'] %in% colnames(identifications)) 
     identifications <- .fix.il.peptide(identifications)

  ## Separate protein columns (focus on peptide-spectrum matches)
  PC <- setdiff(.PEPTIDE.COLS, SC['PEPTIDE'])
  protein.colnames <- which(colnames(identifications) %in% c(SC['PEPTIDE'],PC))
  pept.n.prot <- unique(identifications[,protein.colnames])
  identifications <- unique(identifications[,-which(colnames(identifications) %in% PC)])

  ## Merge identifications
  if (max(table(identifications[,SC['SPECTRUM']])) > 1) {
    if (SC['DISSOCMETHOD'] %in% colnames(identifications) && 
        length(unique(identifications[,SC['DISSOCMETHOD'] ]))>1) {
      identifications <- ddply(identifications,'dissoc.method',.merge.identifications)
      identifications <- .merge.quant.identifications(identifications)
    } else {
      identifications <- .merge.identifications(identifications)
    }
  }
  identifications <- .remove.duplications(identifications)
 
  # Create ProteinGroup
  proteinGroup <- ProteinGroup(merge(pept.n.prot,identifications,by="peptide"),...)
  
  return(proteinGroup)
}

setMethod("ProteinGroup",signature(from="data.frame",template="NULL",proteinInfo="ANY"),
          function(from,template,proteinInfo) ProteinGroup(from,proteinInfo=proteinInfo))

setMethod("ProteinGroup",signature(from="data.frame",template="missing",proteinInfo="ANY"),
          function(from,proteinInfo) {
      message("Creating ProteinGroup ... ",appendLF=FALSE)

      from <- .factor.as.character(from)
      if ('accession' %in% colnames(from)) 
        colnames(from)[colnames(from)=='accession'] <- 'protein'
        
      required.cols <- c("spectrum","peptide","modif","protein")
      required.cols.present <- sapply(required.cols,function(x) x %in% colnames(from))
      if (!all(required.cols.present)) {
        stop("Not all required columns ['spectrum','peptide','modif','accession' or 'protein'] present:\n",
             paste(required.cols[!required.cols.present],collapse=" and ")," missing")
      }
      if (!'start.pos' %in% colnames(from)) {
        warning('start.pos not in supplied data.frame, setting to NA')
        from$start.pos <- NA
      }
      
      ## Substitute Isoleucins with Leucins (indistinguishable by Masspec)

     if (!.PEPTIDE.COLS['REALPEPTIDE'] %in% colnames(from)) 
       from <- .fix.il.peptide(from)

      subset.s <- function(my.df,j) {
        if (is(my.df,"data.table"))
          my.df[,j,with=FALSE]
        else
          my.df[,j]
      }

      spectrumToPeptide <- .as.vect(unique(subset.s(from,c("spectrum","peptide"))))
      ##spectrumToPeptide <- setNames(from[["peptides"]],from[["spectrum"]])

      spectrumId <- unique(subset.s(from,setdiff(colnames(from),c("protein","start.pos","aa.before","aa.after"))))
      peptideInfo <- unique(subset.s(from,intersect(c("protein","peptide","start.pos","aa.before","aa.after","modif","real.peptide"),colnames(from))))
      peptideInfo <- peptideInfo[order(peptideInfo[["protein"]],
                                       peptideInfo[["start.pos"]],
                                       peptideInfo[["peptide"]]),]
      from <- unique(subset.s(from,c("protein","peptide")))
    
      # use numbers to not exceed the maximum length
      from$peptn <- as.integer(as.factor(from[["peptide"]]))
      from$protn <- as.integer(as.factor(from[["protein"]]))
      
      # with which peptide-combination are proteins detected?
      combn.peptides <- tapply(from[["peptn"]],from[["protein"]],
                               function(s) paste0(sort(unique(s)),collapse="-"))
      
      # group proteins with same peptide combinations together
      # - those are indistiguishable
      combn.proteins <- tapply(names(combn.peptides),combn.peptides,
                               function(x) paste0(x,collapse=","))
      
      # merge data frames
      f.peptides <- combn.peptides[ from[["protein"]] ]
      from[["protein.g"]] <- combn.proteins[ f.peptides ]
                                                
      pep.n.prots <- unique(subset(from,,c("peptide","protein.g")))
      all.protein.g <- table(pep.n.prots[['protein.g']])
      prots.to.prot <- as.matrix(unique(subset.s(from,c("protein.g","protein"))))
      prots.group <- c()
      
      # GROUPING:
      # 1) determine proteins w/ specific peptides [reporterProteins]
      ssp <- table(pep.n.prots[["peptide"]])==1
      sel.ssp <- pep.n.prots[["peptide"]] %in% names(ssp)[ssp]
      # take first those proteins which explain most specific peptides
      tt <- table(pep.n.prots[["protein.g"]][sel.ssp])
      protein.reporterProteins <- names(sort(tt,decreasing=TRUE))
      protein.groupMembers <- setdiff(names(all.protein.g),protein.reporterProteins)
      # proteins to consider for grouping
      pg.length <- length(all.protein.g)

      pgt <- new.env()
      pgt[["protein.group.table"]] <- data.frame(reporter.protein=rep("",pg.length),
                                        protein.g=rep("",pg.length),
                                        'n.reporter-specific'=0,
                                        'n.group-specific'=0,
                                        'n.unspecific'=0,
                                        stringsAsFactors=FALSE,check.names=FALSE)
      pgt[["my.idx"]] <- 1
      pgt[["pep.n.prots"]] <- pep.n.prots
      pgt[["prots.to.consider"]] <- new.env()
      for (my.protein.g in protein.groupMembers) {
        pgt[["prots.to.consider"]][[my.protein.g]] <- 
          subset(pep.n.prots,protein.g == my.protein.g)[["peptide"]]
      }

      # 2) get proteins, which are contained by those detected proteins
      #    (i.e. no specific and no overlapping peptides)
      for (reporter.protein in protein.reporterProteins) {
        .build.protein.group(pgt,reporter.protein)
      }

      # proteins left in pep.n.prots.toconsider are those proteins
      # with overlapping peptides only. 
      # 3. They might have no reporter-specific peptides, but group-specific peptides

      grouped.proteins <- pgt[["protein.group.table"]][["protein.g"]]
      grouped.peptides <- subset(pep.n.prots,protein.g %in% protein.reporterProteins)[["peptide"]]
      avail.peptides.n.prots <- subset(pep.n.prots,!peptide %in% grouped.peptides)
      ungrouped.peptides <- subset(pep.n.prots,!peptide %in% grouped.peptides)[["peptide"]]

      while (length(ls(env=pgt[["prots.to.consider"]])) > 0) {
        sorted.proteins <- sort(unlist(eapply(pgt[["prots.to.consider"]],length)))
        my.protein.g <- names(sorted.proteins)[length(sorted.proteins)]
       
        remove(list=my.protein.g,envir=pgt[["prots.to.consider"]])
        .build.protein.group(pgt,my.protein.g)
      }

      pgt[["protein.group.table"]]$is.reporter <-
        pgt[["protein.group.table"]]$protein.g == pgt[["protein.group.table"]]$reporter.protein
      
      # define peptide specificity - pep.n.prots.toconsider are ignored here for now
      tmp <- merge(as.data.frame(pgt[["protein.group.table"]],stringsAsFactors=FALSE),
                   pep.n.prots,by.x="protein.g",by.y="protein.g")
      
      peptideSpecificity <- .factor.as.character(ddply(tmp,"peptide",function(d) {
        data.frame(specificity=
                   ifelse(length(unique(d[,"protein.g"]))==1,REPORTERSPECIFIC,
                          ifelse(length(unique(d[,"reporter.protein"]))==1,GROUPSPECIFIC,
                                 UNSPECIFIC)),
                   n.shared.proteins=length(unique(d[,"protein.g"])),
                   n.shared.groups=length(unique(d[,"reporter.protein"])),
                   stringsAsFactors=FALSE)
      }
      ))

      # isoforms (handled for Uniprot only, ATM)
      proteins <- as.character(sort(unique(prots.to.prot[,"protein"])))
      isoforms <- data.frame(database = .get.database(proteins),
                             proteinac.w.splicevariant = proteins,
                             proteinac.wo.splicevariant = proteins,
                             splicevariant = NA,stringsAsFactors=FALSE)
      
      sel.uniprot <- isoforms$database == .UNIPROTDATABASE
      pos.isoforms <- sel.uniprot & grepl("^[^-]*\\-[0-9]*$",proteins)
      isoforms$proteinac.wo.splicevariant[pos.isoforms] <- 
         sub("-[^-]*$","",proteins[pos.isoforms])
      isoforms$splicevariant[pos.isoforms] <- sub(".*-","",proteins[pos.isoforms])
      rownames(isoforms) <- proteins

      message("done")
      
      return(
          new("ProteinGroup",
              peptideNProtein = as.matrix(pep.n.prots),
              peptideSpecificity = peptideSpecificity,
              proteinGroupTable = pgt[["protein.group.table"]],
              indistinguishableProteins =
                .as.vect(prots.to.prot,col.data='protein.g',col.names='protein'),
              spectrumId = spectrumId,
              spectrumToPeptide = spectrumToPeptide,
              isoformToGeneProduct = isoforms,
              peptideInfo = peptideInfo,
              proteinInfo = proteinInfo)
      )
      
    }
)


.uniprot.pattern.ac <- "^[A-Z][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]-?[0-9]*$"
.uniprot.pattern.id <- "^[A-Z0-9]{2,5}_[A-Z9][A-Z0-9]{2,5}$"
.entrez.pattern.id <- "^gi\\|[0-9]{5,15}$"
.numeric.pattern.ac <- "^[0-9]{5,15}$"

.UNIPROTDATABASE <- "UniprotKB"
.UNIPROTDATABASE_ID <- "UniprotKB IDs" # Maybe: add support to get Protein Info for Uniprot IDs, as some people use them (but they are not stable, see http://www.uniprot.org/faq/6)
.ENTREZDATABASE <- "Entrez Protein"

.get.database <- function(acs) {
  res <- rep("unknown",length(acs))
  sel.uniprot <- grepl(.uniprot.pattern.ac,acs)
  sel.entrez <- grepl(.entrez.pattern.id,acs)
  sel.numeric <- grepl(.numeric.pattern.ac,acs)

  if (any(sel.uniprot))
    res[sel.uniprot] <- .UNIPROTDATABASE

  if (any(sel.entrez))
    res[sel.entrez] <- .ENTREZDATABASE 

  if (any(sel.numeric))
    res[sel.numeric] <- "numeric"

  return(res)  
}




.build.protein.group <- function(pgt,reporter.protein) {

    pep.for.reporter <- subset(pgt[["pep.n.prots"]],protein.g==reporter.protein)[["peptide"]]

    # define proteins which have a subset of the reporter protein peptides
    groupmember.pept <- eapply(pgt[["prots.to.consider"]],function(pep.for.prot) {
      if (length(pep.for.prot) < length(pep.for.reporter) & 
              all(pep.for.prot %in% pep.for.reporter))
        pep.for.prot
      else
        NA
    })
   
    groupmember.pept <- groupmember.pept[!is.na(groupmember.pept)]
    groupmembers <- names(groupmember.pept)
    remove(list=groupmembers,envir=pgt[["prots.to.consider"]])

    my.idxes <- seq(from=pgt[["my.idx"]],to=pgt[["my.idx"]]+length(groupmembers))
    pgt[["protein.group.table"]][my.idxes,'reporter.protein'] <- reporter.protein
    pgt[["protein.group.table"]][my.idxes,'protein.g'] <- c(reporter.protein,groupmembers)
    # define number of reporter-/group-/non-specific peptides for group members

    all.peptides <- pep.for.reporter
    unspecific.pept <- subset(pgt[["pep.n.prots"]],peptide %in% all.peptides &
                              !protein.g %in% c(reporter.protein,groupmembers))[["peptide"]]
    groupspecific.pept <- all.peptides[all.peptides %in% as.character(unlist(groupmember.pept)) & !all.peptides %in% unspecific.pept]
    reporterspecific.pept <- setdiff(pep.for.reporter,c(groupspecific.pept,unspecific.pept))

    pgt[["protein.group.table"]][my.idxes,'n.reporter-specific'] <- 
      c(length(reporterspecific.pept),rep(0,length(groupmembers)))
    pgt[["protein.group.table"]][my.idxes,'n.group-specific'] <- 
      c(length(groupspecific.pept),
        sapply(groupmember.pept,function(x) sum(x %in% groupspecific.pept)))
    pgt[["protein.group.table"]][my.idxes,'n.unspecific'] <- 
      c(length(unspecific.pept),
        sapply(groupmember.pept,function(x) sum(x %in% unspecific.pept)))

    pgt[["my.idx"]] = pgt[["my.idx"]]+length(groupmembers)+1
}

getProteinInfoFromBiomart <- function(x,database="Uniprot") {
  protein.acs <- x@isoformToGeneProduct[names(indistinguishableProteins(x)),"proteinac.wo.splicevariant"]

  protein.info <- data.frame(accession=c(),name=c(),protein_name=c(),
                             gene_name=c(),organism=c())
  if (database == "Uniprot") {
    #require(biomaRt)
    tryCatch({
      mart <- biomaRt::useMart("uniprot_mart",dataset="UNIPROT",
                      host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")
    protein.info <- biomaRt::getBM(attributes=c("accession","name","protein_name","gene_name","organism"),
          filter='accession',values=protein.acs,mart=mart)
    },error=function(e) warning("could not set Biomart"))

  } else {
    stop(sprintf("getProteinInfo for database %s not defined.",database))
  }
  protein.info$gene_name[!is.na(protein.info$gene_name) & protein.info$gene_name == ""] <- NA

  return(protein.info)
}

# for NCBI protein: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=115298678

getProteinInfoFromTheInternet <- function(x) {
  if (is.character(x)) {
    protein.acs <- x
  } else {
    my.df <- x@isoformToGeneProduct
    protein.acs <- unique(my.df[,"proteinac.wo.splicevariant"])
    if ("database" %in% colnames(my.df)) {
      sel.uniprot <- my.df[["database"]] == .UNIPROTDATABASE
      sel.entrez <- my.df[["database"]] == .ENTREZDATABASE
    }
  }

  if (!exists("sel.uniprot")) {
    sel.uniprot <- grepl(.uniprot.pattern.ac,protein.acs)
    sel.entrez <- grepl(.entrez.pattern.id,protein.acs)
  }

  res <- data.frame()
  if (any(sel.uniprot)) 
    res <- getProteinInfoFromUniprot(protein.acs[sel.uniprot])

  if (any(sel.entrez))
    res <- rbind.fill(res,getProteinInfoFromEntrez(protein.acs[sel.entrez]))

  return(res)
}

getProteinInfoFromUniprot <- 
  function(x,splice.by=200, 
           fields = c(accession="id",name="entry%20name",protein_name="protein%20names",
                      gene_name="genes",organism="organism", length="length",
                      sequence="sequence")) {
 
  if (is.character(x))
    protein.acs <- x
  else  
    protein.acs <- x@isoformToGeneProduct[names(indistinguishableProteins(x)),"proteinac.wo.splicevariant"]
  protein.acs <- gsub("%","",protein.acs)
  protein.info <- c()
  i <- 1
  while (i <= length(protein.acs)) {
    cat(".")
    uniprot.url <- paste0("http://www.uniprot.org/uniprot/?query=",
                 paste0("accession:",protein.acs[seq(from=i,to=min(length(protein.acs),i+splice.by-1))],collapse="+OR+"),
                 "&format=tab&compress=no&columns=",
                 paste0(fields,collapse=","))
    if (isTRUE(getOption('isobar.verbose')))
      message("fetching protein info from ",uniprot.url)
    protein.info <- rbind(protein.info,read.delim(url(uniprot.url),stringsAsFactors=FALSE,col.names=names(fields)))
    i <- i + splice.by
  }
  if (!is.null(protein.info) && nrow(protein.info) > 0) {
    protein.info$protein_name <- sapply(strsplit(protein.info$protein_name," (",fixed=TRUE),function(x) x[1])
    protein.info$gene_name <- sapply(strsplit(protein.info$gene_name," "),function(x) x[1])
    protein.info$sequence <- gsub(" ","",protein.info$sequence)
  } else {
    warning("getProteinInfoFromUniprot returned no results for ",protein.acs," accessions")
  }
  return(protein.info)
}

getProteinInfoFromEntrez <- function(x,splice.by=200) {
  require(XML)

  eutils.url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id="
  if (is.character(x))
    protein.acs <- x
  else  
    protein.acs <- x@isoformToGeneProduct[names(indistinguishableProteins(x)),"proteinac.wo.splicevariant"]

  if (all(grepl("^gi|",protein.acs))) {
    new.protein.acs <- sapply(strsplit(acs[sel.entrez],"|",fixed=TRUE),function(x) x[2]) # strip ^gi| 
    protein.conv <- setNames(protein.acs,new.protein.acs)
    protein.acs <- new.protein.acs
  }

  protein.info <- c()
  i <- 1
  while (i < length(protein.acs)) {
    cat(".")

    full.url <- paste0(eutils.url,paste(as.character(protein.acs[seq(from=i,to=min(length(protein.acs),i+splice.by-1))]),collapse=","))
    res <- xmlToList(full.url)
    res.l <- lapply(res, function(x) {
      unlist(setNames(lapply(x[names(x) != 'Id'], function(y) 
        if ('.attrs' %in% names(y) && y$.attrs['Name'] != 'Comment') 
          setNames(y$text,y$.attrs['Name'])),NULL))
    })

    res.l <- res.l[!sapply(res.l,is.null)]
  
    enames.to.isobar <- c(Gi="gi_accession",Caption="name",Title="protein_name")
    res.df <- ldply(res.l, function(x) { 
                    setNames(x[names(enames.to.isobar)],enames.to.isobar)
    })
  
   protein.info <- rbind(protein.info,res.df)
    i <- i + splice.by
  }
  protein.info[,c('protein_name','organism')] <- t(sapply(strsplit(sub("^(.*) \\[(.*)\\]$","\\1\t\\2",protein.info[,'protein_name']),"\t"),function(x) x))
 
  protein.info$gene_name <- NA
  protein.info$.id <- NULL
  protein.info$accession <- protein.conv[protein.info$gi_accession]
  protein.info
}

getProteinInfoFromNextProt <- function(x) {
  fields <- c(accession="id",name="entry%20name",protein_name="protein%20names",
              gene_name="genes",organism="organism",
              length="length",sequence="sequence")

  stop("NOT IMPLEMENTED YET")
   
  protein.acs <- x@isoformToGeneProduct[names(indistinguishableProteins(x)),]
  protein.info  <- ddply(protein.acs,"proteinac.wo.splicevariant",
                         function(y) {
                         url <- sprintf("http://www.nextprot.org/rest/protein/NX_%s?format=json",
                                        unique(y$proteinac.wo.splicevariant))
                         # PARSE JSON, esp isoform sequences
  })
  attr(protein.info,"on.splice.variant") <- TRUE
  return(protein.info)
}

proteinInfoIsOnSpliceVariants <- function(protein.info) {
  return(isTRUE(attr(protein.info,"on.splice.variant")))
}


getProteinInfoFromBioDb <- function(x,...,con=NULL) {
  if (is.null(con)) {
    con <- dbConnect(...)
    do.disconnect <- TRUE
  } else {
    do.disconnect <- FALSE
  }

  if (is.character(x))
    protein.acs <- x
  else  
    protein.acs <- names(indistinguishableProteins(x))
  query <- paste("SELECT primaryac AS accession,id AS name,",
                 "  description AS protein_name,",
                 "  (SELECT string_agg(g.genename,',') FROM genenames g WHERE g.entryid=d.entryid AND g.synonym=FALSE AND g.sourcedb=3) AS gene_name,",
                 "  os AS organism, seqlength as length, sequence",
                 "FROM dbentries d ",
                 "WHERE dbid IN (2,3) AND primaryac IN (",paste0("'",protein.acs,"'",collapse=","),")")
  res <- dbGetQuery(con,query)
  if (do.disconnect)
    dbDisconnect(con)
  attr(res,"on.splice.variant") <- TRUE
  return(res)
}

getPtmInfoFromPhosphoSitePlus <- function(protein.group,file.name=NULL,modif="PHOS",
                                          psp.url="http://www.phosphosite.org/downloads/",
                                          mapping=c(PHOS="Phosphorylation_site_dataset.gz",
                                                    ACET="Acetylation_site_dataset.gz",
                                                    METH="Methylation_site_dataset.gz",
                                                    SUMO="Sumoylation_site_dataset.gz",
                                                    UBI="Ubiquitination_site_dataset.gz")) {
  if (length(modif) > 1) {
    return(do.call(rbind,lapply(modif,function(m) getPtmInfoFromPhosphoSitePlus(protein.group,file.name,m,psp.url,mapping))))
  }

  if (is.null(file.name)) file.name <- mapping[modif]

  if (!file.exists(file.name) && is.null(modif)) stop("provide PhosphoSitePlus file or modif name")
  if (!file.exists(file.name) && !is.null(modif)) {
    download.file(paste0(psp.url,mapping[modif]),mapping[modif])
    file.name <- mapping[modif]
  }

  sites <- read.delim(file.name,
                      sep="\t",header=TRUE,skip=3,stringsAsFactors=FALSE)
  ac.column <- ifelse("ACC_ID" %in% colnames(sites),"ACC_ID","ACC.")
  species.column <- ifelse("ORG" %in% colnames(sites),"ORG","SPECIES")
  residue.column <- ifelse("MOD_RSD" %in% colnames(sites),"MOD_RSD","RSD")

  sites <- sites[gsub("-.*","",sites[,ac.column]) %in% gsub("-.*","",names(indistinguishableProteins(protein.group))),]

  sites$PUBMED_LTP[!is.na(sites$PUBMED_LTP)] <- paste("n.publ ltp:",sites$PUBMED_LTP[!is.na(sites$PUBMED_LTP)])
  sites$PUBMED_MS2[!is.na(sites$PUBMED_MS2)] <- paste("n.publ htp:",sites$PUBMED_MS2[!is.na(sites$PUBMED_MS2)])

  data.frame(.id=sites[,ac.column],
             isoform_ac=sapply(sites[,ac.column],function(ac) ifelse(grepl("-[0-9]$",ac),ac,paste0(ac,"-1"))),
             species=tolower(sites[,species.column]),
             modification.name=tolower(sites[,'MOD_TYPE']),
             description=apply(sites,1,function(x) {
                               y <- tolower(x['MOD_TYPE'])
                               if (nchar(x['IN_DOMAIN']) > 0)
                                 y <- paste0(y," (domain ",x['IN_DOMAIN'],")")
                               y
                             }),
             evidence=apply(sites[,c("PUBMED_LTP","PUBMED_MS2")],1,
                            function(x) { x<-x[!is.na(x)]; paste(x,collapse=";")}),
             position=as.integer(substr(sites[,residue.column],2,nchar(sites[,residue.column]))),
             stringsAsFactors=FALSE)
}

getPtmInfoFromNextprot <- function(protein.group,
                                   nextprot.url="http://www.nextprot.org/rest/entry/NX_XXX/ptm?format=json",
                                   url.wildcard="XXX") {
  protein.acs <- unique(protein.group@isoformToGeneProduct$proteinac.wo.splicevariant)
  require(RJSONIO)
  pb <- txtProgressBar(max=length(protein.acs),style=3)

  nextprot.ptmInfo <- 
    lapply(seq_along(protein.acs),function(ac_i) {
           setTxtProgressBar(pb,ac_i)
           tryCatch(fromJSON(sub(url.wildcard,protein.acs[ac_i],nextprot.url)),
                    error=function(e) {
              if (isTRUE(getOption('isobar.verbose')))
                warning("Could not fetch from ",
                        sub(url.wildcard,protein.acs[ac_i],nextprot.url),":",
                        e$message)
              return()})
  })
  names(nextprot.ptmInfo) <- protein.acs
  sel.length0 <- sapply(nextprot.ptmInfo,length)==0
  if (any(sel.length0))
    warning("Could not fetch neXtProt results for some proteins:",
            " ",.abbrev(protein.acs[sel.length0],4,", "),".")
  nextprot.ptmInfo <- nextprot.ptmInfo[sapply(nextprot.ptmInfo,length)>0]
  ptm.info <- ldply(nextprot.ptmInfo,
                    function(x) 
                      ldply(x,function(y) {
                            y[sapply(y,is.null)] <- NA
                            y$modification.name <- y$modification['name']
                            y$modification.accession <- y$modification['accession']
                            y$modification <- NULL
                            data.frame(y,stringsAsFactors=FALSE)
                      })
                    )
  ptm.info$isoform_ac <- sub("^NX_","",ptm.info$isoform_ac)
  ptm.info$position <- ptm.info$first_position
  ptm.info
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

as.data.frame.ProteinGroup <- function(x,...) as(x,"data.frame")

setAs("data.frame","ProteinGroup",function(from) ProteinGroup(from))
setMethod("as.data.frame",signature(x="ProteinGroup"),
          function(x, row.names=NULL, optional=FALSE, ...) as(x,"data.frame"))

setAs("ProteinGroup","data.frame",function(from) {
  sp.df <- from@spectrumId
  ip.df <- .vector.as.data.frame(indistinguishableProteins(from),
                                 colnames=c("protein","protein.g"))
  
  ## merge spectrum to peptide to indistinguishable proteins (protein.g)
  spectrum.to.proteing <-
    merge(sp.df, as.data.frame(peptideNProtein(from),stringsAsFactors=FALSE))
 
  ## merge with indist. proteins
  spectrum.to.protein <- merge(ip.df,spectrum.to.proteing)
  spectrum.to.protein.w.peptideinfo <- merge(spectrum.to.protein,from@peptideInfo)

  return(spectrum.to.protein.w.peptideinfo[,c("spectrum","peptide","modif","start.pos","protein")])
})

## Export a /concise/ data.frame which is peptide centric:
##  - for each peptide there is one line only
##  - column proteins concatenates all the different protein/protein group ACs
##    - indistinguishable peptides are separated by ',', protein groups by ';'
##  - column n.groups: number of groups a peptide appears
##  - column n.acs: number of acs a peptide appears in (without splice variant!)
##      counts double if the same ACs are in different groups 
##  - column n.variants: number of acs a peptide appears in (with splice variant!)
proteinGroup.as.concise.data.frame <-
      function(from) {
        ip.df <- .vector.as.data.frame(indistinguishableProteins(from),
                                       colnames=c("protein","protein.g"))
        ip.df <- merge(ip.df,from@isoformToGeneProduct,
                       by.x="protein",by.y="proteinac.w.splicevariant")
        pep.n.prot <- as.data.frame(peptideNProtein(from),stringsAsFactors=FALSE)
        in.df <- unique(ddply(ip.df, "protein.g",
                        function(x) 
                          c(n.acs=length(unique(x[,"proteinac.wo.splicevariant"])),
                            n.variants=length(unique(x[,"protein"])))))
        pep.n.prot <- merge(pep.n.prot,ip.df)
        pep.n.prot <- merge(pep.n.prot,proteinGroupTable(from)[,c("protein.g","reporter.protein")])
        res <- ddply(pep.n.prot,"peptide",
                     function(x) {
                       res <- data.frame(n.acs=length(unique(x[,"proteinac.wo.splicevariant"])),
                                         n.variants=length(unique(x[,"protein"])))
                       x <- unique(x[,c("reporter.protein","protein.g")])
                       protein.gs <- unique(x[,'reporter.protein'])
                       res <- cbind(proteins=paste(tapply(x$protein.g,factor(x$reporter.protein),
                                                          paste,collapse=","),collapse=";"),
                                    n.groups=length(protein.gs),
                                    res,stringsAsFactors=FALSE)
                       if (!is.null(attr(from,"protein.group.ids"))) 
                         res  <- cbind(groups=paste(attr(from,"protein.group.ids")[protein.gs],collapse=","),
                                       res,stringsAsFactors=FALSE)
                       res

                     })
        return(unique(res))
      }


.proteinGroupAsConciseDataFrame1 <- 
  function(from,only.reporters=TRUE,show.proteinInfo=TRUE,
           human.protein.acs=TRUE,show.startpos=TRUE,modif.pos=NULL,
           ptm.info=NULL,link.url="http://www.uniprot.org/uniprot/") {

        pep.n.prot <- merge(as.data.frame(peptideNProtein(from),stringsAsFactors=FALSE),
                            from@peptideInfo)
        p.ac <- .protein.acc(reporterProteins(from),from)
        ip.df <- .vector.as.data.frame(indistinguishableProteins(from),
                                       colnames=c("protein","protein.g"))
        if (only.reporters)
          ip.df <- ip.df[ip.df$protein.g %in% reporterProteins(from),]
        ip.df <- merge(ip.df,from@isoformToGeneProduct,
                       by.x="protein",by.y="proteinac.w.splicevariant")
        ip.df <- merge(ip.df,proteinGroupTable(from)[,c("protein.g","reporter.protein")])
        ipp.df <- merge(ip.df,pep.n.prot)

        in.df <- ddply(ipp.df, c("reporter.protein"),
          function(x) {
            # for each 'protein AC' (no splice variant), summarize the splice variants (eg P123-[1-3,5])
            merged.splicevariants <- ddply(x,"proteinac.wo.splicevariant",.summarize.splice.variants,link.url=link.url)

            # for each modified peptide, get one record summarizing its information
            merged.pepmodifs <- ddply(x,c("peptide","modif"),.summarize.pepmodif,modif.pos=modif.pos,ptm.info=ptm.info,from=from)
            
            data.frame(proteinn=paste(merged.splicevariants$ac,collapse=","),
                       link=merged.splicevariants$link[1],
                       merged.pepmodifs,stringsAsFactors=FALSE)
        })
        merge(ip.df,in.df)
  }



# nexprot url: "http://www.nextprot.org/db/entry/NX_"
.proteinGroupAsConciseDataFrame <- 
  function(from,only.reporters=TRUE,show.proteinInfo=TRUE,
           human.protein.acs=TRUE,show.startpos=TRUE,modif.pos=NULL,
           ptm.info=NULL,link.url="http://www.uniprot.org/uniprot/",
           report.only.modified=FALSE) {
        i = 0

        pep.n.prot <- merge(as.data.frame(peptideNProtein(from),stringsAsFactors=FALSE),
                            from@peptideInfo,by="peptide")
        if (!is.null(modif.pos) && isTRUE(report.only.modified))
          pep.n.prot <- pep.n.prot[grepl(modif.pos,pep.n.prot[,'modif']),]

        p.ac <- .protein.acc(reporterProteins(from),from)

        # ip.df will contain information on
        #   protein.g,protein,proteinac.wo.splicevariant,splicevariant,reporter.protein

        ip.df <- .vector.as.data.frame(indistinguishableProteins(from),
                                       colnames=c("protein","protein.g"))
        if (only.reporters)
          ip.df <- ip.df[ip.df[,'protein.g'] %in% reporterProteins(from),]
        ip.df <- merge(ip.df,from@isoformToGeneProduct,
                       by.x="protein",by.y="proteinac.w.splicevariant")
        ip.df <- merge(ip.df,proteinGroupTable(from)[,c("protein.g","reporter.protein")],by="protein.g")
        ipp.df <- merge(ip.df,pep.n.prot)

        in.df <- ddply(ipp.df, c("reporter.protein"),
          function(x) {
            ## for each 'protein AC' (no splice variant), 
            ##  summarize the splice variants (eg P123-[1-3,5])
            merged.splicevariants <- 
              ddply(x,"proteinac.wo.splicevariant",
                    .summarize.splice.variants,link.url=link.url)
            # for each modified peptide, get one record summarizing its information
            merged.pepmodifs <- ddply(x,c("peptide","modif"),.summarize.pepmodif,
                                      modif.pos=modif.pos,ptm.info=ptm.info,from=from)
            data.frame(proteinn=paste(merged.splicevariants[,'ac'],collapse=","),
                       link=merged.splicevariants[1,'link'],
                       merged.pepmodifs,stringsAsFactors=FALSE)
        },.parallel=isTRUE(getOption('isobar.parallel')))

        cols <- c("proteinn","link","start.pos","real.peptide","modif","modif.pos","modif.comment")
        if (show.proteinInfo) 
          pnd <- proteinNameAndDescription(from,ip.df[,'protein.g'])


        my.res <- 
          ddply(merge(ip.df,in.df,by="reporter.protein"),c("peptide","modif"),
                function(x) {
                  #protein.gs <- unique(x[,'reporter.protein'])
                  protein.gs <- unique(x[,'protein.g'])
                  x <- unique(x[,cols[cols %in% colnames(x)]])

                  res <- data.frame(start.pos=.unique.or.collapse(x[,'start.pos'],";"),
                                    proteins=paste(unique(x[,'proteinn']),collapse=";"),
                                    uniprot=paste0("@link=",x[1,'link'],"@x",collapse=";"),
                                    stringsAsFactors=FALSE)

                  ## n.acs=length(unique(x[,"proteinac.wo.splicevariant"])),
                  ## n.variants=length(unique(x[,"protein"])),
                  ## n.variants=length(protein.gs),

                  if ('real.peptide' %in% colnames(x))
                    res <- 
                      cbind(res,
                            real.peptide=.unique.or.collapse(x[,'real.peptide'],";"),
                            stringsAsFactors=FALSE)

                  if (!is.null(modif.pos)) {
                    null.comments <- x[,'modif.comment'] == ""
                    res <- cbind(res,
                                 modif.pos=ifelse(!any(x[,'modif.pos']!=0),"",
                                                  paste(x[,'modif.pos'],collapse=";")),
                                 comment=ifelse(all(null.comments),"",
                                                paste(x[,'modif.comment'][!null.comments],
                                                collapse="\n")),
                                 stringsAsFactors=FALSE)
                  }

                  if (show.proteinInfo && nrow(pnd) > 0) 
                    res <- cbind(res,
                                 lapply(pnd[ip.df[,'protein.g'] %in% protein.gs,],
                                        .paste_unique,collapse=";"))

                  if (!is.null(attr(from,"from.ids"))) 
                    res  <- cbind(groups=paste(attr(from,"from.ids")[protein.gs],
                                               collapse=","),
                                  res,stringsAsFactors=FALSE)

                  res

                },.parallel=isTRUE(getOption('isobar.parallel')))

        # res[,'peptide'] <- .convertModifToPos(res[,'peptide'],res[,'modif'])
        return(unique(my.res))
}

.summarize.pepmodif <- function(x,modif.pos,ptm.info,from) {
  res <- data.frame(peptide=x$peptide[1],start.pos=paste(x$start.pos,collapse=";"),
                    stringsAsFactors=FALSE)
  if ('real.peptide' %in% colnames(x)) 
    res <- cbind(res,real.peptide=.unique.or.collapse(x[,'real.peptide']),stringsAsFactors=FALSE)
  
  if (!is.null(modif.pos)) 
    res <- cbind(res,.get.modif.pos(x,modif.pos,ptm.info,from))
  res
}

.summarize.splice.variants <- function(x,link.url) {
  res <- c(ac=unique(x$proteinac.wo.splicevariant),
           link=paste0(link.url,unique(x$proteinac.wo.splicevariant)))
  if (!all(is.na(x$splicevariant))) {
    if (length(unique(x$splicevariant))==1) { res['ac'] <- x$protein[1] 
    } else { res['ac'] <- sprintf("%s-[%s]",unique(x$proteinac.wo.splicevariant),
                              number.ranges(as.integer(x$splicevariant)))
    }
  }
  res
}

.get.modif.pos <- function(x,modif.pos,ptm.info,from) {
   pepseq <- strsplit(x$peptide[1],"")[[1]]
   # get modification position foreach protein (in peptide) from modification string
   modification.positions.foreach.protein <- 
     .convertModifToPos(x$modif,modif.pos,collapse=NULL,simplify=FALSE,and.name=!is.null(names(modif.pos)))

   modif.posi <- t(mapply(.get.modif.pos.for.ac,x$protein,x$splicevariant,modification.positions.foreach.protein,x$start.pos,
                          MoreArgs=list(from=from,pepseq=pepseq,ptm.info=ptm.info)))

   my.gene <- sprintf("%s %s",modif.posi[,3],.string.number.ranges(modif.posi[,4]))
   if (is.na(modif.posi[1,1]) || all(modif.posi[,1]==modif.posi[1,1])) {
     my.modif.posi <- modif.posi[1,1]
   } else {
     my.modif.posi <- paste(modif.posi[,1],collapse=";")
   }


   cbind(modif=unique(x$modif),
         modif.pos=my.modif.posi,
         modif.comment=ifelse(any(!is.na(modif.posi[,2])),
                              paste(modif.posi[!is.na(modif.posi[,2]),2],collapse="\n"),""),
         stringsAsFactors=FALSE)
} 

.get.modif.pos.for.ac <- function(ac,sv,pep.pos.n.modif,start.pos,from,pepseq,ptm.info) {
  gene_name <- proteinInfo(from,protein.ac=ac,select="gene_name",do.warn=FALSE,collapse=",")
  if (length(pep.pos.n.modif) == 0) {
    return(c(NA,NA,ifelse(length(gene_name)==0,NA,gene_name),sv))
  }

  if (is.data.frame(pep.pos.n.modif)) {
    pep.pos <- pep.pos.n.modif[,1]
    modif <- pep.pos.n.modif[,2]
  } else {
    pep.pos <- pep.pos.n.modif
    modif <- NULL
  }
  if (all(pep.pos == 0))
    residue <- 'Nterm'
  else
    residue <- pepseq[pep.pos]

  poss <- start.pos + pep.pos - 1
  if (!is.null(ptm.info) && all(c("isoform_ac","position") %in% colnames(ptm.info))) {
    comments <- sapply(poss,function(pp) {
                       if (grepl("-[0-9]$",ac)) 
                         sel <- ptm.info[,"isoform_ac"]==ac
                       else
                         sel <- ptm.info[,"isoform_ac"]==paste(ac,"-1",sep="")
                       sel  <- sel  & ptm.info[,"position"]==pp

                       if (any(sel)) {
                         res <- apply(ptm.info[sel,],1,
                                      function(pi) paste(sprintf("%s pos %g: %s",pi["isoform_ac"],
                                                                 as.integer(pi["position"]),pi["description"])))
                         paste(res,collapse="\n")
                       } else {
                         NA
                       }
                    })
    #known.pos <- sapply(poss,function(pp) any(ptm.info[,"isoform_ac"]==ac & ptm.info[,"position"]==pp))
    known.pos <- !is.na(comments)
  } else {
    comments <- NA
    known.pos <- rep(FALSE,seq_along(start.pos))
  }
  null.comments <- is.na(comments)

  if (!is.null(modif))
    poss <- paste0(poss,modif)
  if (sum(known.pos) > 0)
    poss[known.pos] <- paste0(residue[known.pos],poss[known.pos],"*")
  poss[!known.pos] <- paste0(residue[!known.pos],poss[!known.pos])

  c(paste(poss,collapse="&"),
    ifelse(all(null.comments),NA,paste(comments[!null.comments],collapse="\n")),
    ifelse(length(gene_name)==0,NA,gene_name),
    sv)
}



proteinID <- function(protein.group,protein.g=reporterProteins(protein.group)) {
  proteinInfo(protein.group,protein.g,do.warn=FALSE,collapse=",")
}
proteinDescription <- function(protein.group,protein.g=reporterProteins(protein.group)) {
  proteinInfo(protein.group,protein.g,select="protein_name",do.warn=FALSE,collapse=",")
}
proteinGeneName <- function(protein.group,protein.g=reporterProteins(protein.group)) {
  proteinInfo(protein.group,protein.g,select="gene_name",do.warn=FALSE,collapse=",")
}

proteinNameAndDescription <- function(protein.group,protein.g=reporterProteins(protein.group),collapse=FALSE) {
  if (collapse)
    data.frame(ID=.paste_unique(proteinID(protein.group,protein.g),collapse=","),
               Description=.paste_unique(proteinDescription(protein.group,protein.g),collapse=";"),
               Gene=.paste_unique(proteinGeneName(protein.group,protein.g),collapse=","),
               stringsAsFactors=FALSE)
  else
    data.frame(ID=proteinID(protein.group,protein.g),
               Description=proteinDescription(protein.group,protein.g),
               Gene=proteinGeneName(protein.group,protein.g),
               stringsAsFactors=FALSE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setGeneric("peptideNProtein",function(x) standardGeneric("peptideNProtein"))
setGeneric("peptideSpecificity", function(x) standardGeneric("peptideSpecificity") )
setGeneric("peptideInfo", function(x) standardGeneric("peptideInfo"))
setGeneric("peptides",function(x,protein,...) standardGeneric("peptides"))
setGeneric("indistinguishableProteins",
           function(x,protein,protein.g) standardGeneric("indistinguishableProteins"))
setGeneric("reporterProteins",function(x, require.reporter.specific=FALSE) standardGeneric("reporterProteins"))
setGeneric("proteinGroupTable",function(x) standardGeneric("proteinGroupTable"))
setGeneric("spectrumToPeptide", function(x) standardGeneric("spectrumToPeptide"))
setGeneric("proteinInfo",function(x,protein.g,protein.ac,...) standardGeneric("proteinInfo"))
setGeneric("proteinInfo<-",function(x,value) standardGeneric("proteinInfo<-"))

setMethod("peptideSpecificity", "ProteinGroup", function(x) x@peptideSpecificity)
setMethod("peptideInfo", "ProteinGroup", function(x) x@peptideInfo)
setMethod("peptideNProtein",   "ProteinGroup", function(x) x@peptideNProtein)
setMethod("proteinGroupTable",      "ProteinGroup", function(x) x@proteinGroupTable )

setMethod("indistinguishableProteins",signature(x="ProteinGroup",protein="missing",protein.g="missing"),
          function(x) x@indistinguishableProteins)
setMethod("indistinguishableProteins",signature(x="ProteinGroup",protein="missing",protein.g="character"),
          function(x,protein.g) {
            sel <- x@indistinguishableProteins %in% protein.g
            if (!any(sel)) {
              warning(protein.g," is no protein group - see reporterProteins",
                      " for a list of all reporter groups")
              return(NA)
            }
            names(x@indistinguishableProteins)[sel]
          }
)
setMethod("indistinguishableProteins",signature(x="ProteinGroup",protein="character",protein.g="missing"),
          function(x,protein) {
            protein.groups <- x@indistinguishableProteins[protein]
            if (length(protein.groups) == 0) {
              warning(protein," is in not in the list")
              return(NA)
            }
            return(protein.groups)
          }
)

setMethod("spectrumToPeptide",   "ProteinGroup", function(x) x@spectrumToPeptide)

setMethod("peptides",signature(x="ProteinGroup",protein="missing"),
    function(x,...) peptides(x,reporterProteins(x),...)
)

setMethod("peptides",signature(x="ProteinGroup",protein="character"),
    function(x,protein,
             specificity=c(UNSPECIFIC,GROUPSPECIFIC,REPORTERSPECIFIC),
             columns=c('peptide'),set=union,drop=TRUE,
             groupspecific.if.same.ac=FALSE,do.warn=TRUE,
             modif=NULL) {

      if (!all(specificity %in% c(UNSPECIFIC,GROUPSPECIFIC,REPORTERSPECIFIC)))
        stop("specificity should be one or more of ['",REPORTERSPECIFIC,"','",GROUPSPECIFIC,"','",REPORTERSPECIFIC,"']")

      pnp <- peptideNProtein(x)
      ps <- peptideSpecificity(x)
      pi <- merge(peptideInfo(x),peptideSpecificity(x),by="peptide")

      if (groupspecific.if.same.ac) {
        group.proteins <- subset(proteinGroupTable(x),reporter.protein==protein,"protein.g",drop=TRUE)
        group.proteins.all <- indistinguishableProteins(x,protein.g=group.proteins)
        protein.acs <- subset(x@isoformToGeneProduct,
                              proteinac.w.splicevariant %in% group.proteins.all,
                              "proteinac.wo.splicevariant",drop=TRUE)
        if (length(unique(protein.acs)) == 1)
          specificity <- c(specificity,GROUPSPECIFIC)
      }
      
      peptides <- lapply(protein, function(p) pnp[pnp[,"protein.g"] == p,"peptide"])
      peptides <- Reduce(set,peptides)

      sel <- pi$peptide %in% peptides                                          # is in protein?
      sel <- sel & pi$peptide %in% ps$peptide[ps$specificity %in% specificity] # has specificity?
      if (!is.null(modif)) sel <- sel & .has.modif(pi$modif,modif)             # has modification?

      peptides <- pi[sel,columns,drop=drop]
      if (length(peptides) == 0 && do.warn)
        warning("No peptide for protein ",protein," with specificity ",paste(specificity,collapse=","))
      return(unique(peptides))
    }
)

.has.modif <- function(modifstring,modif) {
  sapply(strsplit(modifstring,":"),function(m) any(m %in% modif))
}

setMethod("reporterProteins","ProteinGroup",
    function(x,require.reporter.specific=FALSE) {
      if (is.null(require.reporter.specific) || !require.reporter.specific) {
        return(x@proteinGroupTable$protein.g[x@proteinGroupTable$is.reporter])
      } else {
        return(x@proteinGroupTable$protein.g[
                                             x@proteinGroupTable$is.reporter &
                                             x@proteinGroupTable[,'n.reporter-specific'] > 0 ])
      }
    }
)

setMethod("proteinInfo",signature(x="ProteinGroup",protein.g="missing",protein.ac="missing"),function(x) x@proteinInfo)
setMethod("proteinInfo",signature(x="ProteinGroup",protein.g="missing",protein.ac="character"),
    function(x,protein.ac,select="name",collapse=", ",simplify=TRUE,do.warn=TRUE) {
      protein.info <- proteinInfo(x)
      if (length(protein.info) == 0) {
        if (isTRUE(do.warn))
          warning("proteinInfo is empty")
        return()
      }
      if (!all(select %in% colnames(protein.info))) 
        warning("Column(s) ",
                paste(select[!select %in% colnames(protein.info)],collapse=", "),
                " not in proteinInfo (available:",paste(colnames(protein.info),collapse=", "))

      if ((is.null(protein.info) || nrow(protein.info) == 0) && do.warn)
        warning("protein info is NULL! Set for ProteinGroup.")

      res <- sapply(protein.ac,function(p) {
        if (!all(select %in% colnames(protein.info))) 
          return(NA)
        
        if (proteinInfoIsOnSpliceVariants(protein.info))
          protein.acs <- x@isoformToGeneProduct[protein.ac,"proteinac.wo.splicevariant"]
        sel <- protein.info$accession %in% protein.ac
        if (!any(sel)) {
          if (do.warn) warning("No protein info for ",p)
          return(rep(NA,length(select)))
        }
        protein.infos <- protein.info[sel,select]
        if (simplify) {
          if (length(select) > 1)
            apply(protein.infos,2,function(pi) paste(sort(unique(pi)),collapse=", "))
          else 
            paste(sort(unique(protein.infos)),collapse=", ")
        } else {
      if (length(select) == 1)
            names(protein.infos) <- protein.info$accession[sel]
          protein.infos
        }
      },simplify=simplify)
      if (simplify & length(select)>1) t(res)
      else res
    }
)

setMethod("proteinInfo",signature(x="ProteinGroup",protein.g="character",protein.ac="missing"),
    function(x,protein.g,select="name",collapse=", ",simplify=TRUE,do.warn=TRUE) {
      protein.info <- proteinInfo(x)
      if (length(protein.info) == 0) {
        if (isTRUE(do.warn))
          warning("proteinInfo is empty")
        return()
      }

      if (!all(select %in% colnames(protein.info))) 
        warning("Column(s) ",
                paste(select[!select %in% colnames(protein.info)],collapse=", "),
                " not in proteinInfo (available:",paste(colnames(protein.info),collapse=", "))

      if ((is.null(protein.info) || nrow(protein.info) == 0) && do.warn)
        warning("protein info is NULL! Set for ProteinGroup.")

      res <- sapply(protein.g,function(p) {
        if (!all(select %in% colnames(protein.info))) 
          return(NA)
        
        if (proteinInfoIsOnSpliceVariants(protein.info))
          protein.acs <- indistinguishableProteins(x,protein.g=p)
        else
          protein.acs <- x@isoformToGeneProduct[indistinguishableProteins(x,protein.g=p),"proteinac.wo.splicevariant"]
        sel <- protein.info$accession %in% protein.acs
        if (!any(sel)) {
          if (do.warn) warning("No protein info for ",p)
          return(rep(NA,length(select)))
        }
        protein.infos <- protein.info[sel,select]
        if (simplify) {
          if (length(select) > 1)
            apply(protein.infos,2,function(pi) paste(sort(unique(pi)),collapse=", "))
          else 
            paste(sort(unique(protein.infos)),collapse=", ")
        } else {
      if (length(select) == 1)
            names(protein.infos) <- protein.info$accession[sel]
          protein.infos
        }
      },simplify=simplify)
      if (simplify & length(select)>1) t(res)
      else res
    }
)
setReplaceMethod("proteinInfo","ProteinGroup",
                 function(x,value) {
                   x@proteinInfo <- value
                   if ('sequence' %in% colnames(value) && all(is.na(x@peptideInfo[,'start.pos']))) {
                     message("Recaculating peptide start position based on sequence")
                     x@peptideInfo <- calcPeptidePosition(x@peptideInfo,value)
                   }
                   x
                 })

setGeneric("reporter.protein",function(x,protein.g) standardGeneric("reporter.protein"))
setMethod("reporter.protein",signature("ProteinGroup","character"),
          function(x,protein.g) {
            proteinGroupTable(x)[proteinGroupTable(x)[,'protein.g']==protein.g, 'reporter.protein']
          })

setGeneric("protein.g",function(x,pattern,variables=c("AC","name"),...) standardGeneric("protein.g"))
setMethod("protein.g",signature("ProteinGroup","character","ANY"),
          function(x,pattern,variables=c("AC","name"),...) {
  ip <- indistinguishableProteins(x)
  protein.info <- proteinInfo(x)
  result <- c()
  for (p in pattern) {
    protein.acs <- c()
    if ("AC" %in% variables)
      result <- c(result,ip[grep(p,names(ip),...)])
    if ("name" %in% variables) {
      if (length(protein.info) != 0L) {
        protein.acs <- unique(c(
           protein.info$accession[grep(p,protein.info$name,...)],
           protein.info$accession[grep(p,protein.info$gene_name,...)],
           protein.info$accession[grep(p,protein.info$protein_name,...)]
        ))
      }
      
      i.to.gp <- x@isoformToGeneProduct
      sel.isoforms <- i.to.gp[,"proteinac.wo.splicevariant"] %in% protein.acs
      protein.acs.w.isoforms <-
         i.to.gp[sel.isoforms,"proteinac.w.splicevariant"]
      ## TODO: test for proteinInfoIsOnSpliceVariants
      protein.gs <- ip[names(ip) %in% c(protein.acs.w.isoforms,protein.acs)]
      result <- c(result,protein.gs)
   }
 }
 result <- unique(as.character(result))
 if (length(result) == 0)
   warning("Could not find protein group identifier for ",pattern)
 return(result)
})

setGeneric("protein.ac",function(x,protein.g) standardGeneric("protein.ac"))
setMethod("protein.ac",signature("ProteinGroup","character"),
  function(x,protein.g) {
   unique(x@isoformToGeneProduct[indistinguishableProteins(x,protein.g=protein.g),"proteinac.wo.splicevariant"])
}) 

setMethod("protein.ac",signature("ProteinGroup","missing"),
  function(x) {
   unique(x@isoformToGeneProduct[,"proteinac.wo.splicevariant"])
})



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show","ProteinGroup",function(object) {
      cat("ProteinGroup object\n")
      cat(sprintf("%8d peptides are detected\n",nrow(peptideSpecificity(object))))
      cat(sprintf("%8d protein groups with specific peptides\n",length(reporterProteins(object,require.reporter.specific=TRUE))))
    }
)

get.pep.group <- function(x,protein) {
  pep.specificity <- peptideSpecificity(x)
  protein.group.table <- proteinGroupTable(x)

  message(protein)
  speci <- rep(3,length(pep.specificity$peptide))
  speci[pep.specificity$specificity==GROUPSPECIFIC] <-  2
  speci[pep.specificity$specificity==REPORTERSPECIFIC] <- 1
  names(speci) <- pep.specificity$peptide

  protein.group <- protein.group.table[protein.group.table$reporter.protein==protein,]
  pepnprots <- peptideNProtein(x)
  pp <- as.data.frame(pepnprots[pepnprots[,'protein.g'] %in%
                                protein.group$protein.g,,drop=FALSE],stringsAsFactors=FALSE)

  # get all peptides of protein
  peptides <- unique(pp[,'peptide'])
  t <- table(pp$protein.g)
  proteins <- names(sort(t,decreasing=TRUE))

  # to which proteins are peptides assigned?
  pep.groups <-
    ddply(pp[pp[,'protein.g'] %in% proteins,],"peptide",
          function(s) data.frame(peptide=unique(s$peptide),
                                 protein.g=s['protein.g'],
                                 group=paste(sort(unique(s$protein.g)),collapse="-"),
                                 speci=speci[unique(s$peptide)],
                                 stringsAsFactors=FALSE
                                 )
          )
  pep.groups <- unique(pep.groups[,c("peptide","group","speci")])
  pep.groups[order(pep.groups$speci,pep.groups$peptide),]
}

## create a extended proteinGroup table
## TODO: check function usage
.extend.protein.group <- function(ibspectra,calc.ratio) {
  pg <- proteinGroup(ibspectra)
  protein.group <- proteinGroupTable(pg)
  pnp <- peptideNProtein(pg)
  t <- table(pnp[,"peptide"])
  pnp <- pnp[pnp[,"peptide"] %in% names(t)[t==2],]

  protein.group$lratio <- NA
  protein.group$var    <- NA
  protein.group$is.different    <- FALSE
  protein.group$is.smaller    <- 0
  protein.group$is.bigger    <- 0
  for (i in seq_len(nrow(protein.group))) {
    if (protein.group[i,"is.reporter"]) {
      reporter.protein <- protein.group[i,"protein.g"]
      reporter.ratio <- calc.ratio(ibspectra,protein=reporter.protein)
      protein.group[i,c("lratio","var")] <- reporter.ratio[c("lratio","variance")]
    } else {
      ## reporter.ratio set above - protein.group is sequential
      ## subset protein: get peptide shared with (and only with) reporter
      shared.peptides <- pnp[pnp[,"protein.g"]==protein.group[i,"protein.g"],"peptide"]
      if (length(shared.peptides)>0) {
        shared.ratio <- calc.ratio(ibspectra,peptide=shared.peptides)
        protein.group[i,c("lratio","var")] <- shared.ratio[c("lratio","variance")]
      
        shared.issmaller <- shared.ratio["lratio"] + sqrt(shared.ratio["variance"]) < 
            reporter.ratio["lratio"] - sqrt(reporter.ratio["variance"])
        shared.isbigger <- shared.ratio["lratio"] - sqrt(shared.ratio["variance"]) > 
            reporter.ratio["lratio"] + sqrt(reporter.ratio["variance"])

        protein.group[i,"is.different"] <- shared.issmaller | shared.isbigger
        protein.group[i,"is.smaller"] <- shared.issmaller
        protein.group[i,"is.bigger"] <- shared.isbigger

      }
    }
  }
  return(protein.group)
}


groupMemberPeptides <- function(x,reporter.protein.g,
                                ordered.by.pos=TRUE,only.first.pos=TRUE) {
   peptide.info <- unique(x@peptideInfo[,c("protein","peptide","start.pos")])
   group.table <- proteinGroupTable(x)
   indist.proteins <- x@indistinguishableProteins

   reporter.proteins <- names(indist.proteins)[indist.proteins==reporter.protein.g]
   my.group.table <- group.table[group.table$reporter.protein==reporter.protein.g,]

   # 1. get reporter peptides (w/ position and specificity)
   my.peptide.info <- peptides(x,reporter.protein.g,
                               columns=c("peptide","specificity",
                                         "n.shared.groups","n.shared.proteins"),drop=FALSE,do.warn=FALSE)
 
   if (ordered.by.pos) {     
      reporter.peptide.startpos <- 
         peptide.info[peptide.info$peptide %in% my.peptide.info$peptide &
         peptide.info$protein %in% reporter.proteins[1],c("peptide","start.pos")]

      if (only.first.pos) {
         my.peptide.info$start.pos <- 0
         for (peptide.i in seq_len(nrow(my.peptide.info))) {
            sel <- reporter.peptide.startpos$peptide==my.peptide.info$peptide[peptide.i]
            my.peptide.info$start.pos[peptide.i] <- sort(reporter.peptide.startpos$start.pos[sel])[1]
         }
      } else {
          my.peptide.info <- merge(reporter.peptide.startpos,my.peptide.info)
      }
      my.peptide.info <- my.peptide.info[order(my.peptide.info$start.pos),]
   } 
   rownames(my.peptide.info) <- NULL

   # 2. get group member peptides
   group.member.peptides <- do.call(cbind,lapply(my.group.table$protein.g,
      function(protein.g) my.peptide.info$peptide %in% peptides(x,protein.g,do.warn=FALSE)))
   rownames(group.member.peptides) <- my.peptide.info$peptide
   colnames(group.member.peptides) <- my.group.table$protein.g
   return(list(peptide.info=my.peptide.info,group.member.peptides=group.member.peptides))
}


human.protein.names <- function(my.protein.info) {
  my.df <- my.protein.info[,c("accession","splicevariant","gene_name","protein_name")]
  my.df$gene_name <- sanitize(my.df$gene_name)
  collapsed.splicevariant <- ddply(my.df,"accession",function(x) {
        only_one <- nrow(x) == 1
        if (!all(is.na(x$splicevariant)) & any(!is.na(x$splicevariant)))
          x$splicevariant[is.na(x$splicevariant)] <- 1

        x$splicevariant <- number.ranges(x$splicevariant)
        x <- unique(x)
        if (nrow(x) > 1) 
          x <- as.data.frame(lapply(x,function(y) paste(y[!is.na(y)],collapse=",")),stringsAsFactors=FALSE)
        
        if (is.na(x$splicevariant)) {
          x$ac_link <- sprintf("\\uniprotlink{%s}",sanitize(x$accession,dash=FALSE))
          x$ac_nolink <-  x$accession
        } else {
          x$ac_link <- sprintf("\\uniprotlink{%s}$_{%s}$",sanitize(x$accession,dash=FALSE),x$splicevariant);
          if (only_one) {
            x$ac_nolink <- sprintf("%s-%s",x$accession,x$splicevariant)
          } else {
            x$ac_nolink <- sprintf("%s-[%s]",x$accession,x$splicevariant)
          }
        }
        return(unique(x))
    })
    
    collapsed.gene_name <- ddply(
          unique(collapsed.splicevariant[,c("gene_name","ac_link","ac_nolink","protein_name")]),
          "gene_name",function(x) {
        if (is.na(x$gene_name[1]) || x$gene_name[1] == "") {
          x$name_nolink <- x$ac_nolink
        } else {
          x$ac_link_bold <- sprintf("\\textbf{%s}: %s",unique(x$gene_name),x$ac_link)
          x$ac_link <- sprintf("%s {\\small %s}",unique(x$gene_name),x$ac_link)
          x$ac_nolink <- sprintf("%s: %s",unique(x$gene_name),x$ac_nolink)
          x$name_nolink <- sprintf("%s: %s",x$gene_name,x$protein_name)
        }
        return(unique(x))
    })
    return(collapsed.gene_name)
}


my.protein.info <- function(x,protein.g) {
    protein.acs <- indistinguishableProteins(x,protein.g=protein.g)
    if (all(is.na(protein.acs))) return()
    isoforms <- x@isoformToGeneProduct
    res <- data.frame(protein.ac=protein.acs,
                      accession=isoforms[protein.acs,"proteinac.wo.splicevariant"],
                      splicevariant=as.integer(isoforms[protein.acs,"splicevariant"]), 
                      stringsAsFactors=FALSE)

    if (length(proteinInfo(x)) > 0) {
      if (proteinInfoIsOnSpliceVariants(proteinInfo(x))) {
        res <- merge(res,proteinInfo(x),by.x="protein.ac",by.y="accession",all.x=TRUE,all.y=FALSE)
      } else {
        res <- merge(res,proteinInfo(x),by="accession",all.x=TRUE,all.y=FALSE)
      }
    } else {
      res <- cbind(res,gene_name=NA,protein_name=NA)
    }

    return(.moveToFirstCol(res,"accession"))
}

summary.ProteinGroup <- function(object,only.reporters=TRUE,...) {
  peptide.cnt <- table(peptideNProtein(object)[,"protein.g"])
  peptide.spectra.cnt <- table(spectrumToPeptide(object))
  spectra.cnt <- tapply(peptideNProtein(object)[,"peptide"],peptideNProtein(object)[,"protein.g"],function(pep) sum(peptide.spectra.cnt[pep]))
  
  if (only.reporters) 
    proteins <- reporterProteins(object)
  else
    proteins <- names(spectra.cnt)
  
  proteins <- proteins[order(spectra.cnt[proteins],decreasing=TRUE)]
    
  return(data.frame(peptide.cnt=peptide.cnt[proteins],
                    spectra.cnt=spectra.cnt[proteins]))
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quantification functions (dNSAF).
###  based on peptide or spectral count

peptide.count <- function(protein.group,protein.g=reporterProteins(protein.group),
                          specificity=c("reporter-specific","group-specific","unspecific"),...) {
  sapply(protein.g,
         function(p) length(peptides(protein.group,p,specificity=specificity,...)))
}

spectra.count2 <- function(ibspectra,value=reporterProteins(protein.group),type="protein.g",
                          specificity=c("reporter-specific","group-specific","unspecific"),
                          modif=NULL,combine=FALSE,subset=NULL,require.quant=NULL,...) {
  protein.group <- proteinGroup(ibspectra)
  if (!isTRUE(combine)) {
    spectra.count <- sapply(value, function(p) 
                            spectra.count2(ibspectra,p,type,specificity,modif,combine=TRUE,
                                           subset=subset,require.quant=require.quant,...))
    names(spectra.count) <- value
    return(spectra.count)
  }

  pep <- switch(type,
                protein.g=peptides(protein.group,protein=value,specificity=specificity,...),
                peptide=value)

  ## Calculate unique spectrum counts for all proteins
  if (!is.null(subset)) {
    fd <- subset(fData(ibspectra),eval(subset) & peptide %in% pep)
  } else {
    fd <- subset(fData(ibspectra),peptide %in% pep)
  }

  if (!is.null(require.quant)) {
    spectra <- rownames(reporterIntensities(ibspectra,na.rm=TRUE,na.rm.f=require.quant))
    fd <- fd[fd$spectrum %in% spectra,]
  }

  if (is.null(modif)) {
    peptide.spectra.count <- table(fd$peptide)
    sum(peptide.spectra.count)
  } else {
    has.modif <- sapply(strsplit(fd$modif,":"),function(x) any(x %in% modif))
    sum(has.modif)
  }
}

spectra.count <- function(protein.group,protein.g=reporterProteins(protein.group),
                          specificity=c("reporter-specific","group-specific","unspecific"),
                          modif=NULL,...) {
  ## Calculate unique spectrum counts for all proteins
  if (is.null(modif)) {
    peptide.spectra.count <- table(spectrumToPeptide(protein.group))
    spectra.count <- sapply(protein.g, function(p) 
                            sum(peptide.spectra.count[peptides(protein.group,protein=p,
                                                               specificity=specificity,...)]))
  } else {
    si <- protein.group@spectrumId
    has.modif <- sapply(strsplit(si$modif,":"),function(x) any(x %in% modif))
    spectra.count <- sapply(protein.g, function(p) {
                            pep <- peptides(protein.group,protein=p,specificity=specificity,...)
                            sum(si$peptide %in% pep & has.modif) })
  }

  names(spectra.count) <- protein.g
  return(spectra.count)
}

calcPeptidePosition <- function(peptide.info,protein.info,calc.il.peptide) {
  peptide.info[,'start.pos'] <- NA
  if (missing(calc.il.peptide)) 
    calc.il.peptide <- !('real.peptide' %in% colnames(peptide.info) && !all(is.na(peptide.info[,'real.peptide'])))
  if (isTRUE(calc.il.peptide)) 
    peptide.info[,'real.peptide'] <- NA
  peptide.info <- unique(peptide.info)
  last.elem <- nrow(peptide.info)

  for (p.i in seq_len(nrow(peptide.info))) {
    ## WARNING: Does not yet work on protein groups which have no splice variant specified!
    seqq <- protein.info[protein.info[,'accession']==peptide.info[p.i,'protein'],'sequence'][1]
    pep.i <- peptide.info[p.i,]
    pep <- peptide.info[p.i,'peptide']
    if (isTRUE(getOption("isobar.verbose")))
      message(pep,appendLF=FALSE)
    if (is.na(pep) || nchar(pep) == 0) stop("peptide on position ",p.i," is NA or of zero length")
    if (length(seqq) == 0 || is.na(seqq) || nchar(seqq) == 0) {
      peptide.info[p.i,'real.peptide'] <- paste(pep,"(not matched, no protein sequence)")
      if (isTRUE(getOption("isobar.verbose")))
      message(": no prot sequence")
      next
    } 

    il.seqq <- gsub("I","L",seqq) ## to be confomant with converted peptides
    pep <- gsub("I","L",pep)
    matches <- gregexpr(pep,il.seqq,fixed=TRUE)[[1]]

    pp.i <- p.i
    for (m.i in seq_along(matches)) {
       if (m.i > 1) {
         last.elem <- last.elem + 1
         peptide.info <- rbind(peptide.info,pep.i)
         pp.i <- last.elem
       }
       if (matches[m.i] == -1) {
         peptide.info[pp.i,'start.pos'] <- -1
         peptide.info[pp.i,'real.peptide'] <- paste(pep,"(not matched)")
       } else {
         peptide.info[pp.i,'start.pos'] <- matches[m.i]
         peptide.info[pp.i,'real.peptide'] <- 
                 substr(seqq,matches[m.i],matches[m.i]+attr(matches,"match.length")[m.i]-1)
       }
    }
    if (isTRUE(getOption("isobar.verbose")))
      message(" --> ",peptide.info[pp.i,'real.peptide'])
  }
  if (any(peptide.info[,'start.pos'] == -1,na.rm=TRUE)) {
    sel.bad <- !is.na(peptide.info[,'start.pos']) & peptide.info[,'start.pos'] == -1
    warning("Could not match ",length(unique(peptide.info[sel.bad,'peptide'])),
            " peptides in ",length(unique(peptide.info[sel.bad,'protein']))," proteins")
    peptide.info[sel.bad,'start.pos'] <- NA
  }
  
  peptide.info[order(peptide.info[,'protein'],peptide.info[,'start.pos'],peptide.info[,'peptide']),]
}

# TOFIX: proteinInfo does not contain splice information. With splice sequence, an accurate seqcov could be calculated
sequence.coverage <- function(protein.group,protein.g=reporterProteins(protein.group),
                              specificity=c("reporter-specific","group-specific","unspecific"),
                              simplify=TRUE,...) {

  if (!proteinInfoIsOnSpliceVariants(proteinInfo(protein.group)))
    warning("Protein information is not on splice variants - sequence coverage will be approximate only.")


  if ("length" %in% colnames(proteinInfo(protein.group)) && 
      "start.pos" %in% colnames(protein.group@peptideInfo)) {

    if (all(is.na(protein.group@peptideInfo[,'start.pos']))) {
      protein.group@peptideInfo <- calcPeptidePosition(protein.group@peptideInfo,proteinInfo(protein.group))
    }
    lengths <- proteinInfo(protein.group,protein.g=protein.g,select="length",simplify=FALSE)
    peptides <- peptides(protein.group,protein=protein.g,specificity=specificity,...)
    peptide.info <- subset(unique(protein.group@peptideInfo[,c("protein","peptide","start.pos")]),
                           peptide %in% peptides)
    protein.ac.wo.splice <- isobar:::.as.vect(protein.group@isoformToGeneProduct)

    res <- sapply(protein.g,function(p) {
             seqcov <- sapply(indistinguishableProteins(protein.group,protein.g=p),function(pp) {
               if (proteinInfoIsOnSpliceVariants(proteinInfo(protein.group)))
                 ppp <- pp
               else
                 ppp <- protein.ac.wo.splice[pp]
               .calc.seqcov(lengths[[p]][as.character(ppp)],
                          subset(peptide.info,protein==pp))
             })
             if (isTRUE(simplify))
               .simplify.seqcov(seqcov)
             else
               seqcov
    })
    names(res) <- protein.g
    res
  } else {
    stop("Necessary information for sequence coverage not available")
  }
}

.calc.seqcov <- function(length,peptide.info) {
      if (is.na(length) || all(is.na(peptide.info$start.pos)))
        return(NA)
      seqq <- rep(FALSE,length)
      peptide.info <- peptide.info[!is.na(peptide.info$start.pos),]
      peptide.info$peplength <- nchar(peptide.info$peptide)
      peptide.info$end.pos <- peptide.info$start.pos+peptide.info$peplength-1
      if (any(peptide.info[,'end.pos']>length)) {
        warning(".calc.seqcov: peptides present which span over the protein length, removing them")
        peptide.info <- peptide.info[peptide.info[,'end.pos'] <= length,]
      }
      for (i_r in seq_len(nrow(peptide.info)))
        seqq[seq(from=peptide.info[i_r,"start.pos"],to=peptide.info[i_r,"end.pos"])]  <- TRUE
      return(sum(seqq)/length(seqq))
}

.simplify.seqcov <- function(seqcov) {
      if (length(seqcov) > 1)
        mean(seqcov,na.rm=TRUE)
      else
        seqcov
}


calculate.dNSAF <- function(protein.group,use.mw=FALSE,normalize=TRUE,combine.f=mean) {
  if (is.null(proteinInfo(protein.group)) || length(proteinInfo(protein.group)) == 0)
    stop("slot proteinInfo not set - it is needed for sequence length")
  if (!"length" %in% colnames(proteinInfo(protein.group)))
    stop("no column 'length' in proteinInfo slot")

  ip <- indistinguishableProteins(protein.group)
  if(proteinInfoIsOnSpliceVariants(proteinInfo(protein.group)))
    get.to.acs <- function(my.protein.g) names(ip)[ip==my.protein.g]
  else
    get.to.acs <- function(my.protein.g) unique(protein.ac(protein.group,my.protein.g))

  seqlength <- proteinInfo(protein.group)$length
  names(seqlength) <- proteinInfo(protein.group)$accession
  
  spectrum.counts <- table(spectrumToPeptide(protein.group))
  ## Calculate unique spectrum counts for all proteins
  uspc <- sapply(reporterProteins(protein.group), function(p) 
                 sum(spectrum.counts[peptides(protein.group,protein=p,
                                              specificity=c("reporter-specific","group-specific"))]))
  names(uspc) <- reporterProteins(protein.group)
  
  ## Calculate distribution factors for each shared peptide
  pnp <- peptideNProtein(protein.group)
  pnp <- pnp[pnp[,"protein.g"] %in% reporterProteins(protein.group),]
  t <- table(pnp[,"peptide"])
  shared.peptides <- names(t)[t>1]
  distr.factors <- lapply(shared.peptides,function(pep) {
    proteins <- pnp[pnp[,"peptide"]==pep,"protein.g"]
    uspc.p <- uspc[proteins]
    uspc.p/sum(uspc.p)
  })
  names(distr.factors) <- shared.peptides

  ## Calculate dSAF
  dSAF <- sapply(reporterProteins(protein.group), function(p) {
    shared.pep.p <- peptides(protein.group,protein=p,specificity="unspecific",do.warn=FALSE)
    if (length(shared.pep.p) > 0)
      d.sspc <- sum(sapply(shared.pep.p, function(pep) spectrum.counts[pep]*distr.factors[[pep]][p]))
    else
      d.sspc <- 0

    ## for protein length, we take the sum of lengths of all acs
    protein.length <- sum(seqlength[get.to.acs(p)],na.rm=TRUE)
    dSAF <- (uspc[p] + d.sspc) / protein.length
    if (is.nan(dSAF) | !is.finite(dSAF)) dSAF <- NA
    return(dSAF)
  })

  if (use.mw) {
    mw <- .calculate.mw(protein.group,protein.g,combine.f)
    dSAF <- dSAF*mw
  }

  ## normalize to dNSAF
  if (normalize) dNSAF <- dSAF/sum(dSAF,na.rm=TRUE)

  names(dNSAF) <- reporterProteins(protein.group)
  return(dNSAF)
}

## get.seqlength returns amino acid sequence length of Uniprot ACs
.DEPR.get.seqlength <- function(my.acs,uniprot=NULL,uniprot.acs=NULL) { 
  if (is.null(uniprot)) {
    ## Download Uniprot/Swissprot from:
    ## ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

    ## Alternative: use ANUC
    ## choosebank("swissprot")
    ## cat(banknameSocket$details, sep = "\n")

    #library(seqinr)
    uniprot <- read.fasta("uniprot_sprot.fasta",seqtype="AA",as.string=TRUE)
  }
  if (is.null(uniprot.acs))
    uniprot.acs <- do.call(rbind,strsplit(names(uniprot),"|",fixed=TRUE))[,2]
  seqlength <- sapply(uniprot[uniprot.acs %in% my.acs],nchar)
  names(seqlength) <- uniprot.acs[uniprot.acs %in% my.acs]
  return(seqlength)
}


.DEPR.n.obs.peptides <- function(seq,nmc=1,min.length=6) {
  seq.c <- strsplit(seq,"")[[1]]
  pos <- gregexpr("[KR]",seq,perl=TRUE)[[1]]
  pos <- pos[seq.c[pos+1] != "P"]
  pos <- pos[!is.na(pos)]
  
  pep <- sapply(seq_along(c(pos,1)),function(i) {
    idx.from <- ifelse(i == 1,1,pos[i-1] + 1)
    idx.to <- ifelse(i > length(pos), length(seq.c), pos[i])
    
    if (idx.to-idx.from >= min.length)
      substr(seq,idx.from,idx.to)[[1]]
    else
      NA
  })
  pep <- pep[!is.na(pep)]

  for (i in seq(len=nmc)) {
    nmcx <- sapply(seq_along(pos),function(i) {
      idx.from <- ifelse(i == 1,1,pos[i-1] + 1)
      idx.to <- ifelse(i == length(pos), length(seq.c), pos[i+1])
      if (idx.to-idx.from >= min.length)
        substr(seq,idx.from,idx.to)
      else
        NA
    })
    pep <- unique(c(pep,nmcx[!is.na(nmcx)]))
  }

  return(length(pep))
}

n.observable.peptides <- function(...) {
  return(nrow(observable.peptides(...)))
}

# crude calculations:
#   mass of B: (mass_N+mass_D)/2, mass_N=114.04293; mass_D=215.06680
#   mass of Z: (mass_E+mass_Q)/2, mass_E=229.08245; mass_Q=328.13828
#   mass of J: mass of L
observable.peptides <- function(seq,nmc=1,min.length=6,min.mass=600,max.mass=4000,
                                custom=list(code=c("B","Z","J","U"),
                                            mass=c(164.554862,278.61037,213.12392,150.953636)),...) {
  if (is.na(seq) || length(seq)==0 || nchar(seq) == 0)
    return(0)
  seq <- gsub("[ ,]","",seq)
  require(OrgMassSpecR)
  pep <- Digest(seq,missed=nmc,custom=custom,...)
  min.length.ok <- nchar(pep[,"peptide"]) >= min.length
  mass.ok <- pep[,"mz1"] >= min.mass & pep[,"mz3"] <= max.mass
  pep[min.length.ok & mass.ok,]
}

.calculate.mw <- function(protein.group,protein.g,combine.f=mean) {
  require(OrgMassSpecR)
  sequences <- proteinInfo(protein.group)$sequence
  ip <- indistinguishableProteins(protein.group)
  names(sequences) <- proteinInfo(protein.group)$accession

  if(proteinInfoIsOnSpliceVariants(proteinInfo(protein.group)))
    get.to.acs <- function(my.protein.g) names(ip)[ip==my.protein.g]
  else
    get.to.acs <- function(my.protein.g) unique(protein.ac(protein.group,my.protein.g))

  sapply(protein.g, function(my.protein.g) 
         do.call(mean,lapply(get.to.acs(my.protein.g),
                                  function(protein.ac) MolecularWeight(formula = ConvertPeptide(sequences[protein.ac])))))
}

calculate.emPAI <- function(protein.group,protein.g=reporterProteins(protein.group),
                            normalize=FALSE,observed.pep=c("pep","mod.charge.pep"),
                            use.mw=FALSE,combine.f=mean,...,nmc=0,report.all=FALSE) {
  if (is.null(proteinInfo(protein.group)) || length(proteinInfo(protein.group)) == 0)
    stop("slot proteinInfo not set - it is needed for sequence length")
  if (!"sequence" %in% colnames(proteinInfo(protein.group))) {
    stop("no column 'sequence' in proteinInfo slot")
  }
  ip <- indistinguishableProteins(protein.group)

  sequences <- proteinInfo(protein.group)$sequence
  names(sequences) <- proteinInfo(protein.group)$accession
 
  n.observed.peptides <- switch(observed.pep[1],
                                pep = table(peptideNProtein(protein.group)[,"protein.g"])[protein.g],
                                mod.charge.pep = {
                                  pep.cnt <- table(unique(protein.group@spectrumId[,c("peptide","modif","charge")])[,"peptide"])
                                  tapply(peptideNProtein(protein.group)[,"peptide"],peptideNProtein(protein.group)[,"protein.g"],function(x) sum(pep.cnt[x]))[protein.g]
                                })
  if(proteinInfoIsOnSpliceVariants(proteinInfo(protein.group)))
    get.to.acs <- function(my.protein.g) names(ip)[ip==my.protein.g]
  else
    get.to.acs <- function(my.protein.g) unique(protein.ac(protein.group,my.protein.g))

  n.observable.peptides <- sapply(protein.g, function(my.protein.g) 
                                  do.call(combine.f,lapply(get.to.acs(my.protein.g),
                                                           function(protein.ac) n.observable.peptides(sequences[protein.ac],nmc=nmc,...))))
  if (use.mw) 
    mw <- .calculate.mw(protein.group,protein.g,combine.f)
   else 
    mw <- 1

  empai <- 10^(n.observed.peptides/n.observable.peptides) - 1
  if (report.all) 
    return(data.frame(protein.g,n.observed.peptides,n.observable.peptides,mw,
                      empai,protein_mol=empai/sum(empai,na.rm=TRUE),protein_weight=empai*mw/sum(empai*mw,na.rm=TRUE)))

  if (use.mw) empai <- empai*mw
  if (normalize) empai/sum(empai,na.rm=TRUE)

  return(empai)
}

# pretty format group identifiers
#   Input: c("P1234-1,P1234-2","P6543")
#  Output: c("P1234-[1,2]","P6543")
.protein.acc <- function(my.protein.g,protein.group) {

  if (length(my.protein.g) > 1) 
    return(sapply(my.protein.g,function(p) .protein.acc(p,protein.group)))

  protein.acs <- as.character(indistinguishableProteins(protein.group,protein.g=my.protein.g))

  if (length(protein.acs) == 1) return(protein.acs)
  splice.df <- protein.group@isoformToGeneProduct[protein.acs,]
  if (all(is.na(splice.df[,'splicevariant']))) 
    return(paste0(splice.df[,'proteinac.w.splicevariant'],collapse=", "))

  res <- ddply(splice.df,"proteinac.wo.splicevariant",function(y) {
        if (nrow(y) == 1) return(y[,'proteinac.w.splicevariant'])
        if (all(is.na(y$splicevariant))) return(y[1,'proteinac.w.splicevariant'])

        return(sprintf("%s-[%s]",
                       y[1,'proteinac.wo.splicevariant'],
                       number.ranges(y[,'splicevariant'])))
  })
  return(as.character(paste0(res[,'V1'],collapse=", ")))
}


