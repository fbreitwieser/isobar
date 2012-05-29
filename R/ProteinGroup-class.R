### =========================================================================
### ProteinGroup objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("ProteinGroup",
    contains = "VersionedBiobase",
    representation(
                   spectrumToPeptide = "character",
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
        new("VersionedBiobase",versions=c(ProteinGroups="1.0.0")),
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
          function(from,template,proteinInfo) {      
      # TODO: exclude groups from beeing reporters when
      #       there's no reporter-specific peptide ?
      # TODO: add isoformToGeneProduct
      if (ncol(from) == 3) {
          colnames(from) <- c("spectrum","peptide","protein")
      }
      if (ncol(from) == 4) {
          colnames(from) <- c("spectrum","peptide","startpos","protein")
      }
      peptideNProtein <- peptideNProtein(template)[peptideNProtein(template)[,"peptide"] %in% from[,'peptide'],]
      protein.group.table <-
        subset(proteinGroupTable(template),protein.g %in% peptideNProtein[,"protein.g"])
      ipt <- indistinguishableProteins(template)
      indistinguishableProteins <- ipt[ipt %in% peptideNProtein[,"protein.g"]]
      peptideSpecificity <-
        subset(peptideSpecificity(template),peptide %in% from[,"peptide"])
      spectrumToPeptide <- spectrumToPeptide(template)[
                     names(spectrumToPeptide(template)) %in% from[,"spectrum"]]
      
      isoforms <-template@isoformToGeneProduct[names(indistinguishableProteins),]
      peptideInfo <- subset(template@peptideInfo, peptide %in% from[,'peptide'])
      ## TODO: overlappingProteins are missing

      return(
          new("ProteinGroup",
              peptideNProtein = peptideNProtein,
              peptideSpecificity = peptideSpecificity,
              proteinGroupTable = protein.group.table,
              indistinguishableProteins = indistinguishableProteins,
              spectrumToPeptide = spectrumToPeptide,
              isoformToGeneProduct = isoforms,
              proteinInfo = proteinInfo,
              peptideInfo = peptideInfo)
      )
    }
)

readProteinGroup <- function(id.file,...) {
  # TODO: check if this works
  pp <- do.call(rbind,lapply(id.file,function(id.f)
                read.table(id.f,header=T,stringsAsFactors=FALSE,sep="\t")[,c(.SPECTRUM.COLS[c('SPECTRUM','PEPTIDE')],
           .PEPTIDE.COLS['STARTPOS'],
           .PROTEIN.COLS['PROTEINAC'])]
))
#mzident <- lapply(id.file,read.table,header=T,stringsAsFactors=FALSE,sep="\t")
# si <- do.call(rbind,lapply(mzident,function(x) x$spectrum.identifications))
# rownames(si) <- NULL
# si <- si[order(si[,.SPECTRUM.COLS['SPECTRUM']]),]
  
#  err.spectra <- .get.dupl.n.warn(si[,.SPECTRUM.COLS[c('SPECTRUM','PEPTIDE','MODIFSTRING')]],.SPECTRUM.COLS['SPECTRUM'])
#  pp <- do.call(rbind,lapply(mzident,function(x) x$data))
#  pp <- pp[.SPECTRUM.COLS[c('PEPTIDE','SPECTRUM')],
#           .PEPTIDE.COLS['STARTPOS'],
#           !pp[,.SPECTRUM.COLS['SPECTRUM']] %in% err.spectra,c(.PROTEIN.COLS['PROTEINAC'])]
  return(ProteinGroup(from=pp,...))
}

setMethod("ProteinGroup",signature(from="data.frame",template="NULL",proteinInfo="ANY"),
          function(from,template,proteinInfo) ProteinGroup(from,proteinInfo=proteinInfo))

setMethod("ProteinGroup",signature(from="data.frame",template="missing",proteinInfo="ANY"),
          function(from,proteinInfo) {

      from <- .factor.as.character(from)

      if (ncol(from) == 3) {
          colnames(from) <- c("spectrum","peptide","protein")
          from$start.pos <- 0
      } else if (ncol(from) == 4) {
          colnames(from) <- c("spectrum","peptide","start.pos","protein")
      }

      spectrumToPeptide <- .as.vect(unique(from[,c("spectrum","peptide")]))
      peptideInfo <- unique(from[,c("protein","peptide","start.pos")])
      peptideInfo <- peptideInfo[order(peptideInfo[,"protein"],
                                       peptideInfo[,"start.pos"],
                                       peptideInfo[,"peptide"]),]
      from <- unique(from[,c("protein","peptide")])
      
      # use numbers to not exceed the maximum length
      from$peptn <- as.numeric(as.factor(from$peptide))
      from$protn <- as.numeric(as.factor(from$protein))
      
      # with which peptide-combination are proteins detected?
      combn.peptides <- ddply(from,"protein",
                              function(s){paste(sort(unique(s$peptn)),collapse="-")})
      colnames(combn.peptides) <- c("protein","peptides")
      
      # group proteins with same peptide combinations together
      # - those are indistiguishable
      combn.proteins=data.frame(); 
      for (c in unique(combn.peptides$peptides) ) { 
        combn.proteins=rbind(combn.proteins,
            data.frame(
                protein.g=paste(combn.peptides$protein[
                  combn.peptides$peptides %in% c],collapse=","),
                peptides=c)); 
      }
      combn.proteins$peptides <- as.character(combn.proteins$peptides)
      
      # merge data frames
      from <- merge(from,merge(combn.peptides,combn.proteins))
      
      pep.n.prots <- unique(from[,c("peptide","protein.g")])
      pep.n.prots$peptide <- as.character(pep.n.prots$peptide)
      pep.n.prots$protein.g <- as.character(pep.n.prots$protein.g)
      prots.to.prot <- as.matrix(unique(from[,c("protein.g","protein")]))
      prots.group <- c()
      
      
      # grouping:
      # 1) proteins w/ specific peptides [group reporterProteins]
      ssp <- table(as.matrix(pep.n.prots))==1
      sel.ssp <- pep.n.prots[,"peptide"] %in% names(ssp)[ssp]
      # take first those proteins which explain most specific peptides
      t <- table(pep.n.prots[sel.ssp,"protein.g"])
      protein.reporterProteins <- names(sort(t,decreasing=TRUE))
      
      # proteins to consider for grouping
      pep.n.prots.toconsider <-
        pep.n.prots[!pep.n.prots[,"protein.g"] %in% protein.reporterProteins,]
      
      protein.group.table <- matrix(nrow=0,ncol=5)
      colnames(protein.group.table) <- c("reporter.protein","protein.g","n.reporter-specific","n.group-specific","n.unspecific")
      
      # 2) get proteins, which are contained by those detected proteins
      #    (i.e. no specific and no overlapping peptides)
      for (reporter.protein in protein.reporterProteins) {
        protein.group.table <- rbind(protein.group.table,
               .build.protein.group(reporter.protein,pep.n.prots,pep.n.prots.toconsider))
        pep.n.prots.toconsider <- 
               pep.n.prots.toconsider[!pep.n.prots.toconsider[,"protein.g"] %in% protein.group.table[,"protein.g"],]

      }
      # proteins left in pep.n.prots.toconsider are those proteins
      # with overlapping peptides only
      # 3. useful peptides
      used.peptides <- pep.n.prots[pep.n.prots[,"protein.g"] %in% 
                                   protein.group.table[,"protein.g"],"peptide"]
      avail.pepnprots <- !pep.n.prots.toconsider[,"peptide"] %in% used.peptides
      while (any(avail.pepnprots)) {
        proteins.good <-  unique(pep.n.prots.toconsider[avail.pepnprots,"protein.g"])
        sel <- pep.n.prots.toconsider[,"protein.g"] %in% proteins.good
        # take protein which explains most peptides which are left
        t1 <- table(pep.n.prots.toconsider[sel&avail.pepnprots,"protein.g"])
        t2 <- table(pep.n.prots.toconsider[sel,"protein.g"])
        reporter.protein <- proteins.good[order(t2,t1,decreasing=TRUE)[1]]
        #reporter.protein <- proteins.good[1]

        protein.group.table <- rbind(protein.group.table,
               .build.protein.group(reporter.protein,pep.n.prots,pep.n.prots.toconsider))
        pep.n.prots.toconsider <- 
               pep.n.prots.toconsider[!pep.n.prots.toconsider[,"protein.g"] %in% protein.group.table[,"protein.g"],]
        used.peptides <- pep.n.prots[pep.n.prots[,"protein.g"] %in% 
                                   protein.group.table[,"protein.g"],"peptide"]
        avail.pepnprots <- !pep.n.prots.toconsider[,"peptide"] %in% used.peptides
      }
      
      protein.group.table <- as.data.frame(protein.group.table,stringsAsFactors=FALSE)
      protein.group.table$is.reporter <-
        protein.group.table$protein.g == protein.group.table$reporter.protein
      
      # define peptide specificity - pep.n.prots.toconsider are ignored here for now
      tmp <- merge(as.data.frame(protein.group.table,stringsAsFactors=FALSE),
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
      isoforms <- data.frame(proteinac.w.splicevariant = proteins,
                             proteinac.wo.splicevariant = proteins,
                             splicevariant = NA,stringsAsFactors=FALSE)
      pos.isoforms <- grep("^[^-]*-[0-9]*$",proteins)
      isoforms$proteinac.wo.splicevariant[pos.isoforms] <- 
         sub("-[^-]*$","",proteins[pos.isoforms])
      isoforms$splicevariant[pos.isoforms] <- sub(".*-","",proteins[pos.isoforms])
      rownames(isoforms) <- proteins
      
      return(
          new("ProteinGroup",
              peptideNProtein = as.matrix(pep.n.prots),
              peptideSpecificity = peptideSpecificity,
              proteinGroupTable = protein.group.table,
              indistinguishableProteins =
                .as.vect(prots.to.prot,col.data='protein.g',col.names='protein'),
              spectrumToPeptide = spectrumToPeptide,
              overlappingProteins = as.matrix(pep.n.prots.toconsider),
              isoformToGeneProduct = isoforms,
              peptideInfo = peptideInfo,
              proteinInfo = proteinInfo)
      )
      
    }
)

.build.protein.group <- function(reporter.protein,pep.n.prots,
                               pep.n.prots.toconsider) {

    protein.group.table <- matrix(c(reporter.protein,reporter.protein,0,0,0),ncol=5)
    colnames(protein.group.table) <- c("reporter.protein","protein.g","n.reporter-specific","n.group-specific","n.unspecific")
    pep.for.reporter <- pep.n.prots[pep.n.prots[,"protein.g"]==reporter.protein,"peptide"]
   
    groupmembers <- c()    
    pep.n.prots.toconsider <- pep.n.prots.toconsider[pep.n.prots.toconsider[,"protein.g"] != reporter.protein,]

    # proteins which have a subset of those proteins:
    proteins.w.subsetofpeptides <-
      pep.n.prots.toconsider[pep.n.prots.toconsider[,"peptide"] %in% pep.for.reporter,,drop=FALSE]

    # do these proteins have no common peptide w/ other proteins?
    for (subset.protein in unique(proteins.w.subsetofpeptides[,"protein.g"])) {
      #n peptides shared == n peptides total
      n.peptides.shared <- sum(proteins.w.subsetofpeptides[,"protein.g"] == subset.protein)
      n.peptides.total <- sum(pep.n.prots[,"protein.g"] == subset.protein)
      if (n.peptides.shared == n.peptides.total) {
        protein.group.table <- rbind(protein.group.table,c(reporter.protein,subset.protein,0,0,0))
        pep.n.prots.toconsider <- pep.n.prots.toconsider[!pep.n.prots.toconsider[,"protein.g"] %in% subset.protein,]
        groupmembers <- c(groupmembers,subset.protein)
      }
    }

    other.pep            <- pep.n.prots[pep.n.prots[,"protein.g"] != reporter.protein,"peptide"]
    reporterspecific.pep <- pep.for.reporter[!pep.for.reporter %in% other.pep]
    pep.for.groupmembers <- pep.n.prots[pep.n.prots[,"protein.g"] %in% groupmembers,"peptide"]
    groupspecific.pep    <- pep.for.groupmembers[pep.for.groupmembers %in% other.pep]
    unspecific.pep       <- pep.for.reporter[!pep.for.reporter %in% reporterspecific.pep & 
                                             !pep.for.reporter %in% groupspecific.pep]    

    for (i in seq_len(nrow(protein.group.table))) {
      pep.for.prot <- pep.n.prots[pep.n.prots[,"protein.g"]==protein.group.table[i,"protein.g"],"peptide"]
      protein.group.table[i,"n.reporter-specific"] <- sum(pep.for.prot %in% reporterspecific.pep)
      protein.group.table[i,"n.group-specific"] <- sum(pep.for.prot %in% groupspecific.pep)
      protein.group.table[i,"n.unspecific"] <- sum(pep.for.prot %in% unspecific.pep)
    }

    return(protein.group.table)
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


getProteinInfoFromUniprot <- function(x,splice.by=200) {
  fields <- c(accession="id",name="entry%20name",protein_name="protein%20names",
              gene_name="genes",organism="organism",
              length="length",sequence="sequence")
   
  protein.acs <- x@isoformToGeneProduct[names(indistinguishableProteins(x)),"proteinac.wo.splicevariant"]
  protein.info <- c()
  i <- 1
  while (i < length(protein.acs)) {
    url <- paste("http://www.uniprot.org/uniprot/?query=",
                 paste0("accession:",protein.acs[seq(from=i,to=min(length(protein.acs),i+splice.by-1))],collapse="+OR+"),
                 "&format=tab&compress=no&columns=",
                 paste0(fields,collapse=","))
    protein.info <- rbind(protein.info,read.delim(url,stringsAsFactors=FALSE,col.names=names(fields)))
    i <- i + splice.by
  }
  if (nrow(protein.info) > 0) {
    protein.info$protein_name <- sapply(strsplit(protein.info$protein_name," (",fixed=TRUE),function(x) x[1])
    protein.info$gene_name <- sapply(strsplit(protein.info$gene_name," "),function(x) x[1])
    protein.info$sequence <- gsub(" ","",protein.info$sequence)
  } else {
    warning("getProteinInfoFromUniprot returned no results for ",protein.acs," accessions")
  }
  return(protein.info)
}

getProteinInfoFromBioDb <- function(x,...,con=NULL) {
  if (is.null(con)) {
    con <- dbConnect(...)
    do.disconnect <- TRUE
  } else {
    do.disconnect <- FALSE
  }

  protein.acs <- x@isoformToGeneProduct[names(indistinguishableProteins(x)),"proteinac.wo.splicevariant"]
  query <- paste("SELECT primaryac AS accession,id AS name,",
                 "  description AS protein_name,",
                 "  (SELECT g.genename FROM genenames g WHERE g.entryid=d.entryid AND g.synonym=FALSE AND g.sourcedb=3 LIMIT 1) AS gene_name,",
                 "  os AS organism, seqlength as length, sequence",
                 "FROM dbentries d ",
                 "WHERE dbid IN (2,3) AND primaryac IN (",paste0("'",protein.acs,"'",collapse=","),")")
  res <- dbGetQuery(con,query)
  if (do.disconnect)
    dbDisconnect(con)
  return(res)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("data.frame","ProteinGroup",function(from) ProteinGroup(from))
setMethod("as.data.frame",signature(x="ProteinGroup"),
          function(x, row.names=NULL, optional=FALSE, ...) as(x,"data.frame"))

setAs("ProteinGroup","data.frame",function(from) {
  sp.df <- .vector.as.data.frame(spectrumToPeptide(from),
                                 colnames=c("spectrum","peptide"))
  ip.df <- .vector.as.data.frame(indistinguishableProteins(from),
                                 colnames=c("protein","protein.g"))
  
  ## merge spectrum to peptide to indistinguishable proteins (protein.g)
  spectrum.to.proteing <-
    merge(sp.df, as.data.frame(peptideNProtein(from),stringsAsFactors=FALSE))
 
  ## merge with indist. proteins
  spectrum.to.protein <- merge(ip.df,spectrum.to.proteing)
  spectrum.to.protein.w.peptideinfo <- merge(spectrum.to.protein,from@peptideInfo)

  return(spectrum.to.protein.w.peptideinfo[,c("spectrum","peptide","start.pos","protein")])
})

## Export a /concise/ data.frame which is peptide centric:
##  - for each peptide there is one line only
##  - column proteins concatenates all the different protein/protein group ACs
##    - indistinguishable peptides are separated by ',', protein groups by ';'
##  - column n.groups: number of groups a peptide appears
##  - column n.acs: number of acs a peptide appears in (without splice variant!)
##      counts double if the same ACs are in different groups 
##  - column n.variants: number of acs a peptide appears in (with splice variant!)
setAs("ProteinGroup","data.frame.concise",
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
      })

.paste_unique <- function(x,...,na.rm=TRUE) {
  x <- unique(x)
  x <- x[!is.na(x)]
  paste(x,...)
}

.proteinGroupAsConciseDataFrame <- function(from,only.reporters=TRUE,show.proteinInfo=TRUE,human.protein.acs=TRUE) {
        rp <- reporterProteins(from)
        p.ac <- .protein.acc(rp,proteinInfo(from),indistinguishableProteins(from))
        names(p.ac) <- rp
        ip.df <- .vector.as.data.frame(indistinguishableProteins(from),
                                       colnames=c("protein","protein.g"))
        if (only.reporters)
          ip.df <- ip.df[ip.df$protein.g %in% reporterProteins(from),]
        ip.df <- merge(ip.df,from@isoformToGeneProduct,
                       by.x="protein",by.y="proteinac.w.splicevariant")
        ip.df <- merge(ip.df,proteinGroupTable(from)[,c("protein.g","reporter.protein")])
        in.df <- ddply(ip.df, "reporter.protein",
                        function(x) 
                          data.frame(proteinn=paste(ddply(x,"proteinac.wo.splicevariant",
                            function(x) ifelse(any(!is.na(x$splicevariant)),sprintf("%s-[%s]",unique(x$proteinac.wo.splicevariant),
                                                                                    paste(sort(as.numeric(unique(x$splicevariant))),collapse=",")),
                                               x$proteinac.wo.splicevariant))$V1,collapse=",")
                        ))

        ip.df <- merge(ip.df,in.df)
        pep.n.prot <- as.data.frame(peptideNProtein(from),stringsAsFactors=FALSE)
        pep.n.prot <- merge(pep.n.prot,ip.df)

        res <- ddply(pep.n.prot,"peptide",
                     function(x) {
                       res <- data.frame(n.acs=length(unique(x[,"proteinac.wo.splicevariant"])),
                                         n.variants=length(unique(x[,"protein"])))
                       x <- unique(x[,c("reporter.protein","protein.g","proteinn")])
                       protein.gs <- unique(x[,'reporter.protein'])
                       res <- data.frame(proteins=paste(unique(x$proteinn),collapse=";"))
                       if (show.proteinInfo) 
                         res <- cbind(res,
                                      ID=.paste_unique(proteinInfo(from,protein.gs,do.warn=FALSE,collapse=","),collapse=","),
                                      Description=.paste_unique(proteinInfo(from,protein.gs,"protein_name",do.warn=FALSE,collapse=","),collapse=";"),
                                      Gene=.paste_unique(proteinInfo(from,protein.gs,"gene_name",do.warn=FALSE,collapse=","),collapse=","),stringsAsFactors=FALSE)
                       res <- cbind(res,n.groups=length(protein.gs),stringsAsFactors=FALSE)
                       if (!is.null(attr(from,"from.ids"))) 
                         res  <- cbind(groups=paste(attr(from,"from.ids")[protein.gs],collapse=","),
                                       res,stringsAsFactors=FALSE)
                       res

                     })
        return(unique(res))

}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setGeneric("peptideNProtein",function(x) standardGeneric("peptideNProtein"))
setGeneric("peptideSpecificity", function(x) standardGeneric("peptideSpecificity") )
setGeneric("peptideSpecificity",function(x) standardGeneric("peptideSpecificity"))
setGeneric("peptides",function(x,protein,...) standardGeneric("peptides"))
setGeneric("indistinguishableProteins",
           function(x,protein,protein.g) standardGeneric("indistinguishableProteins"))
setGeneric("reporterProteins",function(x, require.reporter.specific=FALSE) standardGeneric("reporterProteins"))
setGeneric("proteinGroupTable",function(x) standardGeneric("proteinGroupTable"))
setGeneric("spectrumToPeptide", function(x) standardGeneric("spectrumToPeptide"))
setGeneric("proteinInfo",function(x,protein.g,...) standardGeneric("proteinInfo"))
setGeneric("proteinInfo<-",function(x,value) standardGeneric("proteinInfo<-"))

setMethod("peptideSpecificity", "ProteinGroup", function(x) x@peptideSpecificity)
setMethod("peptideNProtein",   "ProteinGroup", function(x) x@peptideNProtein)
setMethod("proteinGroupTable",      "ProteinGroup", function(x) x@proteinGroupTable )

setMethod("indistinguishableProteins","ProteinGroup",
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
              warning(protein," is in no protein group")
              return(NA)
            }
            return(protein.groups)
          }
)

setMethod("spectrumToPeptide",   "ProteinGroup", function(x) x@spectrumToPeptide)

setMethod("peptides",signature(x="ProteinGroup",protein="missing"),
    function(x) peptideSpecificity(x)[,"peptide"]
)

setMethod("peptides",signature(x="ProteinGroup",protein="character"),
    function(x,protein,
             specificity=c(UNSPECIFIC,GROUPSPECIFIC,REPORTERSPECIFIC),
             columns=c('peptide'),set=union,drop=TRUE,
             groupspecific.if.same.ac=FALSE,do.warn=TRUE) {

      pnp <- peptideNProtein(x)
      ps <- peptideSpecificity(x)

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

      peptides <- ps[ps$specificity %in% specificity & 
                     ps$peptide %in% peptides,columns,drop=drop]
      if (length(peptides) == 0 && do.warn)
        warning("No peptide for protein ",protein," with specificity ",paste(specificity,collapse=","))
      return(peptides)
    }
)

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

setMethod("proteinInfo",signature(x="ProteinGroup",protein.g="missing"),function(x) x@proteinInfo)
setMethod("proteinInfo",signature(x="ProteinGroup",protein.g="character"),
    function(x,protein.g,select="name",collapse=", ",simplify=TRUE,do.warn=TRUE) {
      protein.info <- proteinInfo(x)
      sapply(protein.g,function(p) {
        if (is.null(proteinInfo(x)) || nrow(proteinInfo(x)) == 0) {
          if (do.warn)
            warning("protein info is NULL! Set for ProteinGroup.")
          return(NA)
        }
        if (!select %in% colnames(proteinInfo(x))) {
          if (do.warn)
            warning("column ",select," not available.\n",
                    "Available columns:\n","\t",paste(colnames(proteinInfo(x)),collapse="\n\t"))
          return(NA)
        }
        protein.acs <- x@isoformToGeneProduct[indistinguishableProteins(x,protein.g=p),"proteinac.wo.splicevariant"]
        sel <- protein.info$accession %in% protein.acs
        if (!any(sel)) {
          if (do.warn) warning("No protein info for ",p)
          return(NA)
        }
        protein.infos <- protein.info[sel,select]
        if (simplify) {
          paste(sort(unique(protein.infos)),collapse=", ")
        } else {
          names(protein.infos) <- protein.info$accession[sel]
          protein.infos
        }
      })
    }
)
setReplaceMethod("proteinInfo","ProteinGroup",
                 function(x,value) {
                   x@proteinInfo <- value
                   x
                 })

setGeneric("protein.g",function(x,pattern,variables=c("AC","name"),...) standardGeneric("protein.g"))
setMethod("protein.g",signature("ProteinGroup","character","ANY"),
          function(x,pattern,variables=c("AC","name"),...) {
  ip <- indistinguishableProteins(x)
  result <- c()
  for (p in pattern) {
    protein.acs <- c()
    if ("AC" %in% variables)
      result <- c(result,ip[grep(p,names(ip),...)])
    if ("name" %in% variables) {
      pi <- proteinInfo(x)
      if (length(pi) != 0L) {
        protein.acs <- unique(c(
           pi$accession[grep(p,pi$name,...)],
           pi$accession[grep(p,pi$gene_name,...)],
           pi$accession[grep(p,pi$protein_name,...)]
        ))
      }
      
      i.to.gp <- x@isoformToGeneProduct
      sel.isoforms <- i.to.gp[,"proteinac.wo.splicevariant"] %in% protein.acs
      protein.acs.w.isoforms <-
         i.to.gp[sel.isoforms,"proteinac.w.splicevariant"]
      protein.gs <- ip[names(ip) %in% protein.acs.w.isoforms]
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
      cat(sprintf("%8d protein groups with specific peptides\n",length(reporterProteins(object))))
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
   peptide.info <- x@peptideInfo
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
  df <- my.protein.info[,c("accession","splicevariant","gene_name","protein_name")]
  df$gene_name <- sanitize(df$gene_name)
  collapsed.splicevariant <- ddply(df,"accession",function(x) {
        only_one <- nrow(x) == 1
        x$splicevariant <- number.ranges(x$splicevariant)
        x <- unique(x)
        if (nrow(x) > 1) stop("something went wrong: ",x)
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
        if (is.na(x$gene_name) || x$gene_name == "") {
          x$name_nolink <- x$ac_nolink
        } else {
          x$ac_link_bold <- sprintf("\\textbf{%s}: %s",unique(x$gene_name),x$ac_link)
          x$ac_link <- sprintf("%s: %s",unique(x$gene_name),x$ac_link)
          x$ac_nolink <- sprintf("%s: %s",unique(x$gene_name),x$ac_nolink)
          x$name_nolink <- sprintf("%s: %s",unique(x$gene_name),unique(x$protein_name))
        }
        return(unique(x))
    })
    return(collapsed.gene_name)
}


my.protein.info <- function(x,protein.g) {
    protein.acs <- indistinguishableProteins(x,protein.g=protein.g)
    isoforms <- x@isoformToGeneProduct
    res <- data.frame(protein.ac=protein.acs,
                      accession=isoforms[protein.acs,"proteinac.wo.splicevariant"],
                      splicevariant=as.numeric(isoforms[protein.acs,"splicevariant"]), 
                      stringsAsFactors=FALSE)

    if (length(proteinInfo(x)) > 0) 
      res <- merge(res,proteinInfo(x),all.x=TRUE,all.y=FALSE)
    else 
      res <- cbind(res,gene_name=NA,protein_name=NA)


    return(res)
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

spectra.count <- function(protein.group,protein.g=reporterProteins(protein.group),
                          specificity=c("reporter-specific","group-specific","unspecific"),...) {
  peptide.spectra.count <- table(spectrumToPeptide(protein.group))
  ## Calculate unique spectrum counts for all proteins
  spectra.count <- sapply(protein.g, function(p)
                 sum(peptide.spectra.count[peptides(protein.group,protein=p,
                                              specificity=specificity,...)]))
  names(spectra.count) <- protein.g
  return(spectra.count)
}

sequence.coverage <- function(protein.group,protein.g=reporterProteins(protein.group),
                              specificity=c("reporter-specific","group-specific","unspecific"),
                              simplify=TRUE,...) {
  if ("length" %in% colnames(proteinInfo(protein.group)) && 
      "start.pos" %in% colnames(protein.group@peptideInfo)) {
    lengths <- proteinInfo(protein.group,protein.g=protein.g,"length",simplify=FALSE)
    peptides <- peptides(protein.group,protein=protein.g,specificity=specificity,...)

    sapply(protein.g, function(p) {
                     protein.acs <- protein.ac(protein.group,p)
                     seqcov <- sapply(protein.acs,function(pp) {
                            if (is.na(lengths[[p]][pp]))
                              return(NA)
                            seqq <- rep(FALSE,lengths[[p]][pp])
                            pi <- subset(protein.group@peptideInfo,protein==pp&peptide%in%peptides)
                            pi$peplength <- nchar(pi$peptide)
                            pi$end.pos <- pi$start.pos+pi$peplength-1
                            for (i_r in seq_len(nrow(pi)))
                              seqq[seq(from=pi[i_r,"start.pos"],to=pi[i_r,"end.pos"])]  <- TRUE
                            return(sum(seqq)/length(seqq))
                    })
                    if (simplify && length(seqcov) > 1)
                      mean(seqcov,na.rm=TRUE)
                    else
                      seqcov
      })    
  } else {
    stop("Necessary information for sequence coverage not available")
  }
}


calculate.dNSAF <- function(protein.group) {
  if (is.null(proteinInfo(protein.group)) || length(proteinInfo(protein.group)) == 0)
    stop("slot proteinInfo not set - it is needed for sequence length")
  if (!"length" %in% colnames(proteinInfo(protein.group))) {
    stop("no column 'length' in proteinInfo slot")
  }
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
    protein.length <- sum(seqlength[unique(protein.ac(protein.group,p))],na.rm=TRUE)
    dSAF <- (uspc[p] + d.sspc) / protein.length
    if (is.nan(dSAF) | !is.finite(dSAF)) dSAF <- NA
    return(dSAF)
  })

  ## normalize to dNSAF
  dNSAF <- dSAF/sum(dSAF,na.rm=TRUE)
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

n.observable.peptides <- function(seq,nmc=1,min.length=6,min.mass=800,max.mass=4000,...) {
  if (is.na(seq) || length(seq)==0 || nchar(seq) == 0)
    return(0)
  pep <- Digest(seq,missed=nmc,...)
  min.length.ok <- nchar(pep[,"peptide"]) >= min.length
  mass.ok <- pep[,"mz2"] >= min.mass & pep[,"mz3"] <= max.mass
  #pep[min.length.ok & mass.ok,]
  return(sum(min.length.ok & mass.ok))
}

calculate.emPAI <- function(protein.group,protein.g=reporterProteins(protein.group), ...) {
  require(OrgMassSpecR)
  if (is.null(proteinInfo(protein.group)) || length(proteinInfo(protein.group)) == 0)
    stop("slot proteinInfo not set - it is needed for sequence length")
  if (!"sequence" %in% colnames(proteinInfo(protein.group))) {
    stop("no column 'sequence' in proteinInfo slot")
  }

  sequences <- proteinInfo(protein.group)$sequence
  names(sequences) <- proteinInfo(protein.group)$accession
 
  n.observed.peptides <- table(peptideNProtein(protein.group)[,"protein.g"])[protein.g]
  n.observable.peptides <- sapply(protein.g,
       function(protein.g) 
         do.call(sum,lapply(unique(protein.ac(protein.group,protein.g)),
             function(protein.ac) n.observable.peptides(sequences[protein.ac])),...)
  )
  empai <- 10^(n.observed.peptides/n.observable.peptides) - 1
  return(empai)
}
