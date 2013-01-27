### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor and readIBSpectra.
###

.check.columns <- function(identifications,data.ions,data.mass,allow.missing.columns) {
  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  PC <- .PROTEIN.COLS[.PROTEIN.COLS %in% colnames(identifications)]
  missing.cols = c()
  for (col in c('SPECTRUM','PEPTIDE','MODIFSTRING'))
    if (!is.element(.SPECTRUM.COLS[col],SC))
      missing.cols <- c(missing.cols,.SPECTRUM.COLS[col])
  for (col in c('PROTEINAC'))
    if (!is.element(.PROTEIN.COLS[col],PC))
      missing.cols <- c(missing.cols,.PROTEIN.COLS[col])

  if (!is.null(data.ions) && !is.null(data.mass)) {
    if (is.null(rownames(data.ions)) || is.null(rownames(data.mass)))
      stop("provide spectrum ids as rownames for data.ions and data.mass")
    if (!identical(rownames(data.ions),rownames(data.mass)))
      stop("data.ions and data.mass do not have identical rownames!")
  }

  # handle missing columns
  if (length(missing.cols) > 0) {
      msg <- paste("not all required columns in identifications, the following are missing: \n\t",
               paste(missing.cols,collapse="\n\t"))
      if (allow.missing.columns) {
          warning(msg)
          identifications[,missing.cols] <- NA
      } else {
          stop(msg)
      }
  }
  return(identifications)
}

.get.quant <- function(identifications,colname,reporterTagNames) {
    cols <- sprintf(colname,reporterTagNames)
    if (!all(cols %in% colnames(identifications)))
      stop(" Quantitative information missing: Supply columns '",paste(cols,collapse="', '"),"'.")
    qdata <- unique(identifications[,c(.SPECTRUM.COLS['SPECTRUM'],cols)])
    identifications <<- identifications[,-which(colnames(identifications) %in% cols)]
    rownames(qdata) <- qdata[,.SPECTRUM.COLS['SPECTRUM']]
    qdata <- as.matrix(qdata[,-1])
    colnames(qdata) <- reporterTagNames
    return(qdata)
}

.check.assayDataElements <- function(assayDataElements) {

  for (elem in assayDataElements) {
    if (is.null(rownames(elem)))
      stop("provide rownames for all assayDataElements")
    #if (!identical(rownames(elem),rownames(data.ions)))
    #  stop("all assayDataElements must have the same rownames")
    #if (!identical(dim(data.ions),dim(elem)))
    #  stop("all assayDataElements must have the same dim")
  }
}

.get.quant.elems <- function(assayDataElements,data.ions,data.mass,spectra.ids,fragment.precision) {
  if (is.null(assayDataElements))
    assayDataElements <- list()

  assayDataElements$ions <- data.ions
  assayDataElements$mass <- data.mass
  
  if (!identical(spectra.ids,rownames(data.ions))) {
    id.n.quant <- intersect(spectra.ids,rownames(data.mass))
    id.not.quant <- setdiff(spectra.ids,rownames(data.mass))
    quant.not.id <- setdiff(rownames(data.mass),rownames(data))
    print(head(data.ions))
    print(head(spectra.ids))

    if (length(id.n.quant)==0) stop("No spectra could be matched between identification and quantifications")
    
    if (length(quant.not.id) > 0)
      warning(length(quant.not.id)," spectra could not be mapped from quant to id: ",
              paste(quant.not.id[1:min(length(quant.not.id),10)],collapse=";"))
    if (length(id.not.quant) > 0)
      warning(length(id.not.quant)," spectra could not be mapped from id to quant: ",
              paste(id.not.quant[1:min(length(id.not.quant),10)],collapse=";"))

    na.matrix <- matrix(NA,nrow=length(spectra.ids),ncol=ncol(data.mass),
                        dimnames=list(spectra.ids,colnames(data.mass)))

    for (elem in names(assayDataElements)) {
      tmp <- na.matrix
      tmp[id.n.quant,] <- assayDataElements[[elem]][id.n.quant,]
      assayDataElements[[elem]] <- tmp
    }
  }

  # filter based on fragment precision
  if (!is.null(fragment.precision)) {
    reporterMasses <- reporterTagMasses(.Object)
    min.masses <- reporterMasses - fragment.precision/2
    max.masses <- reporterMasses + fragment.precision/2
    for (i in seq_len(nrow(assayDataElements$mass))) {
      bad.mass <- assayDataElements$mass[i,]<min.masses |  assayDataElements$mass[i,] > max.masses
      if (any(bad.mass,na.rm=TRUE)) {
        for (elem in names(assayDataElements)) {
          assayDataElements[[elem]][i,bad.mass] <- NA
        }
      }
    }
  }

  assayDataElements$ions[which(assayDataElements$ions==0)] <- NA
  assayDataElements$mass[which(assayDataElements$mass==0)] <- NA
  if (all(apply(is.na(assayDataElements$mass),2,all))) stop("all masses are NA")
  if (all(apply(is.na(assayDataElements$ions),2,all))) stop("all intensities are NA")

  return(assayDataElements)
}

.get.varMetaData <- function(identifications) {
  nn <- .SPECTRUM.COLS %in% colnames(identifications)

  label.desc <- c(PEPTIDE='peptide sequence',
      MODIFSTRING='modifications of peptide',
      CHARGE='peptide charge state',
      THEOMASS='theoretical peptide mass',
      EXPMASS='experimental peptide mass',
	    PRECURSOR.ERROR='precursor error',
      PARENTINTENS='parent ion intensity',
      RT='retention time',
      DISSOCMETHOD='dissociation METHOD',
      SPECTRUM='spectrum title',
      SPECTRUM.QUANT='title of spectrum used for quantiation',
      PRECURSOR.PURITY="precursor purity",
      RAWFILE="raw file",NMC="nmc",
	    DELTASCORE="delta score",DELTASCORE.PEP="delta score peptide",
      SCANS="scans",
      SCANS.FROM="scans from",SCANS.TO="scans to",
	    MASSDELTA.ABS="massdelta (abs)",MASSDELTA.PPM="massdelta (ppm)",
      SEARCHENGINE='protein search engine',
      SCORE='protein search engine score',
	    SCORE.MSGF='MSGF+ score',

      SCORE.MASCOT="Mascot search score",
      SCORE.PHENYX="Phenyx search score",
      SCORE.PHOSPHORS='PhosphoRS pepscore',
      PROB.PHOSPHORS="PhosphoRS probability",

      PHOSPHO.SITES="phosphorylation sites",
      USEFORQUANT='use spectrum for quantification',
      SEQPOS='PTM seqpos',
      SITEPROBS='PhosphoRS site.probs',
      FILE='file',SAMPLE='sample',NOTES='notes'
      )
  if (!all(names(.SPECTRUM.COLS) %in% names(label.desc)))
    stop("Not all SPECTRUM COLS have a label description:\n\t",
         paste(names(.SPECTRUM.COLS)[!names(.SPECTRUM.COLS) %in% names(label.desc)],collapse="\n\t"))

  VARMETADATA=data.frame(labelDescription=label.desc[names(.SPECTRUM.COLS)],
                         row.names=.SPECTRUM.COLS)
	    
  return(VARMETADATA[nn,,drop=FALSE])
}

.merge.identifications.full <- function(identifications) {
  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  ## Substitute Isoleucins with Leucins (indistinguishable by Masspec)
  identifications[,.PEPTIDE.COLS['REALPEPTIDE']] <- identifications[,SC['PEPTIDE']]
  identifications[,SC['PEPTIDE']] <- gsub("I","L",identifications[,SC['PEPTIDE']])

  ## Separate protein columns (focus on peptide-spectrum matches)
  PC <- c(.PROTEIN.COLS['PROTEINAC'],.PEPTIDE.COLS['STARTPOS'],.PEPTIDE.COLS['REALPEPTIDE'])
  protein.colnames <- which(colnames(identifications) %in% c(SC['PEPTIDE'],PC))
  pept.n.prot <- unique(identifications[,protein.colnames])
  identifications <- unique(identifications[,-which(colnames(identifications) %in% PC)])

  ## Merge identifications
  if (max(table(identifications[,SC['SPECTRUM']])) > 1)
    identifications <- .merge.identifications(identifications)
  if (anyDuplicated(identifications[,SC['SPECTRUM']])) 
    stop('No divergent spectra identifications should be left.')
  
  ## Merge back the protein mapping
  merge(pept.n.prot,identifications,by="peptide")
}

setMethod("initialize","IBSpectra",
    function(.Object,identifications=NULL,data.ions=NULL,data.mass=NULL,
             proteinGroupTemplate=NULL,fragment.precision=NULL,
             assayDataElements=list(),allow.missing.columns=FALSE,
             write.excluded.to=NULL,...) { 

  if (is.null(identifications))
    return(callNextMethod(.Object,...))

  reporterTagNames <- reporterTagNames(.Object)
  identifications <- .factor.to.chr(identifications)

  ## Check that obligatory columns are present
  identifications <- .check.columns(identifications, data.ions, data.mass, allow.missing.columns)

  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  ## Substitute Isoleucins with Leucins (indistinguishable by Masspec)
  identifications[,.PEPTIDE.COLS['REALPEPTIDE']] <- identifications[,SC['PEPTIDE']]
  identifications[,SC['PEPTIDE']] <- gsub("I","L",identifications[,SC['PEPTIDE']])

  ## Separate protein columns (focus on peptide-spectrum matches)
  PC <- c(.PROTEIN.COLS['PROTEINAC'],.PEPTIDE.COLS['STARTPOS'],.PEPTIDE.COLS['REALPEPTIDE'])
  protein.colnames <- which(colnames(identifications) %in% c(SC['PEPTIDE'],PC))
  pept.n.prot <- unique(identifications[,protein.colnames])
  identifications <- unique(identifications[,-which(colnames(identifications) %in% PC)])

  ## Merge identifications
  if (max(table(identifications[,SC['SPECTRUM']])) > 1)
    identifications <- .merge.identifications(identifications)
  if (anyDuplicated(identifications[,SC['SPECTRUM']])) 
    stop('No divergent spectra identifications should be left.')
  rownames(identifications) <- identifications[,SC['SPECTRUM']]
  
  # Create ProteinGroup
  proteinGroup <- ProteinGroup(merge(pept.n.prot,identifications,by="peptide"),template=proteinGroupTemplate)
  
  ## Get intensities and masses in assayDataElements
  if (is.null(data.ions))
    data.ions <- .get.quant(identifications,.PEAKS.COLS['IONSFIELD'],reporterTagNames)
  if (is.null(data.mass))
    data.mass <- .get.quant(identifications,.PEAKS.COLS['MASSFIELD'],reporterTagNames)

  assayDataElements <- .get.quant.elems(assayDataElements,data.ions,data.mass,identifications[,SC['SPECTRUM']],fragment.precision)

  if (!.SPECTRUM.COLS['USEFORQUANT'] %in% colnames(identifications)) {
    # not perfect - better: set spectra of peptides shared between groups to FALSE
    #                            and spectra with no values
    identifications[,.SPECTRUM.COLS['USEFORQUANT']] <- TRUE
  }
    
  fdata <- identifications[,colnames(identifications) %in% .SPECTRUM.COLS]
  featureData <- new("AnnotatedDataFrame",data=fdata,
                     varMetadata=.get.varMetaData(fdata))
	   
  assayData=do.call(assayDataNew, assayDataElements, envir=parent.frame())
  
  ## create ibspectra object
  callNextMethod(.Object,
      assayData=assayData,
      featureData=featureData,
      proteinGroup=proteinGroup,
      reporterTagNames=reporterTagNames,...)
})



setGeneric("readIBSpectra", function(type,id.file,peaklist.file,...)
           standardGeneric("readIBSpectra"))

setMethod("readIBSpectra",
          signature(type="character",id.file="character",peaklist.file="missing"),
    function(type,id.file,...) {
      ll <- lapply(seq_along(id.file),function(i) {
                   message("\treading ",id.file[i])
                   df <- read.table(id.file[i],header=T,sep="\t")
                   df[,.SPECTRUM.COLS['FILE']]  <- id.file[i]
                   if (!is.null(names(id.file)))
                     df[,.SPECTRUM.COLS['SAMPLE']]  <- names(id.file)[i]
                   df
            })
      new(type,data=do.call(rbind,ll),...)
    }
)
setMethod("readIBSpectra",
          signature(type="character",id.file="data.frame",peaklist.file="missing"),
    function(type,id.file,...) new(type,data=id.file,...)
)

setMethod("readIBSpectra",
          signature(type="character",id.file="character",peaklist.file="character"),
    function(type,id.file,peaklist.file,
             mapping.file=NULL,mapping=c(peaklist="even",id="odd"),
             id.file.domap=NULL,id.format=NULL,decode.titles=TRUE,...) {
      
      data <- .read.identifications(id.file,mapping=mapping.file,mapping.names=mapping,
                                    identifications.quant=id.file.domap,
                                    identifications.format=id.format,decode.titles=decode.titles)

      readIBSpectra(type,data,peaklist.file,...)
})


##' readIBSpectra - read IBSpectra object from files
##'
##' <details>
##' @title 
##' @param type [character] IBSpectra type
##' @param id.file [character] file of format ibspectra.csv or
##' mzid. If format is ibspectra.csv, it may also contain quantitative
##' information, and peaklist.file is not required.
##' @param peaklist.file [character] file in MGF format containing
##' quantitative information (m/z and intensity pairs). If NULL,
##' quantitative information is expected to reside in id.file.
##' @param proteinGroupTemplate [ProteinGroup object] uses certain
##' protein grouping as template. Useful when comparing multiple
##' experiments.
##' @param mapping.file [character] CSV file which maps spectrum
##' titles from peaklist file to those of id file. Usefule when
##' quantitative and identification information reside in different
##' but corresponding spectra. E.g. HCD-CID dissociation
##' @param mapping [named character] column number or name for
##' 'peaklist' and 'id' spectra.
##' @param mapping.file.readopts 
##' @param peaklist.format 
##' @param id.format 
##' @return IBSpectra object of type 
##' @author Florian P Breitwieser
setMethod("readIBSpectra",
          signature(type="character",id.file="data.frame",peaklist.file="character"),
    function(type,id.file,peaklist.file,
             proteinGroupTemplate=NULL,annotate.spectra.f=NULL,
             peaklist.format=NULL,scan.lines=0,
             fragment.precision=NULL,fragment.outlier.prob=NULL,...) {
      
      log <- data.frame(key=c(),message=c())

      data <- id.file
      if (is.function(annotate.spectra.f)) {
        data <- annotate.spectra.f(data,peaklist.file)
      }

      # all identified spectrum titles
      id.spectra <- unique(data[,.SPECTRUM.COLS['SPECTRUM']])

      data.ions=c()
      data.mass=c()
      data.titles=c()
      for (peaklist.f in peaklist.file) {
        peaklist.format.f <- peaklist.format
        if (is.null(peaklist.format.f)) {

          if (grepl(".mgf$",peaklist.f,ignore.case=TRUE))peaklist.format.f <- "mgf"
          else if (grepl(".mcn$",peaklist.f,ignore.case=TRUE))peaklist.format.f <- "mcn"
          else if (grepl(".intensities.csv$",peaklist.f,ignore.case=TRUE))peaklist.format.f <- "csv"
          else
            stop(paste0("cannot parse file ",peaklist.f," - cannot deduce format (mgf or mcn)"))          

        }

        if (tolower(peaklist.format.f) == "mgf") {
          .Object <- new(type)
          reporterMasses <- .Object@reporterTagMasses
          reporterTagNames <- .Object@reporterTagNames

          intensities.f <- .read.mgf(peaklist.f,reporterMasses,reporterTagNames,
                                     fragment.precision=fragment.precision,
                                     prob=fragment.outlier.prob,scan.lines=scan.lines)
          if (nrow(intensities.f$ions) == 0) { stop("only NA data in ions/mass") }
          data.titles <- c(data.titles,intensities.f$spectrumtitles)
          data.ions <- rbind(data.ions,intensities.f$ions)
          data.mass <- rbind(data.mass,intensities.f$mass)
        } else if (tolower(peaklist.format.f) == "mcn") {
            if (type != "iTRAQ4plexSpectra")
              stop("mcn format (iTracker) only supports iTRAQ 4plex spectra!")
            mcn.i <- read.table("itraqdta/i-Tracker.mcn",sep=",",skip=2,
                                colClasses=c("character","character",
                  "numeric","numeric","numeric","numeric",rep("NULL",36)),
                col.names=c("spectrum","spectrum.t","114","115","116","117",rep("",36)),
                check.names=FALSE)
            data.titles <- c(data.titles,mcn.i$spectrum)
            data.ions <- c(data.ions,as.matrix(mcn.i[,3:6]))
            data.mass <- c(data.mass,
                matrix(rep(114:117,nrow(mcn.i)),nrow=nrow(mcn.i),byrow=TRUE))
        } else if (tolower(peaklist.format.f) == "csv") {
          message("reading ",peaklist.f)
          res <- read.delim(peaklist.f,stringsAsFactors=FALSE)
          message("done reading ",peaklist.f)
          
          ## FIX::
          #if (!is.null(id.spectra) && length(id.spectra) != 0) 
          #  res <- res[res[,"spectrum"] %in% id.spectra,]

          data.titles <- c(data.titles,res[,"spectrum"])
          data.ions <- rbind(data.ions,as.matrix(res[,.grep_columns(res,"ions$")]))
          data.mass <- rbind(data.mass,as.matrix(res[,.grep_columns(res,"mass$")]))

        }
      }
      if (.SPECTRUM.COLS['SPECTRUM.QUANT'] %in% colnames(data))
        data.titles <- .do.map(data.titles,unique(data[,.SPECTRUM.COLS[c('SPECTRUM','SPECTRUM.QUANT')]]))
      rownames(data.ions)  <- data.titles
      rownames(data.mass)  <- data.titles
      ## TODO: check that all identified spectra are present in intensities

      
      new(type,data=data,data.mass=data.mass,data.ions=data.ions,...)
    }
)

## TODO: return log
#.read.mapping <- function(mapping.file,readopts,mapping) {
#  ## STUB
#  if (is.null(mapping.file)) return(NULL)

#  if (!all(c("peaklist","id") %in% names(mapping))) {
#    stop("readIBSpectra/mapping must be a named vector with the",
#         " names peaklist and id.")
#  }
  ## read mapping file(s)
#  mapping.quant2id <- do.call(rbind,lapply(mapping.file,function(f) {
#                                           readopts$file <- f
#                                           do.call(read.table,readopts)
#                }))

#  cn <-  colnames(mapping.quant2id)
#  colnames(mapping.quant2id)[c(mapping['id'],mapping['peaklist'])] <- c('id','peaklist')
#  colnames(mapping.quant2id)[cn == mapping['peaklist']] <- 'peaklist'
#  log <- rbind(log,data.frame(rep("mapping file",length(mapping.file)),mapping.file,stringsAsFactors=FALSE))
#  if (!is.null(id.file.domap)) {
#    data.m <- .read.idfile(id.file.domap,id.format,decode.titles=decode.titles,log)
#    map.spectrum <- mapping.quant2id[,"id"]
#    names(map.spectrum) <- mapping.quant2id[,"peaklist"]
#    data.m[,"spectrum"] <- map.spectrum[data.m[,"spectrum"]]
#    data <- rbind(data,data.m)
#  }
#}

.do.map <- function(spectrumtitles,mapping.quant2id) {
  #mapped.spectra.pl <-
  #  mapping.quant2id[,2][ mapping.quant2id[,1] %in% id.spectra]

  #write(mapped.spectra.pl,file='id_spectra.csv')
  #.stopiflengthnotequal(id.spectra,mapped.spectra.pl,
  #                      "not all identified spectra could be matched!\n",
  #                      "\tx ... identified spectra\n",
  #                      "\ty ... mapped spectra\n")
 
  spectra.map <- .as.vect(mapping.quant2id,1,2)
  .stopiflengthnotequal(spectrumtitles,
                        spectra.map[spectrumtitles],
                        "not all spectra could be matched!\n",
                        "\tx ... spectra with quant info\n",
                        "\ty ... identified spectra\n")
  spectra.map[spectrumtitles]
}

.get.dupl.n.warn <- function(df,col,msg="ibspectra",write.to=NULL,f=warning) {
  dupl <- .all.duplicate.rows(df,col)
  if (!is.null(write.to))
    write.table(dupl,file=write.to,row.names=FALSE,sep="\t")
  dupl.msg <- paste(apply(dupl,1,paste,collapse="; "),collapse="\n\t")
  f(sprintf("%s> divergent identifications in %s spectra [%s ids]:\n\t%s",
			  msg,length(unique(dupl[,col])),nrow(dupl),dupl.msg))
  return(unique(dupl[,col]))
}


### READ MzID
read.mzid <- function(f) {
  library(XML)
  doc <- xmlInternalTreeParse(f)
  ns=c(x="http://psidev.info/psi/pi/mzIdentML/1.0")

  message("mapping")
  searchdatabase.mapping <- data.frame(
    ref=xpathSApply(doc,"/x:mzIdentML/x:DataCollection/x:Inputs/x:SearchDatabase",
      xmlGetAttr,name="id",namespaces=ns),
    name=xpathSApply(doc,"/x:mzIdentML/x:DataCollection/x:Inputs/x:SearchDatabase",
      xmlGetAttr,name="name",namespaces=ns),stringsAsFactors=FALSE
  )
   
  peptide.mapping <- data.frame(
      ref=xpathSApply(doc,"/x:mzIdentML/x:SequenceCollection/x:Peptide",
        xmlGetAttr,name="id",namespaces=ns),
      peptide.seq=xpathSApply(doc,
        "/x:mzIdentML/x:SequenceCollection/x:Peptide/x:peptideSequence",
        xmlValue,namespaces=ns),
      modification=xpathSApply(doc,"/x:mzIdentML/x:SequenceCollection/x:Peptide",
        fun=function(pep) {
          pep.length <- nchar(xpathSApply(pep,"x:peptideSequence",xmlValue,namespaces=ns))
          modif.df <- data.frame(
            location=as.numeric(xpathSApply(pep,"x:Modification",
              xmlGetAttr,name="location",namespaces=ns)),
            name=xpathSApply(pep,"x:Modification/x:cvParam[@cvRef='UNIMOD']",xmlGetAttr,name="name",namespaces=ns),
            stringsAsFactors=FALSE)
          modif.tmp <- rep("",pep.length+2)
          modif.tmp[modif.df$location+1] = modif.df$name
          return(paste(modif.tmp,collapse=":"))
        },namespaces=ns),stringsAsFactors=FALSE)
    
  protein.mapping <- data.frame(
      ref=xpathSApply(doc,"/x:mzIdentML/x:SequenceCollection/x:DBSequence",
        xmlGetAttr,name="id",namespaces=ns),
      accession=xpathSApply(doc,"/x:mzIdentML/x:SequenceCollection/x:DBSequence",
        xmlGetAttr,name="accession",namespaces=ns),
      sdb.ref=xpathSApply(doc,"/x:mzIdentML/x:SequenceCollection/x:DBSequence",
        xmlGetAttr,name="SearchDatabase_ref",namespaces=ns),
                                stringsAsFactors=FALSE)
  
  message("spectrum mapping")
  records.spectrumIdentifications <- 
    do.call(rbind,xpathSApply(doc,"/x:mzIdentML/x:DataCollection/x:AnalysisData/x:SpectrumIdentificationList/x:SpectrumIdentificationResult",
              namespaces=ns,
              fun=function(sir) {
    ## _SpectrumIdentificationResult_ #
    ## All identification from searching one spectrum
    spectrum.id <- xmlGetAttr(sir,"spectrumID")
    spectrum.title <- xpathSApply(sir,"x:cvParam[@name='spectrum title']",xmlGetAttr,name='value',namespaces=ns)

    ## get spectrum identifications which pass threshold
    #if (length(siis) > 1) stop(paste("more than one SpectrumIdentificationItems for spectrum ",spectrum.title))
    siis <- xpathApply(sir,"x:SpectrumIdentificationItem[@passThreshold='true' and @rank='1']",namespaces=ns,fun=function(sii) {
      ## _SpectrumIdentificationItem_ #
      ## An identification of a single peptide of a specturm.
      ## Only take the one which passes the threshold.

      attr.s <- xmlAttrs(sii)
      pep.mapping <- peptide.mapping[peptide.mapping[,"ref"] == attr.s['Peptide_ref'],]
 
    spectrum.df <-
      data.frame(
                 spectrum=spectrum.title,
                 search.engine = "Mascot",
                 score         = xpathSApply(sii,"x:cvParam[@name='mascot:score']",xmlGetAttr,name="value",namespaces=ns),
                 peptide.ref   = attr.s['Peptide_ref'],
                 peptide       = pep.mapping[,"peptide.seq"],
                 modif         = pep.mapping[,"modification"],
                 theo.mass     = attr.s['calculatedMassToCharge'],
                 exp.mass      = attr.s['experimentalMassToCharge'],
                 stringsAsFactors=FALSE)
#rm(attr.s)
    pe.df <- do.call(rbind,xpathApply(sii,"x:PeptideEvidence",namespaces=ns,
          fun=function(pe) {
     attr.e <- xmlAttrs(pe)
     data.frame(
                 peptide.ev.id=attr.e["id"],
                 peptide.ev.start=attr.e["start"],
                 peptide.ev.end=attr.e["end"],
                 peptide.ev.nmc=attr.e["missedCleavages"])
    }))
    if (nrow(spectrum.df)>1) stop("spectrum df nrow > 1 - might be a problem")
    spectrum.df <- cbind(pe.df,spectrum.df)
    return (spectrum.df)
  })
  }
  ))

  # TODO: check memory leaks

  message("protein mapping")
  records.proteinDetections <- do.call(rbind,xpathApply(doc,
    "/x:mzIdentML/x:DataCollection/x:AnalysisData/x:ProteinDetectionList/x:ProteinAmbiguityGroup",
    namespaces=ns,
    fun=function(pag) {
    ## apply on ProteinAmbiguityGroup
    do.call(rbind,xmlApply(pag,function(pdh) {
      ## apply on ProteinDetectionHypothesis
      data.frame(dbseq.ref=xmlGetAttr(pdh,"DBSequence_ref"),
        peptide.ev.id=
          as.character(getNodeSet(pdh,"x:PeptideHypothesis/@PeptideEvidence_Ref",namespaces=ns))
    )}))
  }))

  free(doc)
  
  merge(records.spectrumIdentifications,records.proteinDetections,by="peptide.ev.id")
}

## read.mgf
.parse.spectrum <- function(input,reporterMasses,fragment.precision,recordNo) { 
  firstChar = substr(input,1,1);
  sel.numeric = firstChar == 1 & substr(input,4,4) == ".";
  sel.alpha = firstChar %in% c("P","T")
      
  d <- unlist(strsplit(input[sel.alpha],"="))
  d1 <- d[seq(2,length(d),by=2)]
  names(d1) <-d[seq(1,length(d),by=2)]
      
  spectrumTitles[recordNo,] <<- d1[c("TITLE","PEPMASS")]
# TODO: do not take closest but most intense peak?
      
  if (any(sel.numeric)) {
    mzi <- do.call("rbind",strsplit(input[sel.numeric],"\\s"))[,1:2]
    if (!is.null(nrow(mzi))) {
      mzi.mass <- as.numeric(mzi[,1])
      mzi.ions <- as.numeric(mzi[,2])
      if (length(mzi.mass)>0) {
        lapply(seq_along(reporterMasses),
               function(i) {
                 m <- abs(mzi.mass-reporterMasses[i])
                 pos <- which(m == min(m))
                 if (length(pos) > 0 & m[pos] < fragment.precision/2) {
                   observedMasses[recordNo,i] <<- mzi.mass[pos]
                   observedIntensities[recordNo,i] <<- mzi.ions[pos]
                 }
               }
               )
      }
    }
  }
}

##' .read.mgf: read isobaric reporter tag masses and intensities from MGFs
##'
##' MGF files list m/z and intensities for each spectrum. .read.mgf
##' extracts m/z and intensity pairs of masses corresponding to isobaric
##' tag masses (within a fragment precision).
##' @title 
##' @param filename [character] file name of one or multiple files in MGF format
##' @param type [character] denoting IBSpectra class - used to get
##'             reporter masses and reporter names.
##' @param spectra [character vector] if defined, only export spectra whose TITLE
##'                is listed. Speeds up the function when only identified spectra
##'                are of interest.
##' @param fragment.precision [numeric] take m/z-intensity pairs whose m/z value are
##'                          in a range of +/- fragment.precision/2 of 'true' mass.
##' @param prob [numeric] Filter out m/z-intensitiy values with the prob/2 most
##'             unprecise m/z values on both sides.
##' @param substitute.dta [boolean] internal. replace TITLEs: s/.dta.[0-9]*$/.dta/
##' @return list(ions, mass, spectrumtitles)
##' @author Florian P Breitwieser
.read.mgf <- function(filename,reporterMasses,reporterNames,spectra=NULL,fragment.precision=0.05,
                      prob=NULL,substitute.dta=FALSE,check.id.ok=FALSE,
                      scan.lines=0) {
  if (is.null(fragment.precision)) { fragment.precision=0.05 }
  message("  reading mgf file ",filename,
          " [fragment precision: ",fragment.precision,"]")

  ## get reporter masses and names from type
  nReporter <- length(reporterMasses)
  min.mass <- min(reporterMasses)-fragment.precision/2
  max.mass <- max(reporterMasses)+fragment.precision/2
  
  ## read (and concatenate if multiple) mgf files
  input <- c()
  for (f in filename) {
    con <- file(f,'r')
    if (scan.lines > 0) {
      while(length(f.input <- readLines(con, n=scan.lines)) > 0){
        input <- c(input,grep("^[A-Z]|^1[12][0-9]\\.",f.input,value=T))
      }
    } else {
      input <- c(input,readLines(con))
    }
    close(con)
  }
  if (substitute.dta)
    input <- sub(".dta.*$",".dta",input)
  
  ## index mgf file content: positions of BEGIN and END IONS
  begin_ions <- which(input=="BEGIN IONS")
  end_ions <- which(input=="END IONS")
  if (length(begin_ions) != length(end_ions))
    stop("mgf file is errorneous, non-matching number",
         " of BEGIN IONS and END IONS tags");

  bnd <- data.frame(begin=begin_ions,
                    end=end_ions)
  
  if (!is.null(spectra)) {
    ## filtering to take only spectra defined in spectra
    ## all.titles <- .trim(sapply(strsplit(grep("TITLE",input,value=T),"="),function(x) x[2] )) # not efficient
    all.titles <- .trim(substring(grep("^TITLE=",input,value=TRUE),7))
    if (length(all.titles) != nrow(bnd)) {
      # sanity check that each spectrum has a title
      stop("title not specified for all spectra!");
    }

    spectra.to.take <- all.titles %in% spectra

    if (sum(spectra.to.take) ==0) {
      stop("No identified spectrum is found in MGF file.\n",
           "  TITLEs in MGF [1:",length(all.titles),"]: \n\t",
           paste(all.titles[1:2],collapse=", "),", ...\n",
           "  identified spectra [1:",length(spectra),"]: \n\t",
           paste(spectra[1:2],collapse=", "),", ...\n")
    }
    bnd <- bnd[spectra.to.take,]
  }
  bnd$recordNo = seq_len(nrow(bnd))
  nSpectra <- nrow(bnd)
  message(nSpectra," spectra")

  ## create list with all spectra (header+mass list) as entries
  all.spectra <- apply(bnd,1,function(x) input[x[1]:x[2]])
  rm(input)

  ## if all spectra are of equal length, apply returns a matrix
  ##  convert to list
  if (is.matrix(all.spectra)) 
    all.spectra <- split(all.spectra,col(all.spectra))
  
  ## extract information from each spectrum
  result <- do.call(rbind,lapply(all.spectra,function(x) {
    header <- .strsplit_vector(x[grep("^[A-Z]",x)],"=")
    numbers <- do.call(rbind,strsplit(x[grep("^1..\\.",x)],"\\s"))
    mzi.mass <- as.numeric(numbers[,1])
    
    rr <- c(title=header["TITLE"])
    
    sel <- mzi.mass > min.mass & mzi.mass < max.mass
    if (any(sel)) {
      mzi.mass <- mzi.mass[sel]
      mzi.ions <- as.numeric(numbers[sel,2])

      rr <- c(rr,do.call(c,lapply(reporterMasses,function(y) {
        m <- abs(y-mzi.mass)
        pos <- which(m == min(m))
      
        if (length(pos) > 0 & m[pos] < fragment.precision/2)
          return(c(mzi.mass[pos],mzi.ions[pos]))
        else
          return(c(NA,NA))
      }))) 
    } else {
      rr <- c(rr,rep(NA,nReporter*2))
    }
    return(rr)
  }))
  rm(all.spectra)
  if (length(result) == 0 || nrow(result) == 0) {
    stop("error reading MGF file - could not parse masses and intensities.\n",
         "Check the mapping id <-> peaklist.")
  }
 
  mass <- apply(result[,seq(from=2,to=ncol(result),by=2)],2,as.numeric)
  ions <- apply(result[,seq(from=3,to=ncol(result),by=2)],2,as.numeric)
     
  ## only select spectra with itraq masses detected
  sel = apply(mass,1,function(x) any(!is.na(x)))
      
  if (!any(sel)) {
    stop("No values could be extracted from MGF ",filename,"\n",
         "with fragment precision ",fragment.precision,".\n")
  }

  if (!is.null(prob) && prob > 0) {
    ## boundaries: remove extreme outliers of mass spectrum
    ##  to test: another quantile.type might be more appropriate for small datasets
    bnd <- apply(mass[sel,],2,quantile,probs=c(prob*0.5,1-prob*0.5),na.rm=TRUE)
    sel.prob <- !is.na(mass) & (mass < matrix(bnd[1,],byrow=T,nrow=nrow(mass),ncol=ncol(mass)) |
                                mass > matrix(bnd[2,],byrow=T,nrow=nrow(mass),ncol=ncol(mass)))

    ions[sel.prob] <- NA
    mass[sel.prob] <- NA
    message("mass boundaries:\n\t",
            paste(colnames(mass),sprintf("%.5f : %.5f",bnd[1,],bnd[2,]),sep="\t",collapse="\n\t"))
  }

  ions <- ions[sel,,drop=FALSE]
  mass <- mass[sel,,drop=FALSE]
 
  spectrumtitles <- .trim(result[sel,1])
  dimnames(ions) <- list(spectrumtitles,reporterNames)
  dimnames(mass) <- list(spectrumtitles,reporterNames)
  rm(result)
  
  return(list(ions=ions, mass=mass,
              spectrumtitles=spectrumtitles))
}


## TODO: log is not returned
.read.idfile <- function(id.file,id.format=NULL,decode.titles=TRUE,log=NULL,...) {
  data <- unique(do.call("rbind",lapply(id.file,function(f) {
    
    if (is.null(id.format)) {
      if (grepl(".mzid$",f,ignore.case=TRUE)) 
        id.format.f <- "mzid"
      else if (grepl(".ibspectra.csv$",f,ignore.case=TRUE) ||
               grepl(".id.csv$",f,ignore.case=TRUE) || 
               grepl(".mascot.csv$",f,ignore.case=TRUE) || 
               grepl(".phenyx.csv$",f,ignore.case=TRUE))
        id.format.f <- "ibspectra.csv"
      else if (grepl(".peptides.csv$",f) || grepl(".peptides.txt$",f)) 
        id.format.f <- "rockerbox"
      else if (grepl(".msgfp.csv$",f) || grepl(".tsv$",f)) 
        id.format.f <- "msgfp tsv"
      else
        stop(paste("cannot parse file ",f," - cannot deduce format based on extenstion (it is not ibspectra.csv, id.csv, peptides.txt or mzid). Please provide id.format to readIBSpectra",sep=""))
    } else {
      id.format.f <- id.format
    }
    
    if (id.format.f == "ibspectra.csv") {
      data <- read.table(f,header=T,stringsAsFactors=F,sep="\t",...)
      log <- rbind(log,c("identification file [id.csv]",f))
    } else if (id.format.f == "mzid") {
      data <- read.mzid(f)
      log <- rbind(log,c("identification file [mzid]",f))
    } else if (id.format.f == "rockerbox") {
      data <- .read.rockerbox(f,...)
    } else if (id.format.f == "msgfp tsv") {
      data <- .read.msgfp.tsv(f,...)
    } else {
      stop(paste0("cannot parse file ",f," - format [",id.format.f,"] not known."))
    }
    return(data)
  })
  ))
  if (decode.titles) {
    data[,.SPECTRUM.COLS['SPECTRUM']] <- unlist(lapply(data[,.SPECTRUM.COLS['SPECTRUM']],URLdecode))
  }
  data[,.SPECTRUM.COLS['SPECTRUM']] <- .trim(data[,.SPECTRUM.COLS['SPECTRUM']])
  return(data)
}

.read.rockerbox <- function(f) {
  data.r <- read.table(f,header=T,stringsAsFactors=F,sep="\t")
  if (!"scan.title" %in% colnames(data.r))
    stop("no scan.title column in ",f,"; please use a Rockerbox version >= 2.0.6")

  ## transform modification (TODO)
  data.r$modif <- data.r$modifications
  ## end transform
  
  ## transform 'all.peptide.matches' to ac and start.pos
  split.acs <- strsplit(data.r$all.protein.matches,"; ")
  names(split.acs) <- data.r$all.protein.matches
  ac.n.startpos <- ldply(split.acs,function(x) {
    y <- do.call(rbind,strsplit(x,split="\\[|\\]"))
    start.pos <- sapply(strsplit(y[,2],"-",fixed=TRUE),
                        function(z) if(all(!is.na(as.numeric(z))) && length(z) == 2) as.numeric(z[1])
                        else stop("not numeric"))
    data.frame(accession=y[,1],start.pos=start.pos)  
  })
  data.r <- merge(data.r,ac.n.startpos,by.x="all.protein.matches",by.y=".id")
  data.r$all.protein.matches <- NULL
  ## end transform

  sel <- names(.ROCKERBOX.COLS) %in% names(c(.SPECTRUM.COLS,.PEPTIDE.COLS))
  data <- data.r[,.ROCKERBOX.COLS[sel]]
  colnames(data) <- c(.SPECTRUM.COLS,.PEPTIDE.COLS)[names(.ROCKERBOX.COLS)[sel]]
  return(data)
}

###############################################################################
## MERGE IDENTIFCATIONS


.dissect.search.engines <- function(identifications) {
  ## scores are merged together
  if (.SPECTRUM.COLS['SCORE'] %in% colnames(identifications)) {
    engines <- strsplit(identifications[,SC['SEARCHENGINE']],"|",fixed=TRUE)
    scores <- strsplit(identifications[,SC['SCORE']],"|",fixed=TRUE)
  } else {
    engine.n.score <- strsplit(identifications[,SC['SEARCHENGINE']],"[ \\|]")
    engines <- lapply(engine.n.score,function(x) x[seq(from=1,to=length(x),by=2)])
    if (length(unique(unlist(engines))>10)) {
      engines <- lapply(engine.n.score,function(x) x[seq(from=1,to=length(x)/2)])
      scores <- lapply(engine.n.score,function(x) as.numeric(x[seq(from=length(x)/2+1,to=length(x))]))
    } else {
      scores <- lapply(engine.n.score,function(x) as.numeric(x[seq(from=2,to=length(x),by=2)]))
    }
  }
  for (engine in unique(unlist(engines))) {
    name <- paste0('score.',tolower(engine))
    e.scores <- mapply(function(e,s) if(any(e==engine)) s[e==engine] else NA,engines,scores)
    identifications[,name] <- as.numeric(e.scores)
  }
  identifications[,SC['SEARCHENGINE']] <- sapply(engines,paste,collapse="&")
  if ('SCORE' %in% names(SC))
    identifications[,SC['SCORE']] <- NULL
  return(identifications)
}

.merge.search.engine.identifications <- function(identifications) {

  cols.for.merging <- intersect(colnames(identifications),setdiff(.SPECTRUM.COLS,.ID.COLS))

  ## Rounding is necessary. Mascot and Phenyx theoretical masses can differ a lot.
  if ('theo.mass' %in% colnames(identifications))
    identifications[,'theo.mass'] <- round(identifications[,'theo.mass'],0)
  if ('retention.time' %in% colnames(identifications))
    identifications[,'retention.time'] <- round(identifications[,'retention.time'],4)
  if ('parent.intens' %in% colnames(identifications))
    identifications[,'parent.intens'] <- round(identifications[,'parent.intens'],4)
  #for (cc in intersect(c('theo.mass','retention.time','parent.intens'),cols.for.merging)) {
  #  if (is.numeric(identifications[,cc]))
  #    identifications[,cc] <- round(identifications[,cc],0)
  #}

  cn.search.engine <- which(colnames(identifications) == .SPECTRUM.COLS['SEARCHENGINE'])
  ids.split <- lapply(split(identifications,identifications[,cn.search.engine]),
		      function(x) x[,-cn.search.engine])
  if (length(ids.split) < 2) stop("Expecting more than two different search engines.")
  if (length(ids.split) > 10) stop("Split data.frame into more than 10 search.engines - seems unrealistic. ",
				   "search.engine values: ",paste(names(ids.split),collapse=", "))

  clean.names <- gsub("[^[:alnum:]]","",tolower(names(ids.split))) ## lower case alpha-numeric names
  score.cols <- paste0("score.",clean.names)

  for (ii in seq_along(ids.split)) { # fix colnames which are duplicate
    cn <- colnames(ids.split[[ii]])
    cn[!cn %in% cols.for.merging] <- paste(cn[!cn %in% cols.for.merging],clean.names[ii],sep=".")
    colnames(ids.split[[ii]]) <- cn
  }

  ## Merge identification data frames
  ids.merged <- merge(ids.split[[1]],ids.split[[2]],by=cols.for.merging,all=TRUE)
  if (length(ids.split) > 2) {
    for (ii in seq(from=3,to=length(ids.split))) {
      ids.merged <- merge(ids.merged,ids.split[[ii]],by=cols.for.merging,all=TRUE)
    }
  }

  ## Handle spectra with different identifications with different search engines
  tt <- table(ids.merged[,'spectrum'])
  spectra.ok <- ids.merged[,'spectrum'] %in% names(tt)[tt==1]
  resolved.ids <- .resolve.differing.identifications(ids.merged[!spectra.ok,],score.cols)

  identifications <- rbind(ids.merged[spectra.ok,],resolved.ids)
  if (any(table(identifications[,'spectrum'])>1))
    stop("Merging not successful, duplicated psms!")

  identifications[,.SPECTRUM.COLS['SEARCHENGINE']] <- 
	  apply(identifications[,score.cols],1,function(x) { paste(names(ids.split)[!is.na(x)],collapse="&") })
  return(identifications)
}

.resolve.differing.identifications <- function(identifications,score.cols) {
  max.uniq <- function(tt) sum(tt==max(tt)) == 1
  take.max <- function(tt) names(tt)[which.max(tt)]
  allequal <- function(x) all(x == x[1])
  by.y <- function(x,ind,fun=sum) sapply(by(x,ind,fun),function(x) return(x) )
  skipna <- function(x) unlist(x[!is.na(x)])
  n.skipped <- 0
  n.max.ids <- 0
  n.modif.pos.dif <- 0

  resolved.identifications <- ddply(identifications,'spectrum', function(x) {
	n.ids <- rowSums(!is.na(x[,score.cols]),na.rm=TRUE)    ## number of identifications for each psm
	if (max.uniq(n.ids)) {
	  n.max.ids <<- n.max.ids + 1
	  return(x[which.max(n.ids),,drop=FALSE])
	}
        
	## resolve peptide differnces
        if (!allequal(x[,'peptide'])) {                ## different peptides identified
	  pep.ids <- by.y(n.ids,x[,'peptide'],sum)
	  if (max.uniq(pep.ids)) {                     ##   any peptide has been seen more often than others?
	    x <- x[x[,'peptide'] == take.max(pep.ids),,drop=FALSE]
	  }  else {
	    n.skipped <<- n.skipped + 1
	    return(NULL)
	  }
	}
	
	## all equal peptides, different modification positions
        if (!allequal(x[,'peptide'])) {
	  message("ERROR while merging: psm peptides should be equal, but they aren't.")
	  print(x)
	  stop()
	}

	## resolve modification differences
        if (!allequal(x[,'modif'])) {                ## different modifs identified
	  modif.ids <- by.y(n.ids,x[,'modif'],sum)
	  if (!max.uniq(modif.ids))
	    n.modif.pos.dif <<- n.modif.pos.dif + 1
	  x <- x[x[,'modif'] == take.max(modif.ids)[1],,drop=FALSE]
	} else if (!allequal(x[,'charge'])) {         ## different charge
          x[1,score.cols] <- as.numeric(sapply(x[,score.cols],skipna))
	  x[1,'charge'] <- max(x[,'charge'])
	  x <- x[1,,drop=FALSE]
	} else { 
          message("Why is the modif equal? ")
	  print(rbind(x,is.equal=apply(x,2,allequal)))
	  stop()
	}
	
	return(x)
  })
  message(" resolving differing ",length(unique(identifications[,'spectrum']))," identifications: " )
  message("   ",n.max.ids," resolved because more evidence was available for one alternative")
  message("   ",n.skipped," removed because of differing peptide ids")
  message("   ",n.modif.pos.dif," kept, as only the modification position was different")
  return(resolved.identifications)
}

.merge.quant.identifications <- function(identifications) {

  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  score.colname <- .SPECTRUM.COLS[c('PROB.PHOSPHORS','SCORE.MASCOT','SCORE.PHENYX','SCORE.MSGF')]
  score.colname <- score.colname[score.colname %in% colnames(identifications)]

  message("Merging identifications from different dissociation methods.")
  ddply(identifications,.SPECTRUM.COLS['SPECTRUM.QUANT'],function(x) {
                           if (nrow(x) == 1) return(x)
                           my.args <- as.list(x[,score.colname,drop=FALSE])
                           my.args <- lapply(my.args,round,digits=2) ## take two significant digits before taking next score into account as top hitter
                           my.args$decreasing=TRUE
                           max.hit <- do.call(order,my.args)[1]
                           if (!all(x[,SC['PEPTIDE']] == x[1,SC['PEPTIDE']]))
                             return(NULL)

                           x[max.hit,SC['DISSOCMETHOD']] = paste0("[",x[max.hit,SC['DISSOCMETHOD']],"]")
                           if (all(x[,SC['MODIFSTRING']] == x[max.hit,SC['MODIFSTRING']]) && 'DISSOCMETHOD' %in% names(SC)) {
                             if ("cid" %in% x[,SC['DISSOCMETHOD']])
                               x[max.hit,SC['SPECTRUM']] <- x[x[,SC['DISSOCMETHOD']]=="cid",SC['SPECTRUM']][1]

                             x[max.hit,SC['DISSOCMETHOD']] <- paste(x[,SC['DISSOCMETHOD']],collapse="&")
                             for (sc in score.colname[c(2,3)])
                               x[max.hit,sc] <- paste(x[,sc],collapse="&")
                           }

                           return(x[max.hit,])
  })

}

.merge.identifications <- function(identifications) {
  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]

  ## Merge results of different search engines / on different spectra
  if ('SEARCHENGINE' %in% names(SC) &&                                          ## search engine column is present
      length(unique(identifications[,SC['SEARCHENGINE']])) > 1 &&               ## there are multiple search engines defines
      !any(grepl('score.',colnames(identifications))))                          ## the individual score.ENGINENAME are not present
  {                     
    if (any(grepl("|",identifications[,SC['SEARCHENGINE']],fixed=TRUE))) {
      identifications <- .dissect.search.engines(identifications)               ## dissect merged id columns
    } else {
      identifications <- .merge.search.engine.identifications(identifications)  ## merge identifications
    }
  }

  ## Merge identifications from different dissociation methods (e.g. CID and HCD)
  ##   TODO: Can this be done by the code above?
  if ('SPECTRUM.QUANT' %in% names(SC)) {
    tt <- table(identifications[,SC['SPECTRUM.QUANT']])
    if (any(tt>1)) {
      spectra.ok <- identifcations[,'SPECTRUM.QUANT'] %in% names(tt)[tt==1]
      ids.quant.merged <- .merge.quant.identifications(identification[!spectra.ok])
      identifications <- rbind(identifications[spectra.ok,],ids.quant.merged)
    }
  }

  if (max(table(identifications[,SC['SPECTRUM']])) > 1)
    stop("Identifications still have duplicate spectra!")

  if ('SEARCHENGINE' %in% names(SC)) {
    message("Identification details:")
    tt <- table(identifications[,SC['SEARCHENGINE']])
    print(data.frame(perc=sprintf("%.2f %%",tt[order(tt)]/sum(tt)*100),n=sort(tt)))
  }

  return(identifications)
}

## end MERGE IDENTIFICATIONS
###############################################################################


.read.identifications <- function(identifications,
                                  mapping=NULL,mapping.names=c(quantification.spectrum="hcd",identification.spectrum="cid"),
                                  identifications.quant=NULL,decode.titles=TRUE,identifications.format=NULL) {

  ## load identifications (is either character or data.frame
  if (is.character(identifications) && file.exists(identifications))
    identifications <- .read.idfile(identifications,identifications.format,decode.titles)
  if (!is.data.frame(identifications)) stop("identifications must be a data.frame or valid file name")

  ## load mapping (either character or data.frame)
  if (!is.null(mapping)) {
    if (is.null(mapping.names)) stop("mapping.names must be given alongside with mapping")
    if (!'identification.spectrum' %in% names(mapping.names) || !'quantification.spectrum' %in% names(mapping.names)) 
      stop("mapping.names must be given alongside with mapping. names(mapping.names) must be 'identification.spectrum' and 'quantification.spectrum'")
    if (is.character(mapping) && file.exists(mapping))
      mapping <- do.call(rbind,lapply(mapping,read.table,sep=",",header=TRUE,stringsAsFactors=FALSE))
    if (!is.data.frame(mapping)) stop("mapping must be a data.frame or valid file name")
    if (!all(mapping.names %in% colnames(mapping))) stop("mapping.names must be column names of mapping")

    quant2id <- .as.vect(mapping,mapping.names['identification.spectrum'],mapping.names['quantification.spectrum'])
    id2quant <- .as.vect(mapping,mapping.names['quantification.spectrum'],mapping.names['identification.spectrum'])

    identifications[,.SPECTRUM.COLS['SPECTRUM.QUANT']] <- id2quant[identifications[,.SPECTRUM.COLS['SPECTRUM']]]
    identifications[,.SPECTRUM.COLS['DISSOCMETHOD']] <- ifelse(is.null(identifications.quant),"hcd","cid") # temp fix
    #identifications[,.SPECTRUM.COLS['DISSOCMETHOD']] <- names(mapping.names)[mapping.names=='identification.spectrum']

    if (!is.null(identifications.quant)) {
      ## load identifications.quant (either character or data.frame)
      if (is.character(identifications.quant) && file.exists(identifications.quant)) 
        identifications.quant <- .read.idfile(identifications.quant,identifications.format,decode.titles)
      if (!is.data.frame(identifications.quant)) stop("identifications.quant must be a data.frame or valid file name")
      identifications.quant[,.SPECTRUM.COLS['SPECTRUM.QUANT']] <- identifications.quant[,.SPECTRUM.COLS['SPECTRUM']]
      identifications.quant[,.SPECTRUM.COLS['DISSOCMETHOD']] <- "hcd"
      #identifications.quant[,.SPECTRUM.COLS['DISSOCMETHOD']] <- names(mapping.names)[mapping.names=='quantification.spectrum']
      identifications <- rbind(identifications,identifications.quant)
    }
  }
  return(identifications)
}


## read MSGF+ tab-separated identification files
.read.msgfp.tsv <- function(filename,filter.rev.hits=TRUE) {
  data <- read.delim(filename, sep="\t", stringsAsFactors=FALSE)
  ib.df <- data.frame(spectrum=data[,'Title'],.convert.msgfp.pepmodif(data[,'Peptide']),
		      scan.from=data[,'ScanNum'],dissoc.method=tolower(data[,'FragMethod']),
		      precursor.error=data[,'PrecursorError.ppm.'],charge=data[,'Charge'],
		      search.engine="MSGF+",score=data[,'MSGFScore'],stringsAsFactors=FALSE)
  ib.protnpep <-  .convert.msgfp.protein(data[,'Protein'],ib.df[,'peptide'],filter.rev.hits=filter.rev.hits)
  merge(ib.df,ib.protnpep,by='peptide',all=TRUE)
}

.convert.msgfp.pepmodif <- function(peptide,modif.masses=
                                      c(iTRAQ4plex=144.102,Cys_CAM=57.021,Oxidation=15.995,
                                        iTRAQ4_ACET="144.102+42.011",ACET=42.011,
                                        METH=14.016,BIMETH=28.031,TRIMETH=42.047
                                        )) {
  modif.masses.r <- setNames(c(names(modif.masses)),c(paste0("+",modif.masses)))
  modifs <- strsplit(paste0(peptide,"+"),"[A-Z]")
  #sort(table(unlist(modifs)))
  modif.p <- sapply(modifs,function(x) { 
                                         x[length(x)] <- sub(".$","",x[length(x)]);
                                         x[x!=""] <- modif.masses.r[x[x!=""]];
                                         x
                                      })
  if (any(is.na(unlist(modif.p))))
    stop("NA in modifications - did not map")
  #sel <- sapply(modif.p,function(x) any(is.na(x)))
  #modifs[sel]
  
  pep.seq <- gsub("[+0-9\\.]","",peptide)
  if (any(nchar(pep.seq)+1 != sapply(modif.p,length)))
    stop("some modif sequences are not of correct length!")
  #cbind(peptide,pep.seq,nchar(pep.seq),sapply(modif.p,length))
  
  modif.i <- mapply(function(pep,modif) {
    pep <- c("Nterm",pep)
    sel.m <- modif!="" & !grepl("_",modif)
    modif[sel.m] <- paste0(modif[sel.m],"_",pep[sel.m])
    paste(modif,collapse=":")
  } ,strsplit(pep.seq,""),modif.p)

  
  
  return(data.frame(peptide=pep.seq,modif=paste0(modif.i,":"),stringsAsFactors=FALSE))
}


.convert.msgfp.protein <- function(protein,peptide,filter.rev.hits=TRUE) {
  # layout: sp|Q60848-1|HELLS_MOUSE(pre=R,post=K);sp|Q60848-2|HELLS_MOUSE(pre=R,post=K)
  # reverse hits: XXX_sp|...
  # TODO: extract pre and post AAs

  protein.acs <- lapply(strsplit(protein,"[;|]"),
                        function(y) { 
                          acs <- y[seq(from=2,to=length(y),by=3)] 
                          if (filter.rev.hits) {
                            is.rev <- grepl("^XXX_",y[seq(from=1,to=length(y),by=3)])
                            acs[!is.rev]
                          } else
                            acs
                        })
  if (filter.rev.hits) {
    fwd.prots <- sapply(protein.acs,length) > 0
    protein.acs <- protein.acs[fwd.prots]
    peptide <- peptide[fwd.prots]
  }

  as.data.frame(do.call(rbind,lapply(seq_along(peptide),
          function(i) cbind(peptide=peptide[i],accession=protein.acs[[i]]) )),stringsAsFactors=FALSE)
}


