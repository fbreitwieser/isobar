### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor and readIBSpectra.
###

.check.columns <- function(identifications,data.ions=NULL,data.mass=NULL,allow.missing.columns=FALSE) {
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
               paste(missing.cols,collapse="\n\t"),
               "\nData:\n")
      #paste(capture.output(print(head(identifications))))
      if (allow.missing.columns) {
         warning(msg)
         print(head(identifications))
         identifications[,missing.cols] <- NA
      } else {
         message(msg)
         print(head(identifications))
         stop()
      }
  }
  message(" data.frame columns OK")
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

.get.quant.elems <- function(assayDataElements,data.ions,data.mass,spectra.ids,fragment.precision,reporterTagMasses) {
  if (is.null(assayDataElements))
    assayDataElements <- list()

  assayDataElements$ions <- data.ions
  assayDataElements$mass <- data.mass
 
  
  if (!identical(spectra.ids,rownames(data.ions))) {
    id.n.quant <- intersect(spectra.ids,rownames(data.mass))
    id.not.quant <- setdiff(spectra.ids,rownames(data.mass))
    quant.not.id <- setdiff(rownames(data.mass),spectra.ids)

    if (length(id.n.quant)==0) stop("No spectra could be matched between identification and quantifications.",
                                    "If the search was preformed with Mascot, set readIBSpectra argument decode.titles=TRUE,\n",
                                    " or in the properties file for report generation:\n",
                                    "   readIBSpectra.args=list(decode.titles=TRUE)")
    
    if (length(quant.not.id) > 0)
      message(" for ",length(quant.not.id)," spectra ",
              "[",round(length(quant.not.id)*100/nrow(data.mass),2),"%],",
              " quantitative information is available,\n",
              "   but no peptide-spectrum match. Spectrum titles: \n\t",
              paste(quant.not.id[1:min(length(quant.not.id),2)],collapse=",\n\t"),", ...")
    if (length(id.not.quant) > 0)
      message(" for ",length(id.not.quant)," spectra ",
              "[",round(length(id.not.quant)*100/length(spectra.ids),2),"%]",
              " with assigned peptides,\n",
              "   no reporter intensities are available. spectrum titles: \n\t",
              paste(id.not.quant[1:min(length(id.not.quant),2)],collapse=",\n\t"),", ...")

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
    
    min.masses <- reporterTagMasses - fragment.precision/2
    max.masses <- reporterTagMasses + fragment.precision/2
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
  if (all(apply(is.na(assayDataElements$mass),2,all))) stop("Unable to extract any reporter m/z values. Try setting 'fragment.precision' higher (current value: ",fragment.precision,"), or to NULL if no filtering is desired.")
  if (all(apply(is.na(assayDataElements$ions),2,all))) stop("Unable to extract any reporter intensities.")

  return(assayDataElements)
}

.get.varMetaData <- function(identifications) {
  nn <- .SPECTRUM.COLS %in% colnames(identifications)

  label.desc <- c(PEPTIDE='peptide sequence',
      MODIFSTRING='modifications of peptide',
      CHARGE='peptide charge state',
      THEOMASS='theoretical peptide mass',
      EXPMASS='experimental peptide mass',
      EXPMOZ='experimental mass-to-charge ratio',
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
      SCAFFOLD.PEPPROB="Scaffold: Peptide Probability",
      SEQUEST.XCORR="sequest:xcorr",
      SEQUEST.DELTACN="sequest:deltacn",

      IS.DECOY="is decoy peptide-spectrum match?",
      P.VALUE="p value",
      MSGF.RAWSCORE="MSGF raw score",
      MSGF.DENOVOSCORE="MSGF de novo score",
      MSGF.SPECEVALUE="MSGF specturm EValue",
      MSGF.EVALUE="MSGF EValue",
      MSGF.QVALUE="MSGF QValue",
      MSGF.PEPQVALUE="MSGF PepValue",

      PHOSPHO.SITES="phosphorylation sites",
      USEFORQUANT='use spectrum for quantification',
      SEQPOS='PTM seqpos',
      SITEPROBS='PhosphoRS site.probs',
      PEP.SITEPROBS='PhosphoRS site.probs inside peptide',
      FILE='file',SAMPLE='sample',NOTES='notes'
      )
  if (!all(names(.SPECTRUM.COLS) %in% names(label.desc))) {
    sel.bad <- !names(.SPECTRUM.COLS) %in% names(label.desc)
    warning("Not all SPECTRUM COLS have a label description:\n\t",
            paste(names(.SPECTRUM.COLS)[sel.bad],collapse="\n\t"))
    label.desc[names(.SPECTRUM.COLS)[sel.bad]] <- .SPECTRUM.COLS[sel.bad]
  }


  VARMETADATA=data.frame(labelDescription=label.desc[names(.SPECTRUM.COLS)],
                         row.names=.SPECTRUM.COLS)
	    
  return(VARMETADATA[nn,,drop=FALSE])
}

.remove.duplications <- function(identifications) {
  identifications <- unique(identifications)
  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  if (max(table(identifications[,SC['SPECTRUM']])) > 1) {
    t <- table(identifications[,SC['SPECTRUM']])
    warning(sum(t>1)," spectra have diverging identifications after the merging, removing them.")
    print(head(identifications[identifications[,SC['SPECTRUM']] %in% names(t)[t>1],]))

    identifications <- identifications[!identifications[,SC['SPECTRUM']] %in% names(t)[t>1],]
  }

  if ('SEARCHENGINE' %in% names(SC)) {
    message("  Identification details:")
    tt <- sort(table(identifications[,SC['SEARCHENGINE']]))
    stats <- data.frame(perc=sprintf("%.2f %%",tt/sum(tt)*100),n=tt)
    if ('SCORE' %in% names(SC)) {
      scores <- identifications[,SC['SCORE']]
      score.stats <- do.call(rbind,lapply(names(tt),function (se) {
        my.scores <- scores[identifications[,SC['SEARCHENGINE']]==se]
        summary(my.scores,na.rm=TRUE)}))
      stats <- cbind(stats,score.stats)
    } 
    print(stats)
    
  }
  identifications
}

.merge.identifications.full <- function(identifications, ...) {
  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  ## Substitute Isoleucins with Leucins (indistinguishable by Masspec)
  if (!.PEPTIDE.COLS['REALPEPTIDE'] %in% colnames(identifications)) 
     identifications <- .fix.il.peptide(identifications)

  ## Separate protein columns (focus on peptide-spectrum matches)
  PC <- unique(c(.SPECTRUM.COLS['PEPTIDE'],.PROTEIN.COLS,.PEPTIDE.COLS))
  protein.colnames <- colnames(identifications)[colnames(identifications) %in% c(PC)]
  pept.n.prot <- unique(identifications[,protein.colnames])
  identifications <- unique(identifications[,-which(colnames(identifications) %in% setdiff(PC,.SPECTRUM.COLS['PEPTIDE']))])

  ## Merge identifications
  if (max(table(identifications[,SC['SPECTRUM']])) > 1) {
    if (SC['DISSOCMETHOD'] %in% colnames(identifications) && 
        length(unique(identifications[,SC['DISSOCMETHOD'] ]))) {
      identifications <- dlply(identifications,'dissoc.method',.merge.identifications,...)
      # rbind separately to assure equal column names
      identifications <- do.call(rbind,identifications)
      identifications <- .merge.quant.identifications(identifications)
    } else {
      identifications <- .merge.identifications(identifications, ...)
    }
  }
  identifications <- .remove.duplications(identifications)
  
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
  if (!.PEPTIDE.COLS['REALPEPTIDE'] %in% colnames(identifications)) 
     identifications <- .fix.il.peptide(identifications)

  ## Separate protein columns (focus on peptide-spectrum matches)
  PC <- unique(c(.SPECTRUM.COLS['PEPTIDE'],.PROTEIN.COLS,.PEPTIDE.COLS))
  protein.colnames <- colnames(identifications)[colnames(identifications) %in% c(PC)]
  pept.n.prot <- unique(identifications[,protein.colnames])
  identifications <- unique(identifications[,-which(colnames(identifications) %in% setdiff(PC,.SPECTRUM.COLS['PEPTIDE']))])
  ## Merge identifications
  if (max(table(identifications[,SC['SPECTRUM']])) > 1) {
    message("merging identifications")
    if (SC['DISSOCMETHOD'] %in% colnames(identifications) && 
        length(unique(identifications[,SC['DISSOCMETHOD'] ]))) {
      identifications <- dlply(identifications,'dissoc.method',.merge.identifications,...)
      # rbind separately to assure equal column names
      identifications <- do.call(rbind,identifications)
      identifications <- .merge.quant.identifications(identifications)
    } else {
      identifications <- .merge.identifications(identifications, ...)
    }
  }
  identifications <- .remove.duplications(identifications)
 
  # Create ProteinGroup
  proteinGroup <- ProteinGroup(merge(pept.n.prot,identifications,by="peptide"),template=proteinGroupTemplate)
  
  ## Get intensities and masses in assayDataElements
  if (is.null(data.ions))
    data.ions <- .get.quant(identifications,.PEAKS.COLS['IONSFIELD'],reporterTagNames)
  if (is.null(data.mass))
    data.mass <- .get.quant(identifications,.PEAKS.COLS['MASSFIELD'],reporterTagNames)

  if (!all(rownames(data.ions)==rownames(data.mass),na.rm=TRUE))
    stop(sum(rownames(data.ions)==rownames(data.mass))," rownames are not equal between ions and mass")
  if (any(is.na(rownames(data.ions))))
    stop(sum(is.na(rownames(data.ions)))," rownames in data.ions are NA")
  if (any(is.na(rownames(data.mass))))
    stop(sum(is.na(rownames(data.ions)))," rownames in data.ions are NA")

  assayDataElements <- .get.quant.elems(assayDataElements,data.ions,data.mass,
                                        identifications[,SC['SPECTRUM']],fragment.precision,
                                        reporterTagMasses(.Object))

  if (!.SPECTRUM.COLS['USEFORQUANT'] %in% colnames(identifications)) {
    # not perfect - better: set spectra of peptides shared between groups to FALSE
    #                            and spectra with no values
    identifications[,.SPECTRUM.COLS['USEFORQUANT']] <- TRUE
  }  

  fdata <- identifications[,colnames(identifications) %in% .SPECTRUM.COLS]

  rownames(fdata) <- fdata[,'spectrum']
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
    function(type,id.file,identifications.format=NULL,
             sep="\t",decode.titles=FALSE,trim.titles=FALSE,...) {
      new(type,
          identifications=.read.idfile(id.file,sep=sep,
                                       identifications.format=identifications.format,
                                       decode.titles=decode.titles,trim.titles=trim.titles),...)
    }
)
setMethod("readIBSpectra",
          signature(type="character",id.file="data.frame",peaklist.file="missing"),
    function(type,id.file,...) new(type,identifications=id.file,...)
)

setMethod("readIBSpectra",
          signature(type="character",id.file="character",peaklist.file="character"),
    function(type,id.file,peaklist.file,sep="\t",
             mapping.file=NULL,mapping=c(quantification.spectrum = "hcd",identification.spectrum = "cid"),
             id.file.domap=NULL,identifications.format=NULL,decode.titles=FALSE,...) {
      
      id.data <- .read.identifications(id.file,sep=sep,
                                       mapping=mapping.file,mapping.names=mapping,
                                       identifications.quant=id.file.domap,
                                       identifications.format=identifications.format,decode.titles=decode.titles)

      readIBSpectra(type,id.data,peaklist.file,...)
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
##' @param identifications.format 
##' @return IBSpectra object of type 
##' @author Florian P Breitwieser
setMethod("readIBSpectra",
          signature(type="character",id.file="data.frame",peaklist.file="character"),
    function(type,id.file,peaklist.file,
             annotate.spectra.f=NULL,
             peaklist.format=NULL,scan.lines=0,
             fragment.precision=NULL,fragment.outlier.prob=NULL,...) {
      
      if (is.function(annotate.spectra.f)) {
        id.file <- annotate.spectra.f(id.file,peaklist.file)
      }

      .Object <- new(type)
      quant <- .read.peaklist(peaklist.file,peaklist.format,
                              .Object@reporterTagMasses,.Object@reporterTagNames,
                              scan.lines,fragment.precision,fragment.outlier.prob,
                              id.data=id.file) 
     
      new(type,identifications=id.file,quant[[2]],data.ions=quant[[1]],...)
    }
)

# returns a list with intensities and mass matrices
.read.peaklist <- function(peaklist.file,peaklist.format,
                           reporterMasses,reporterTagNames,
                           scan.lines,fragment.precision,fragment.outlier.prob,
                           id.data=NULL) {
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
      message("  reading peaklist ",peaklist.f," ...",appendLF=FALSE)
      res <- read.delim(peaklist.f,stringsAsFactors=FALSE,quote="")
      message(" done")
      if (!"spectrum" %in% colnames(res)) stop("'spectrum' column is missing from peaklist CSV")
      if (sum(.grep_columns(res,"ions$")) < length(reporterTagNames)) 
        stop("Not all neccessary intensity columns [",paste0(reporterTagNames,"_ions"),"] in peaklist CSV")
      if (sum(.grep_columns(res,"mass$")) < length(reporterTagNames)) 
        stop("Not all neccessary reporter mass columns [",paste0(reporterTagNames,"_mass"),"] in peaklist CSV")
      
      data.titles <- c(data.titles,res[,"spectrum"])
      data.ions <- rbind(data.ions,as.matrix(res[,.grep_columns(res,"ions$")]))
      data.mass <- rbind(data.mass,as.matrix(res[,.grep_columns(res,"mass$")]))
    }
  }

  if (!is.null(id.data) && .SPECTRUM.COLS['SPECTRUM.QUANT'] %in% colnames(id.data)) {
    data.titles.orig <- data.titles
    data.titles <- .do.map(data.titles,unique(id.data[,.SPECTRUM.COLS[c('SPECTRUM','SPECTRUM.QUANT')]]))
    sel.na <- is.na(data.titles)
    if (any(is.na(data.titles))) {
      message(" for ",sum(sel.na)," of ",length(data.titles)," spectra,",
              " quantitative information is available,\n",
              "   but no peptide-spectrum match. Spectrum titles: \n\t",
              paste(data.titles.orig[sel.na][1:2],collapse=",\n\t"),", ...")
      data.ions <- data.ions[!sel.na,]
      data.mass <- data.mass[!sel.na,]
      data.titles <- data.titles[!sel.na]
    }
  }

  rownames(data.ions)  <- data.titles
  rownames(data.mass)  <- data.titles
  ## TODO: check that all identified spectra are present in intensities

  colnames(data.ions) <- reporterTagNames
  colnames(data.mass) <- reporterTagNames

  if (any(is.na(rownames(data.ions))))
    stop(sum(is.na(rownames(data.ions)))," rownames in data.ions are NA")
  if (any(is.na(rownames(data.mass))))
    stop(sum(is.na(rownames(data.ions)))," rownames in data.ions are NA")
 
  return(list(data.ions,data.mass))
}


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
read.mzid <- function(filename) {
  library(XML)
  doc <- xmlInternalTreeParse(filename)
  ns <- c(x=xmlNamespace(xmlRoot(doc))[[1]])

  root <- ifelse (isTRUE(ns == "http://psidev.info/psi/pi/mzIdentML/1.0"), "/x:mzIdentML", "/x:MzIdentML")
  peptidesequence.name <- ifelse (isTRUE(ns == "http://psidev.info/psi/pi/mzIdentML/1.0"), "peptideSequence", "PeptideSequence")
  dbsequenceref.attrname <- ifelse (isTRUE(ns == "http://psidev.info/psi/pi/mzIdentML/1.0"), "DBSequence_ref", "dBSequence_ref")
  peptideevref.attrname <- ifelse (isTRUE(ns == "http://psidev.info/psi/pi/mzIdentML/1.0"), "PeptideEvidence_ref", "peptideEvidence_ref")
  peptideref.attrname <- ifelse (isTRUE(ns == "http://psidev.info/psi/pi/mzIdentML/1.0"), "Peptide_ref", "peptide_ref")

  searchdatabase.mapping <- data.frame(
    ref=xpathSApply(doc,paste0(root,"/x:DataCollection/x:Inputs/x:SearchDatabase"),
      xmlGetAttr,name="id",namespaces=ns),
    name=xpathSApply(doc,paste0(root,"/x:DataCollection/x:Inputs/x:SearchDatabase"),
      xmlGetAttr,name="name",namespaces=ns),stringsAsFactors=FALSE
  )
   
  modification.mapping.df <- t(do.call(cbind,
                                  xpathApply(doc,paste0(root,"/x:AnalysisProtocolCollection/x:SpectrumIdentificationProtocol/x:ModificationParams/x:SearchModification"),
                                             function(modif) {

    res <- xmlAttrs(modif)
    modparams <- getNodeSet(modif,"x:ModParam",namespaces=ns)
    return (switch(as.character(length(modparams)),
            "0" = c(res,xmlAttrs(modif[['cvParam']])),
            "1" = c(res,xmlAttrs(modparams[[1]]),xmlAttrs(modparams[[1]][['cvParam']])),
            stop("Expecting zero or one ModParam Node in ",
                 "/mzIdentML/AnalysisProtocolCollection/SpectrumIdentificationProtocol/ModificationParams/SearchModification")))
  },namespaces=ns)))

  unknown.modif <- modification.mapping.df[,'name']=='unknown modification' | 
                     modification.mapping.df[,'accession']=='MS:1001460'
  if (any(unknown.modif)) {
    if (! 'value' %in% colnames(modification.mapping.df)) {
      stop('unknown modifications [as defined by nameor accession], and no column value to override.')
    }
    modification.mapping.df[unknown.modif,'name'] <- modification.mapping.df[unknown.modif,'value']
  }

  modif.map <- .as.vect(unique(modification.mapping.df[,c('massDelta','name')]))

  message("peptide and protein mapping")
  peptide.mapping <- do.call(rbind,xpathApply(doc,paste0(root,"/x:SequenceCollection/x:Peptide"),namespaces=ns,
         function(pep) {
           peptide.ref <- xmlGetAttr(pep,name='id')
           peptide.seq <- xmlValue(pep[[peptidesequence.name]])

           loc.n.delta <- xpathSApply(pep,"x:Modification",function(m) c(location=xmlGetAttr(m,name="location"),
                                                                         massDelta=xmlGetAttr(m,name="monoisotopicMassDelta")),
                                      namespaces=ns)

          
           modifstring <- character(nchar(peptide.seq)+2)

           if (length(loc.n.delta) > 0) {
             modif.names <- modif.map[loc.n.delta[2,]]
             if (any(is.na(modif.names)))
               modif.names <- xpathSApply(pep,"x:Modification/x:cvParam[@cvRef='UNIMOD']",xmlGetAttr,name='name',namespaces=ns)

             if (length(modif.names) < ncol(loc.n.delta)) 
               stop("modif names and delta have different size")
 
             modifstring[as.numeric(loc.n.delta[1,])+1] <- modif.names
           }
           
           c(peptide.ref=peptide.ref,peptide=peptide.seq,modif=paste(modifstring,collapse=":"))
         }))

  protein.mapping <- do.call(rbind,xpathApply(doc,paste0(root,"/x:SequenceCollection/x:DBSequence"),function(dbs) {
    c(dbseq.ref = xmlGetAttr(dbs,name='id'),
      accession = xmlGetAttr(dbs,name='accession'),
      length = xmlGetAttr(dbs,name='length'),
      sdb.ref = xmlGetAttr(dbs,name='SearchDatabase_ref'))
      #sequence = xmlValue(dbs[['seq']]))
  },namespaces=ns))
  records.proteinDetections <- do.call(rbind,xpathApply(doc,
    paste0(root,"/x:DataCollection/x:AnalysisData/x:ProteinDetectionList/x:ProteinAmbiguityGroup"),
    namespaces=ns,
    fun=function(pag) {
    ## apply on ProteinAmbiguityGroup
    do.call(rbind,xmlApply(pag,function(pdh) {
      ## apply on ProteinDetectionHypothesis
      cbind(dbseq.ref=xmlGetAttr(pdh,dbsequenceref.attrname),
            peptide.ev.ref=xpathSApply(pdh,"x:PeptideHypothesis",xmlGetAttr,name=peptideevref.attrname,namespaces=ns)
    )}))
  }))

  # map PeptideEvidence attribute names to isobar columns
  pe.attr.names <- c(id="peptide.ev.ref",start="start.pos",end="end.pos",pre="aa.before",post="aa.after",
                     missedCleavages="nmc",isDecoy="is.decoy",DBSequence_Ref="dbseq.ref",dBSequence_Ref="dbseq.ref")
  
  message("spectrum mapping")
  spectrumIdentifications <-
    xpathApply(doc,paste0(root,"/x:DataCollection/x:AnalysisData/x:SpectrumIdentificationList/x:SpectrumIdentificationResult"),
              namespaces=ns,
              fun=function(sir) {
    ## _SpectrumIdentificationResult_ #
    ## All identification from searching one spectrum
    spectrum.id <- xmlGetAttr(sir,"spectrumID")
    spectrum.title.nodes <- getNodeSet(sir,"x:cvParam[@name='spectrum title']",namespaces=ns)
    if (length(spectrum.title.nodes) == 1)
      spectrum.id <- xmlGetAttr(spectrum.title.nodes[[1]],name='value')

    ## get spectrum identifications which pass threshold
      ## _SpectrumIdentificationItem_ #
      ## An identification of a single peptide of a specturm.
      ## Only take the one which passes the threshold.
    sii <- getNodeSet(sir,"x:SpectrumIdentificationItem[@passThreshold='true' and @rank='1']",namespaces=ns)
    if (length(sii) == 0 || length(sii) > 1) return (NULL)
    #"more than one match for [ x:SpectrumIdentificationItem[@passThreshold='true' and @rank='1'] ]")
    # TODO: There will be always two times the same score with peptides with just an I/L difference
    
    sii <- sii[[1]]
 
     scores <- unlist(xpathApply(sii,"x:cvParam",function(cvp) {
                           switch(xmlGetAttr(cvp,name='name'),
                                  "Scaffold: Peptide Probability"=c(scaffold.pepprob=xmlGetAttr(cvp,name="value")),
                                  "mascot:score"=c(score.mascot=xmlGetAttr(cvp,name="value")),
                                  "mascot:expectation value"=c(mascot.evalue=xmlGetAttr(cvp,name="value")),
                                  "sequest:xcorr"=c(sequest.xcorr=xmlGetAttr(cvp,name="value")),
                                  "sequest:deltacn"=c(sequest.deltacn=xmlGetAttr(cvp,name="value")),
                                  NULL
                            )
                 },namespaces=ns))

      if ("PeptideEvidence" %in% names(sii)) {
        pe.attr <- xmlAttrs(sii[["PeptideEvidence"]])
        pe.res <- pe.attr.names[names(pe.attr)[names(pe.attr) %in% names(pe.attr.names)]]
      } else {
        pe.res <- c(peptide.ev.ref=xmlGetAttr(sii[["PeptideEvidenceRef"]],name="peptideEvidence_ref"))
      }

      c(spectrum=spectrum.id,
        peptide.ref   = xmlGetAttr(sii,name=peptideref.attrname),
        theo.mass     = xmlGetAttr(sii,name='calculatedMassToCharge'),
        exp.mass      = xmlGetAttr(sii,name='experimentalMassToCharge'),
        scores,pe.res)
  })

  si.names <- unique(unlist(sapply(spectrumIdentifications,names)))
  records.spectrumIdentifications <- do.call(rbind,lapply(spectrumIdentifications,function(x) x[si.names]))
  records.spectrumIdentifications <- as.data.frame(records.spectrumIdentifications,stringsAsFactors=FALSE)

  peptide.ev.mapping <- xpathApply(doc,paste0(root,"/x:SequenceCollection/x:PeptideEvidence"),namespaces=ns,xmlAttrs)
  if (length(peptide.ev.mapping) > 0) {
    records.spectrumIdentifications <- merge(records.spectrumIdentifications,
                                             as.data.frame(do.call(rbind,peptide.ev.mapping),stringsAsFactors=FALSE))
  }

  spectra.n.peptides <- merge(as.data.frame(records.spectrumIdentifications,stringsAsFactors=FALSE),
                              as.data.frame(peptide.mapping,stringsAsFactors=FALSE),
                              by='peptide.ref')

  free(doc)

  message("merging results")
  protein.mapping1 <- merge(as.data.frame(records.proteinDetections,stringsAsFactors=FALSE),
                            as.data.frame(protein.mapping,stringsAsFactors=FALSE),by="dbseq.ref")
  merge(spectra.n.peptides,
        protein.mapping1,by="peptide.ev.ref")
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
  message("  ",nSpectra," spectra in MGF file.")

  ## create list with all spectra (header+mass list) as entries
  all.spectra <- apply(bnd,1,function(x) input[x[1]:x[2]])
  rm(input)

  ## if all spectra are of equal length, apply returns a matrix
  ##  convert to list
  if (is.matrix(all.spectra)) 
    all.spectra <- split(all.spectra,col(all.spectra))

  ## extract information from each spectrum
  result <- llply(all.spectra,function(x) {
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
        if (length(pos) > 1){
          # (Jacques) added to address cases where 2 peaks (within the required precision) are at the same distance of the reporter theoretical mass
          # Pick most intense
          max.intens <- max(mzi.ions[pos])
          pos <- pos[which(mzi.ions[pos]==max.intens)]
          pos <- pos[1] # in case same intensity as well!
        }
        if (length(pos) > 0 & m[pos] < fragment.precision/2)
          return(c(mzi.mass[pos],mzi.ions[pos]))
        else
          return(c(NA,NA))
      }))) 
    } else {
      rr <- c(rr,rep(NA,nReporter*2))
    }
    return(rr)
  },.parallel=isTRUE(getOption('isobar.parallel')))

  result <- do.call(rbind,result)

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
    message("\tmass boundaries:\n\t",
            paste(colnames(mass),sprintf("%.5f : %.5f",bnd[1,],bnd[2,]),sep="\t",collapse="\n\t"))
  }

  ions <- ions[sel,,drop=FALSE]
  mass <- mass[sel,,drop=FALSE]
 
  spectrumtitles <- .trim(result[sel,1])

  if (ncol(ions) != length(reporterNames) || ncol(mass) != length(reporterNames)) {
    stop("ions or mass matrix have wrong dimension.")
  }

  dimnames(ions) <- list(spectrumtitles,reporterNames)
  dimnames(mass) <- list(spectrumtitles,reporterNames)
  rm(result)
  
  return(list(ions=ions, mass=mass,
              spectrumtitles=spectrumtitles))
}

.read.idfile.df <- function(filename,identifications.format,sep,...) {
    if (is.data.frame(filename))
      return(filename)
    if (!is.character(filename))
      stop("filename should be a data.frame or character")
    if (!file.exists(filename))
      stop("idfile ",filename," defined, but does not exist")
    
    is.format <- function(y,ext)
      identical(identifications.format,y) || 
        any(sapply(ext,function(iext) grepl(paste0(iext,"$"),filename)))
    
    ext.def <- 
    list(mzid      = list(ext=c("mzid"),f=read.mzid),
         rockerbox = list(ext=c("peptides.[ct]sv"),f=.read.rockerbox),
         msgf      = list(ext=c("msgfp.csv","tsv"),f=.read.msgfp.tsv),
         ibspectra = list(ext=c("csv"),
                          f=function(x) read.table(x,header=TRUE,sep=sep,
                                                   stringsAsFactors=FALSE,quote="",...)))

    for (etype in names(ext.def)) {
      if (is.format(etype,ext.def[[etype]][["ext"]])) {
        message("  reading id file ",filename," [type: ",etype,"] ...",appendLF=FALSE)
        res <- ext.def[[etype]][["f"]](filename)
        message(" done")
        return(res)
      }
    }

    stop(paste("cannot parse file ",filename," - cannot deduce format based on extension ",
               "(it is not ibspectra.csv, id.csv, peptides.txt or mzid). ",
               "Please provide id.format to readIBSpectra",sep=""))
}


## TODO: log is not returned
.read.idfile <- function(id.file,identifications.format=NULL,sep="\t",
                         decode.titles=FALSE,trim.titles=FALSE,log=NULL,all.cols=FALSE,...) {
  if (!is.data.frame(id.file)) {
    if (!(is.list(id.file) || is.character(id.file))) 
      stop("id.file argument of .read.idfile should be a data frame, character or list")

    if (all(sapply(id.file,is.character)))
      id.data <- lapply(id.file,.read.idfile.df,sep=sep,
                        identifications.format=identifications.format,...)
    else if (all(sapply(id.file,is.data.frame)))
      id.data <- id.file
    else
      stop()

    rm(id.file)
  
    id.colnames <- lapply(seq_along(id.data),function(s.i) colnames(id.data[[s.i]]))
    colnames.equal <- all(sapply(id.colnames,
                                 function(cn) identical(cn,id.colnames[[1]])))
    if (!colnames.equal) {
      message(" id file colnames are not equal:")
      message(paste(sapply(seq_along(id.data),
                           function(s.i) paste0("    ",s.i,": [",
                                                paste(id.colnames[[s.i]],collapse=","),
                                                "]")),collapse="\n"))
      all.id.colnames <- unique(unlist(id.colnames))
      if (isTRUE(all.cols)) {
        id.data <- do.call(rbind,lapply(id.data,function(i.d) {
          id.d[,all.id.colnames[!all.id.colnames %in% colnames(id.d)]] <- NA
          id.d[,all.id.colnames]
        }))
      } else {
        intersect.colnames <- id.colnames[[1]]
        for (s.i in seq(from=2,to=length(id.colnames))) {
          intersect.colnames <- intersect(intersect.colnames,id.colnames[[s.i]])
        }
        message(" taking intersection: ",paste(intersect.colnames,collapse="; "))
        id.data <- do.call(rbind,lapply(id.data,function(i.d) i.d[,intersect.colnames]))
      }
    } else {
        id.data <- do.call(rbind,id.data)
    }
  }

  #.check.columns(id.data)

  if (!is.character(id.data[,.SPECTRUM.COLS['SPECTRUM']]))
    id.data[,.SPECTRUM.COLS['SPECTRUM']] <- as.character(id.data[,.SPECTRUM.COLS['SPECTRUM']])

  if ('accession' %in% colnames(id.data)) {
    if (all(grepl(sprintf("%s\\|%s",.uniprot.pattern.ac,.uniprot.pattern.id),id.data[,.PROTEIN.COLS['PROTEINAC']]))) {
      split.ac <- strsplit(id.data[,.PROTEIN.COLS['PROTEINAC']],"|",fixed=TRUE)
      id.data[,.PROTEIN.COLS['PROTEINAC']] <- sapply(split.ac,function(x) x[1]) 
      id.data[,.PROTEIN.COLS['NAME']] <- sapply(split.ac,function(x) ifelse(length(x) == 2, x[2], NA))
    }
  }

  if (decode.titles)
    id.data[,.SPECTRUM.COLS['SPECTRUM']] <- sapply(id.data[,.SPECTRUM.COLS['SPECTRUM']],URLdecode)

  if (trim.titles)
    id.data[,.SPECTRUM.COLS['SPECTRUM']] <- .trim(id.data[,.SPECTRUM.COLS['SPECTRUM']])
  return(id.data)
}

.read.rockerbox <- function(filename) {
  data.r <- read.table(filename,header=T,stringsAsFactors=F,sep="\t")
  if (!"scan.title" %in% colnames(data.r))
    stop("no scan.title column in ",filename,"; please use a Rockerbox version >= 2.0.6")

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

  sel <- names(.ROCKERBOX.MAPPING.COLS) %in% names(c(.SPECTRUM.COLS,.PEPTIDE.COLS))
  data.r <- data.r[,.ROCKERBOX.MAPPING.COLS[sel]]
  colnames(data.r) <- c(.SPECTRUM.COLS,.PEPTIDE.COLS)[names(.ROCKERBOX.MAPPING.COLS)[sel]]
  return(data.r)
}

###############################################################################
## MERGE IDENTIFCATIONS


.dissect.search.engines <- function(identifications) {

  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]

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
  score.columns <- paste0('score.',tolower(unique(unlist(engines))))
  for (engine in unique(unlist(engines))) {
    name <- paste0('score.',tolower(engine))
    e.scores <- mapply(function(e,s) if(any(e==engine)) s[e==engine] else NA,engines,scores)
    identifications[,name] <- as.numeric(e.scores)
  }

  identifications[,SC['SEARCHENGINE']] <- sapply(engines,paste,collapse="&")
  tt <- table(identifications[,SC['SPECTRUM']])
  if (max(tt) > 1) {
    message("  resolving duplicated ids ...")
    id.good <- identifications[identifications$spectrum %in% names(tt)[tt==1],]
    id.bad <- identifications[identifications$spectrum %in% names(tt)[tt>1],]
    id.bad$n.se <- rowSums(!is.na(id.bad[,score.columns]))
    id.bad <- ddply(id.bad,'spectrum',function(x) {
      x[which.max(x$n.se),]
    })
    id.bad$n.se <- NULL
    identifications <- rbind(id.good,id.bad)
  }

  if ('SCORE' %in% names(SC))
    identifications[,SC['SCORE']] <- NULL

  return(identifications)
}

.merge.search.engine.identifications <- function(identifications,...) {

  ## Columns to consolidate after merging
  ##   functions are used to resolve or remove the psm state
  COLS.TO.CONSOLIDATE <- list('peptide'=.consolidate.peptide.ids,
			      'modif'=.consolidate.modification.pos,
			      'charge'=.consolidate.charge,  ## sometimes, charge state 0 is reported when too high
			      'theo.mass'=.mean.na.rm,        ## can differ esp. between Mascot and Phenyx
			      'retention.time'=.mean.na.rm,
			      'parent.intens'=.mean.na.rm)

  COLS.TO.CONSOLIDATE <- COLS.TO.CONSOLIDATE[names(COLS.TO.CONSOLIDATE) %in% colnames(identifications)]

  ## Columns used for merging. Typically 'spectrum' and columns which are the same regardless of search engine
  COLS.FOR.MERGING <- setdiff(.SPECTRUM.COLS,c(.ID.COLS,names(COLS.TO.CONSOLIDATE)))
  COLS.FOR.MERGING <- c(COLS.FOR.MERGING,'is.decoy')
  COLS.FOR.MERGING <- intersect(colnames(identifications),COLS.FOR.MERGING)

  #'delta.score','delta.score.pep','delta.score.notpep','n.pep','n.loc'

  ## Split identifications by search engine
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
    cn[!cn %in% COLS.FOR.MERGING] <- paste(cn[!cn %in% COLS.FOR.MERGING],clean.names[ii],sep=".")
    colnames(ids.split[[ii]]) <- cn
  }

  ## Merge identification data frames
  ids.merged <- merge(ids.split[[1]],ids.split[[2]],by=COLS.FOR.MERGING,all=TRUE)
  if (length(ids.split) > 2) {
    for (ii in seq(from=3,to=length(ids.split))) {
      ids.merged <- merge(ids.merged,ids.split[[ii]],by=COLS.FOR.MERGING,all=TRUE)
    }
  }

  if (!.SPECTRUM.COLS['NOTES'] %in% colnames(ids.merged))
    ids.merged[,.SPECTRUM.COLS['NOTES']] <- ''

  ## Consolidate spectra with different identifications with different search engines
  for (colname in names(COLS.TO.CONSOLIDATE)) {
    message("Consolidating ",colname, " [",date(),"]")
    ids.merged <- .resolve.conflicts(ids=ids.merged,resolve.f=COLS.TO.CONSOLIDATE[[colname]],
  	         		     colname=colname,resolve.colnames=paste0(colname,".",clean.names))
  }

  ids.merged <- unique(ids.merged)
  tt <- table(ids.merged[,.SPECTRUM.COLS['SPECTRUM']])
  #spectra.ok <- ids.merged[,'spectrum'] %in% names(tt)[tt==1]
  #resolved.ids <- .resolve.differing.identifications(ids.merged[!spectra.ok,],score.cols)

  #identifications <- rbind(ids.merged[spectra.ok,],resolved.ids)
  if (any(tt>1)) {
    warning("Merging not successful, ",sum(tt>1)," duplicated psms, removing them!")
    ids.merged <- ids.merged[!ids.merged[,.SPECTRUM.COLS['SPECTRUM']] %in% names(tt)[tt>1],]
  }

  ids.merged[,.SPECTRUM.COLS['SEARCHENGINE']] <- 
	  apply(ids.merged[,score.cols],1,function(x) { paste(names(ids.split)[!is.na(x)],collapse="&") })
  return(ids.merged)
}

.na.rm <- function(x) x[!is.na(x)]
.mean.na.rm <- function(x) mean(x,na.rm=TRUE)

# Remove differing peptide ids
.consolidate.peptide.ids <- function(x) NA

# Take a non-zero charge
.consolidate.charge <- function(x) {
	x <- !is.na(x)
	ifelse(any(x > 0),x[x>0][1],0)
}

# Take 'first' modification
.consolidate.modification.pos <- function(x) .na.rm(x)[1]

.searchengine.cols <- function(ids,colname) {

}

.resolve.conflicts <- function(ids, resolve.f, colname, resolve.colnames, keep.cols = FALSE) {

  if (is.character(resolve.colnames))
    resolve.colnames <- which(colnames(ids) %in% resolve.colnames)

  # set result column
  ids[,colname] <- .return.equal.or.na(ids[,resolve.colnames])

  # return resolved conflicting and non-conflicting data
  good.spectra <- !is.na(ids[,colname])
  if (!any(good.spectra))
    stop("No good spectra which don't have to be resolved - something went wrong")
  if (!all(good.spectra)) {
    print(isobar:::.sum.bool(good.spectra))
    ids <- rbind(ids[good.spectra,],
		 .do.resolve.conflicts(ids[!good.spectra,],colname,resolve.colnames, resolve.f))
  }

  if (keep.cols)
    ids
  else
    ids[-resolve.colnames,]
}

.do.resolve.conflicts <- function(ids, colname, resolve.colnames, resolve.f) {
  if (length(ids) == 0)
    return(NULL) 
  if (length(ids) == 1)
    return(ids)

  ids[,colname] <- apply(ids[,resolve.colnames],1,resolve.f)
  if (any(is.na(ids[,colname])))
    message("removing ",sum(is.na(ids[,colname]))) 
  ids <- ids[!is.na(ids[,colname]),]
  if (nrow(ids) == 0)
	  return(NULL)

  ids[,.SPECTRUM.COLS['NOTES']] <- ifelse(ids[,.SPECTRUM.COLS['NOTES']] == '','',
					  paste0(ids[,.SPECTRUM.COLS['NOTES']],"\n"))

  ids[,.SPECTRUM.COLS['NOTES']] <- paste0(ids[,.SPECTRUM.COLS['NOTES']],
					  apply(ids[,resolve.colnames],1,function(x) {
					    x <- x[!is.na(x)]
					    paste(colname,"differing:\n\t",
                                                   paste(names(x),x,sep=": ",collapse="\n\t"))}))
  return(ids)
}

.resolve.modifications <- function(df,colname, modif.cols, standard.modif) {
  message("Resolving modifications")
  if (is.character(modif.cols))
    modif.cols <- which(colnames(df) %in% modif.cols)

  adply(df, 1, function(x) {

    x <- as.data.frame(x)
    modifs <- unlist(x[,modif.cols])
    modifs <- modifs[!is.na(modifs)]
    if (length(modifs) == 0) {
	    print(x)
	    stop ("all modifs are NA!")
    }

    ## coherent modifications
    x$note <- ""
    if (length(modifs) == 1 || all(modifs == modifs[1])) {
      stop("Should be handled before")
      x[,colname] <- modifs[1]
      return(x[,-modif.cols])
    }

    ## resolve incoherent modifs
    note <- paste0("differing modification positions: \n\t",
		   paste(names(modifs),modifs,sep=": ",collapse="\n\t"))

    ## any modification position is more often present?
    tt <- table(modifs)
    if (.max.uniq(tt)) {
      note <- paste0(note,". Taking ", .take.max(tt) )
      stop("TODO: Implement")
    }

    ## take 'first' modification position
    x[,colname] <- modifs[1]
    x[,'note'] <- note
    return(x)
  })  
}

.max.uniq <- function(tt) sum(tt==max(tt)) == 1
.take.max <- function(tt) names(tt)[which.max(tt)]

.resolve.differing.identifications <- function(identifications,score.cols) {
  allequal <- function(x) all(x == x[1])
  by.y <- function(x,ind,fun=sum) sapply(by(x,ind,fun),function(x) return(x) )
  skipna <- function(x) unlist(x[!is.na(x)])
  n.skipped <- 0
  n.max.ids <- 0
  n.modif.pos.dif <- 0

  resolved.identifications <- ddply(identifications,'spectrum', function(x) {
	n.ids <- rowSums(!is.na(x[,score.cols]),na.rm=TRUE)    ## number of identifications for each psm
	if (.max.uniq(n.ids)) {
	  n.max.ids <<- n.max.ids + 1
	  return(x[which.max(n.ids),,drop=FALSE])
	}
        
	## resolve peptide differnces
        if (!allequal(x[,'peptide'])) {                ## different peptides identified
	  pep.ids <- by.y(n.ids,x[,'peptide'],sum)
	  if (.max.uniq(pep.ids)) {                     ##   any peptide has been seen more often than others?
	    x <- x[x[,'peptide'] == .take.max(pep.ids),,drop=FALSE]
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

	  modifs <- gsub("Oxidation_M","")

	  if (!.max.uniq(modif.ids))
	    n.modif.pos.dif <<- n.modif.pos.dif + 1

	  #print(x)
	  modif <- paste(x[,.SPECTRUM.COLS['SEARCHENGINE']],x[,'modif'],sep=": ",collapse=' & ')
	  x <- x[x[,'modif'] == .take.max(modif.ids)[1],,drop=FALSE]
	  x[,'modif'] <- modif
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

  tt <- table(identifications[,SC['SPECTRUM.QUANT']])
  spectra.ok <- identifications[,SC['SPECTRUM.QUANT']] %in% names(tt)[tt==1]

  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  score.colname <- .SPECTRUM.COLS[c('SCORE.MASCOT','SCORE.PHENYX','SCORE.MSGF','PROB.PHOSPHORS')]
  score.colname <- score.colname[score.colname %in% colnames(identifications)]

  message("Merging ",sum(!spectra.ok)," identifications from different dissociation methods.")
  ids.quant.merged <- ddply(identifications[!spectra.ok,],.SPECTRUM.COLS['SPECTRUM.QUANT'],function(x) {
      if (nrow(x) == 1) return(x)

      my.args <- as.list(x[,score.colname,drop=FALSE])
      my.args <- lapply(my.args,round,digits=2) ## take two significant digits before taking next score into account as top hitter
      my.args$decreasing=TRUE

      max.hit <- do.call(order,my.args)[1]

      # dismiss peptide-spectrum matches which are different between search engines
      if (!all(x[,SC['PEPTIDE']] == x[1,SC['PEPTIDE']]))
        return(NULL)

      x[max.hit,SC['DISSOCMETHOD']] = paste0("[",x[max.hit,SC['DISSOCMETHOD']],"]")
      if (all(x[,SC['MODIFSTRING']] == x[max.hit,SC['MODIFSTRING']]) && 'DISSOCMETHOD' %in% names(SC)) {
      if ("cid" %in% x[,SC['DISSOCMETHOD']])
      x[max.hit,SC['SPECTRUM']] <- x[x[,SC['DISSOCMETHOD']]=="cid",SC['SPECTRUM']][1]

      x[max.hit,SC['DISSOCMETHOD']] <- paste(x[,SC['DISSOCMETHOD']],collapse="&")
      for (sc in score.colname[score.colname != 'pepprob'])
        x[max.hit,sc] <- gsub("NA","",paste(x[,sc],collapse="&"))
      }
      return(x[max.hit,,drop=FALSE])
  })

  colnames(ids.quant.merged)[colnames(ids.quant.merged)=='SPECTRUM.QUANT'] <- 'spectrum.quant'
  ids.quant.merged <- ids.quant.merged[,colnames(identifications)]
  rbind(identifications[spectra.ok,],ids.quant.merged)
}

.merge.identifications <- function(identifications,...) {
  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]

  ## Merge results of different search engines / on different spectra
  if ('SEARCHENGINE' %in% names(SC) &&                               ## search engine column is present
      length(unique(identifications[,SC['SEARCHENGINE']])) > 1 &&    ## there are multiple search engines defines
      !any(grepl('^score\\.',colnames(identifications))))              ## the individual score.ENGINENAME are not present
  {                     
    if (any(grepl("|",identifications[,SC['SEARCHENGINE']],fixed=TRUE))) {
      identifications <- .dissect.search.engines(identifications)               ## dissect merged id columns
    } else {
      identifications <- .merge.search.engine.identifications(identifications,...)  ## merge identifications
    }
  }

  return(identifications)
}

## end MERGE IDENTIFICATIONS
###############################################################################


.read.identifications <- function(identifications,...,
                                  mapping=NULL,mapping.names=c(quantification.spectrum="hcd",identification.spectrum="cid"),
                                  identifications.quant=NULL) {

  ## load identifications (is either character or data.frame
  if (is.character(identifications) && all(sapply(identifications,file.exists)))
    identifications <- .read.idfile(identifications,...)

  ## load mapping (either character or data.frame)
  if (!is.null(mapping)) {
    if (is.character(mapping) && file.exists(mapping))
      mapping <- do.call(rbind,lapply(mapping,read.table,sep=",",header=TRUE,stringsAsFactors=FALSE,quote=""))
    if (!is.data.frame(mapping)) stop("mapping must be a data.frame or valid file name")
    if (ncol(mapping) > 2) stop("only one column in mapping")

    if (is.null(mapping.names)) { 
      if (ncol(mapping) > 2) stop("more than one column in mapping data - mapping.names are therefore required")
      message("trying to assess right mapping")
      cn <- apply(mapping,2,function(x) sum(x %in% identifications[,.SPECTRUM.COLS['SPECTRUM']]))
      mapping.names <- c(quantification.spectrum=names(cn)[which.min(cn)],
                         identification.spectrum=names(cn)[which.max(cn)])
    }
    if (!all(c('identification.spectrum','quantification.spectrum') %in% names(mapping.names)))
      stop("mapping.names must be given alongside with mapping. names(mapping.names) must be 'identification.spectrum' and 'quantification.spectrum'")
    if (!all(mapping.names %in% colnames(mapping))) stop("mapping.names must be column names of mapping")

    quant2id <- .as.vect(mapping,mapping.names['identification.spectrum'],mapping.names['quantification.spectrum'])
    id2quant <- .as.vect(mapping,mapping.names['quantification.spectrum'],mapping.names['identification.spectrum'])

    identifications[,.SPECTRUM.COLS['SPECTRUM.QUANT']] <- id2quant[identifications[,.SPECTRUM.COLS['SPECTRUM']]]
    identifications[,.SPECTRUM.COLS['DISSOCMETHOD']] <- ifelse(is.null(identifications.quant),"hcd","cid") # temp fix
    #identifications[,.SPECTRUM.COLS['DISSOCMETHOD']] <- names(mapping.names)[mapping.names=='identification.spectrum']

    if (!is.null(identifications.quant)) {
      ## load identifications.quant (either character or data.frame)
      if (is.character(identifications.quant) && all(sapply(identifications.quant,file.exists)))
        identifications.quant <- .read.idfile(identifications.quant,...)
      if (!is.data.frame(identifications.quant)) stop("identifications.quant must be a data.frame or valid file name")
      identifications.quant[,.SPECTRUM.COLS['SPECTRUM.QUANT']] <- identifications.quant[,.SPECTRUM.COLS['SPECTRUM']]
      identifications.quant[,.SPECTRUM.COLS['SPECTRUM']] <- quant2id[ identifications.quant[,.SPECTRUM.COLS['SPECTRUM.QUANT']] ]
      identifications.quant[,.SPECTRUM.COLS['DISSOCMETHOD']] <- "hcd"
      #identifications.quant[,.SPECTRUM.COLS['DISSOCMETHOD']] <- names(mapping.names)[mapping.names=='quantification.spectrum']
      #if (.SPECTRUM.COLS['DATABASE'] %in% colnames(identifications) &&
      #    !.SPECTRUM.COLS['DATABASE'] %in% colnames(identifications.quant)) 
      #  identifications[,.SPECTRUM.COLS['DATABASE']] <- NULL

      intersect.ids <- intersect(colnames(identifications),colnames(identifications.quant))
      cn.i.intersect <- colnames(identifications)[!colnames(identifications) %in% intersect.ids]
      cn.iq.intersect <- colnames(identifications.quant)[!colnames(identifications.quant) %in% intersect.ids]
      if (length(cn.i.intersect) > 0 || length(cn.iq.intersect)>0)
        stop("identifications and identifications.quant dataframes are different:",
             ifelse(length(cn.i.intersect),paste0("\n\t[",paste0(cn.i.intersect,collapse=";"),
                                                  "] only appear in identifications"),""),
             ifelse(length(cn.iq.intersect),paste0("\n\t[",paste0(cn.iq.intersect,collapse=";"),
                                                   "] only appear in identifications.quant"),""))

      identifications <- rbind(identifications[,intersect.ids],identifications.quant[,intersect.ids])
    }
  }
  return(identifications)
}


## read MSGF+ tab-separated identification files
.read.msgfp.tsv <- function(filename,filter.rev.hits=FALSE) {
  if (is.data.frame(filename)) 
    ib.data <- filename
  else
    id.data <- read.delim(filename, sep="\t", stringsAsFactors=FALSE,quote="")

  if (! "PepQValue" %in% colnames(id.data)) {
    stop("No q-values found in MSGF+ result files.\nPlease repeat the MSGF+ search using a\ntarget-decoy search (option -tda 1)")
  }

  ib.df <- data.frame(Protein=id.data[,'Protein'],
                      spectrum=id.data[,'Title'],spectrum.quant=id.data[, 'Title'],
                      .convert.msgfp.pepmodif(id.data[,'Peptide']),
                      scan.from=id.data[,'ScanNum'],dissoc.method=tolower(id.data[,'FragMethod']),
                      precursor.error=id.data[,'PrecursorError.ppm.'],charge=id.data[,'Charge'],
                      search.engine="MSGF+",score=-log10(id.data[,'SpecEValue']),
                      stringsAsFactors=FALSE)

  sel <- names(.MSGF.MAPPINGCOLS) %in% names(c(.SPECTRUM.COLS,.PEPTIDE.COLS))
  data.r <- id.data[,.MSGF.MAPPINGCOLS[sel]]
  colnames(data.r) <- c(.SPECTRUM.COLS,.PEPTIDE.COLS)[names(.MSGF.MAPPINGCOLS)[sel]]
  
  ib.df <- cbind(ib.df,data.r)

  ib.protnpep <-  .convert.msgfp.protein(unique(ib.df[,c('Protein','peptide')]),filter.rev.hits=filter.rev.hits)
  id.data$Protein <- NULL
  merge(unique(ib.protnpep),ib.df,by='peptide',all=TRUE)
}

.convert.msgfp.pepmodif <- function(peptide,modif.masses=
                                      c(iTRAQ4plex=144.102,Cys_CAM=57.021,Oxidation=15.995,
                                        iTRAQ4_ACET="144.102+42.011",ACET=42.011,
                                        TMT6plex=229.163,
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


.convert.msgfp.protein <- function(protein.n.peptide,filter.rev.hits=FALSE) {
  # layout: sp|Q60848-1|HELLS_MOUSE(pre=R,post=K);sp|Q60848-2|HELLS_MOUSE(pre=R,post=K)
  # reverse hits: XXX_sp|...
  # TODO: extract pre and post AAs

  split.protein <- strsplit(protein.n.peptide[,1],"[;|]")
  res <- lapply(seq_along(split.protein), function(i) { 
                y <- split.protein[[i]]
                cbind(peptide=protein.n.peptide[i,2],
                      database=y[seq(from=1,to=length(y),by=3)],
                      accession=y[seq(from=2,to=length(y),by=3)])
                        })
  res <- as.data.frame(do.call(rbind,res),stringsAsFactors=FALSE)
  res$is.decoy <- grepl("^XXX_",res[,'database'])
  print(c(is.decoy=.sum.bool(res$is.decoy)))
  
  if (isTRUE(filter.rev.hits))
    res <- res[!res$is.decoy,]

  return(res)
}


