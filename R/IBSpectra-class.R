### =========================================================================
### IBSpectra objects
### -------------------------------------------------------------------------
###
### Class definition.

setClass("IBSpectra",
         contains = "eSet",
         representation(
             proteinGroup = "ProteinGroup",
             reporterTagNames = "character",
             reporterTagMasses = "numeric",
             isotopeImpurities ="matrix",
             log = "matrix",
             "VIRTUAL"),
         prototype = 
            prototype(log=matrix(c(format(Sys.time(), "%F %T %Z"),"IBSpectra Log"),
                                 ncol=2,
                                 dimnames=list("init",c("Timestamp","Message"))))
)

setClass("iTRAQSpectra",
    contains = "IBSpectra",
    representation("VIRTUAL")
)

setClass("iTRAQ4plexSpectra",
    contains = "iTRAQSpectra",
    prototype = prototype(
        reporterTagNames = as.character(114:117),
        reporterTagMasses = c(114.1112,115.1083,116.1116,117.1150),
        isotopeImpurities = matrix(1/100*c(
                c( 92.90,  5.90,  0.20,  0.00 ),
                c(  2.00, 92.30,  5.60,  0.10 ),
                c(  0.00,  3.00, 92.40,  4.50 ),
                c(  0.00,  0.10,  4.00, 92.30 )),
            nrow=4,byrow=TRUE)
    )
)

setClass("iTRAQ8plexSpectra",
    contains = "iTRAQSpectra",
    prototype = prototype(
      reporterTagNames = as.character(c(113:119,121)),
      reporterTagMasses = c(113.1078,114.1112,115.1082,116.1116,
                         117.1149,118.1120,119.1153,121.1220),
      isotopeImpurities = diag(nrow=8)
    )
)

setClass("TMTSpectra",
    contains = "IBSpectra",
    representation("VIRTUAL")
)

setClass("TMT2plexSpectra",
    contains = "TMTSpectra",
    prototype = prototype(
      reporterTagNames = c("126","127"),
      reporterTagMasses = c(126.127725,127.131079),
      isotopeImpurities = diag(nrow=2)
    )
)

setClass("TMT6plexSpectra",
    contains = "TMTSpectra",
    prototype = prototype(
      reporterTagNames = as.character(126:131),
      reporterTagMasses = c(126.127725,127.131079,128.134433,
                         129.137787,130.141141,131.138176),
       isotopeImpurities = diag(nrow=6)
       )
    )


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.IBSpectra.slots <- function(object) {
    n <- length(object@reporterTagNames)
    if (length(object@reporterTagMasses) != n) 
        stop("[IBSpectra: validation] Slot reportMasses [",length(object@reporterTagMasses),"]",
             " has not the same length as reporterTagNames [",n,"]!")

    if (!all(dim(object@isotopeImpurities) == n))
        stop("[IBSpectra: validation] isotopeImpurities has dimensions ",
             paste(dim(object@isotopeImpurities),collapse="x"),"!",
             " Expected is ",n,"x",n,". Set it to diag(",n,") when no known.")

       if (length(classLabels(object)) > 0 && length(classLabels(object)) != n) 
        stop("[IBSpectra: validation] Slot classLabels [",length(classLabels(object)),"]",
             " has not the same length as reporterTagNames [",n,"]!")

    return(TRUE)
}

.valid.IBSpectra <- function(object) {
	.valid.IBSpectra.slots(object)
}

setValidity("IBSpectra",.valid.IBSpectra)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Definitions.
###

# Definition of column names corresponding to spectrum,
# peptide, and protein level data.

.PEAKS.COLS <- c(MASSFIELD="X%s_mass",IONSFIELD="X%s_ions")

.ID.COLS <- c(SEARCHENGINE="search.engine",SCORE="score",SCORE.MASCOT="score.mascot",SCORE.PHENYX="score.phenyx",SCORE.MSGF="score.msgf")

.PTM.COLS <- c(SCORE.PHOSPHORS='pepscore',PROB.PHOSPHORS='pepprob',SEQPOS='seqpos',SITEPROBS='site.probs',PHOSPHO.SITES='phospho.sites')

.SPECTRUM.COLS <- c(PEPTIDE="peptide",MODIFSTRING="modif",CHARGE="charge",
                   THEOMASS="theo.mass",EXPMASS="exp.mass",PRECURSOR.ERROR="precursor.error",
                   PARENTINTENS="parent.intens",RT="retention.time",
                   SPECTRUM="spectrum",SPECTRUM.QUANT="spectrum.quant",
                   .ID.COLS,USEFORQUANT="use.for.quant",
                   .PTM.COLS,
                   DISSOCMETHOD="dissoc.method",
                   PRECURSOR.PURITY="precursor.purity",
                   SCANS="scans",SCANS.FROM="scan.from",SCANS.TO="scan.to",
                   RAWFILE="raw.file",NMC="nmc",DELTASCORE="delta.score",DELTASCORE.PEP="delta.score.pep",
                   MASSDELTA.ABS="massdelta.abs",MASSDELTA.PPM="massdelta.ppm",
                   SAMPLE="sample",FILE="file",NOTES="notes")

.PEPTIDE.COLS <- c(PROTEINAC="accession",STARTPOS="start.pos",
                  REALPEPTIDE="real.peptide",AA.BEFORE="aa.before",AA.AFTER="aa.after")

.PROTEIN.COLS <- c(PROTEINAC="accession",PROTEINAC_CONCISE="accessions",
                   NAME="name",PROTEIN_NAME="protein_name",
                   GENE_NAME="gene_name",ORGANISM="organism")

.ROCKERBOX.COLS <- c(PROTEINAC="accession",STARTPOS="start.pos",MODIFSTRING="modif",
                     QN="mascot.query.number",RANK="rank",SCANS="scan.number.s.",RT="retention.time",
                     RAWFILE="raw.file",PEPTIDE="sequence", "phosphosequence",NMC="miscleavages",
                     SCORE="score",DELTASCORE="deltascore",PHOSPHO.DELTASCORE.PEP="phospho_deltascore",
                     EXPMASS="peptide.mr",MASSDELTA.ABS="mass.delta..abs.",MASSDELTA.PPM="mass.delta..ppm.",
                     ROCKERBOX_MOD="modifications",PROTEINACS="all.protein.matches",SPECTRUM="scan.title")


writeIBSpectra <- function(ibspectra,file,sep="\t",row.names=FALSE,...) {
  write.table(as.data.frame(ibspectra),file=file,sep=sep,row.names=row.names,...)
  message("finished writing ibspectra to ",file)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("IBSpectra","data.frame",
    function(from) {

      # prepare ProteinGroup data.frame
      pg.df <- as(proteinGroup(from),"data.frame")
      pg.df <- pg.df[,c("protein","peptide","spectrum","start.pos")]
      colnames(pg.df)[colnames(pg.df) == "protein"] <- .PROTEIN.COLS['PROTEINAC']
      colnames(pg.df)[colnames(pg.df) == "peptide"] <- .SPECTRUM.COLS['PEPTIDE']
      colnames(pg.df)[colnames(pg.df) == "spectrum"] <- .SPECTRUM.COLS['SPECTRUM']

      # prepare fData data.frame
      ri <-reporterIntensities(from)
      colnames(ri) <- paste0("X",reporterTagNames(from),"_ions")
      rm <- reporterMasses(from)
      colnames(rm) <- paste0("X",reporterTagNames(from),"_mass")

      fdata.df <- cbind(fData(from),rm,ri)
      fdata.df[,.SPECTRUM.COLS['PEPTIDE']] <- NULL

      # merge data.frames
      res <- unique(merge(pg.df,fdata.df,by=.SPECTRUM.COLS['SPECTRUM'],all.y=TRUE,sort=FALSE))

      # return in order
      res[order(res[,.PROTEIN.COLS['PROTEINAC']],res[,.PEPTIDE.COLS['STARTPOS']],res[,.SPECTRUM.COLS['PEPTIDE']]),
          c(intersect(c(.PROTEIN.COLS,.PEPTIDE.COLS,.SPECTRUM.COLS),colnames(res)),
            colnames(rm),colnames(ri))]
    }
)

ibSpectra.as.concise.data.frame  <- function(from) {

      # prepare ProteinGroup data.frame
      pg.df <- proteinGroup.as.concise.data.frame(proteinGroup(from))
      colnames(pg.df)[colnames(pg.df) == "proteins"] <- .PROTEIN.COLS['PROTEINAC_CONCISE']
      colnames(pg.df)[colnames(pg.df) == "peptide"] <- .SPECTRUM.COLS['PEPTIDE']

      # prepare fData data.frame
      ri <-reporterIntensities(from)
      colnames(ri) <- paste0("X",reporterTagNames(from),"_ions")
      rm <- reporterMasses(from)
      colnames(rm) <- paste0("X",reporterTagNames(from),"_mass")
      fdata.df <- cbind(fData(from),rm,ri)

      # merge data.frames
      res <- unique(merge(pg.df,fdata.df,by=.SPECTRUM.COLS['PEPTIDE'],all=TRUE,sort=FALSE))

      # return in order
      res[order(res[,.PROTEIN.COLS['PROTEINAC_CONCISE']],res[,.SPECTRUM.COLS['PEPTIDE']]),
          c(intersect(c(.PROTEIN.COLS,.PEPTIDE.COLS,"n.groups","n.acs","n.variants",.SPECTRUM.COLS),colnames(res)),
            colnames(rm),colnames(ri))]
    }

.IBSpectraAsConciseDataFrameNew <- function(from,show.phospho.position=FALSE) {
  # prepare ProteinGroup data.frame
  indist.proteins <- indistinguishableProteins(proteinGroup(from))
      pg.df <- proteinGroup.as.concise.data.frame(from)
  colnames(pg.df)[colnames(pg.df) == "proteins"] <- .PROTEIN.COLS['PROTEINAC_CONCISE']
  colnames(pg.df)[colnames(pg.df) == "peptide"] <- .SPECTRUM.COLS['PEPTIDE']

  # prepare fData data.frame
  ri <-reporterIntensities(from)
  colnames(ri) <- paste0("X",reporterTagNames(from),"_ions")
  rm <- reporterMasses(from)
  colnames(rm) <- paste0("X",reporterTagNames(from),"_mass")
  fdata.df <- cbind(fData(from),rm,ri)

  # merge data.frames
  res <- unique(merge(pg.df,fdata.df,by=.SPECTRUM.COLS['PEPTIDE'],all=TRUE,sort=FALSE))

  if (show.phospho.position) {
    res.nice <- data.frame(Sequence=.convertPeptideModif(res$peptide,res$modif),
                           Phospho.Position=.convertModifToPos(res$modif,"PHOS"))
  } else {
    res.nice <- data.frame(Sequence=res$peptide)
  }
  res.nice <- cbind(res.nice,
                    Modification=res$modif,
                    AC=.protein.acc(res[,"accession"],ip=indist.proteins),
                    ID=proteinInfo(protein.group,res[,"accession"],do.warn=FALSE),
                    n=sapply(res[,"accession"],function(p) {length(names(indist.proteins)[indist.proteins == p])}))
  res <- res[,!.grep_columns(res,c("peptide","modif","accession"))]
  # return in order
  res.nice <- cbind(res.nice,res[,c(intersect(.SPECTRUM.COLS,colnames(res)),
                                    colnames(rm),colnames(ri))])
  res.nice[order(res$AC,res$Sequence,res$modif),]
}

setMethod("as.data.frame",signature(x="IBSpectra"), 
		function(x, row.names=NULL, optional=FALSE, ...) as(x,"data.frame"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

IBSpectraTypes <- function() { 
  unlist(sapply(
          names(getClass("IBSpectra")@subclasses),
          function(c) { if (!isVirtualClass(c)) c; } )
  )
}

setGeneric("proteinGroup", function(x) standardGeneric("proteinGroup"))
setGeneric("isotopeImpurities", function(x) standardGeneric("isotopeImpurities"))
setGeneric("proteinGroup<-", function(x,value) standardGeneric("proteinGroup<-"))
setGeneric("isotopeImpurities<-", function(x,value) standardGeneric("isotopeImpurities<-"))


# reporterTagNames in package Biobase are defunct now
setGeneric("reporterTagNames",function(object) standardGeneric("reporterTagNames"))
setMethod("reporterTagNames","IBSpectra",function(object) object@reporterTagNames)
setGeneric("reporterTagMasses",function(object) standardGeneric("reporterTagMasses"))
setMethod("reporterTagMasses","IBSpectra",function(object) object@reporterTagMasses)
setMethod("proteinGroup","IBSpectra", function(x) x@proteinGroup)
setMethod("isotopeImpurities","IBSpectra", function(x) x@isotopeImpurities)

setReplaceMethod("proteinGroup","IBSpectra",function(x,value) {
      x@proteinGroup <- value
      x
    }
)
setReplaceMethod("isotopeImpurities",signature("IBSpectra"),
    function(x,value) {
       x@isotopeImpurities <- value
       validObject(x)
       x
    }
)


# TODO: add replaceMethods for reporterMasses and reporterIntensities
#       with a selection of peptide or protein

setGeneric("reporterData", function(x,...) standardGeneric("reporterData"))
setGeneric("reporterData<-", function(x,...,value) standardGeneric("reporterData<-"))
setGeneric("reporterIntensities", function(x,...)
           standardGeneric("reporterIntensities"))
setGeneric("reporterIntensities<-", function(x,...,value)
           standardGeneric("reporterIntensities<-"))
setGeneric("reporterMasses", function(x,...) standardGeneric("reporterMasses"))
setGeneric("reporterMasses<-", function(x,...,value)
           standardGeneric("reporterMasses<-"))

setMethod("reporterData","IBSpectra",
    function(x,element="ions",na.rm=FALSE,na.rm.f='any',...) {
      sel <- spectrumSel(x,...)
      data <- assayDataElement(x,element)[sel,,drop=FALSE]

      if (na.rm & length(data) > 0)
        return(data[!apply(is.na(data),1,na.rm.f),,drop=FALSE])
      else
        return(data)
    }
    )

setReplaceMethod("reporterData","IBSpectra",
    function(x,element="ions",...,value) {
      sel <- spectrumSel(x,...)
      assayDataElement(x,element)[sel,] <- value
      x
    }
)

setMethod("reporterIntensities","IBSpectra",
		function(x,...) reporterData(x,...,element="ions"))

setReplaceMethod("reporterIntensities",signature("IBSpectra"),
    function(x,...,value) {
	 	reporterData(x,...,element="ions") <- value
      x
    }
)

setMethod("reporterMasses","IBSpectra",
		function(x,...) reporterData(x,...,element="mass"))

setReplaceMethod("reporterMasses",signature("IBSpectra"),
    function(x,...,value) {
	 	reporterData(x,...,element="mass") <- value
      x
    }
)


writeData <- function(x, element="ions", file=paste(element,"csv",sep="."), quote=FALSE,
        sep="\t", ...){
      col.names <- sampleNames(x)
      row.names <- spectrumTitles(x)
      write.table(reporterData(x,element="ions"), file=file, quote=quote, sep=sep,
          col.names=col.names, row.names=row.names, ...)
    }


# TODO: use protein.g instead of protein (when protein group identifier is meant)
setGeneric("spectrumSel", function(x,peptide,protein,...) standardGeneric("spectrumSel"))
setMethod("spectrumSel",signature(x="IBSpectra",peptide="missing",protein="missing"),
    function(x) rep(TRUE,nrow(fData(x))))

setMethod("spectrumSel",signature(x="IBSpectra",peptide="data.frame",protein="missing"),
    function(x,peptide,...) spectrumSel(x,as.matrix(peptide),...) )

setMethod("spectrumSel",signature(x="IBSpectra",peptide="matrix",protein="missing"),
    function(x,peptide,modif=NULL,spectrum.titles=FALSE,use.for.quant.only=TRUE,do.warn=TRUE) {
        if (length(peptide) == 0) {
          warning("0L peptide provided")
          return(FALSE)
        }
        if (ncol(peptide) != 2 && do.warn)
          warning("don't know how to handle matrix with ",ncol(peptide)," columns! expecting 2.")
        
        sel <- fData(x)[,.SPECTRUM.COLS['PEPTIDE']]  %in% peptide[,1] & 
               fData(x)[,.SPECTRUM.COLS['MODIFSTRING']]  %in% peptide[,2]

        if (use.for.quant.only && .SPECTRUM.COLS['USEFORQUANT'] %in% colnames(fData(x)))
          sel <- sel & fData(x)[,.SPECTRUM.COLS['USEFORQUANT']] 
        
        for (m in modif)
          sel <- sel & grepl(m,fData(x)[,.SPECTRUM.COLS['MODIFSTRING']])

        if (!any(sel) && do.warn) warning("No spectra for peptide ",paste(peptide,collapse="; modif = "))
        if (spectrum.titles)
          return(rownames(fData(x))[sel])
        else
  	      return(sel)
    }
)

setMethod("spectrumSel",signature(x="IBSpectra",peptide="character",protein="missing"),
    function(x,peptide,modif=NULL,spectrum.titles=FALSE,use.for.quant.only=TRUE,do.warn=TRUE) {
        if (length(peptide) == 0) {
          warning("0L peptide provided")
          return(FALSE)
        }
        sel <- fData(x)[,.SPECTRUM.COLS['PEPTIDE']]  %in% peptide
        for (m in modif)
          sel <- sel & grepl(m,fData(x)[,.SPECTRUM.COLS['MODIFSTRING']])

        if (use.for.quant.only && .SPECTRUM.COLS['USEFORQUANT'] %in% colnames(fData(x)))
          sel <- sel & fData(x)[,.SPECTRUM.COLS['USEFORQUANT']] 
        if (!any(sel) && do.warn) warning("No spectra for peptide ",peptide)
         if (spectrum.titles)
          return(rownames(fData(x))[sel])
        else
  	      return(sel)
    }
)

setMethod("spectrumSel",signature(x="IBSpectra",peptide="missing",protein="character"),
    function(x,protein,specificity=REPORTERSPECIFIC,modif=NULL,spectrum.titles=FALSE,use.for.quant.only=TRUE,
             do.warn=TRUE,...) {
      
      peptides <- peptides(x=proteinGroup(x),protein=protein,specificity=specificity,do.warn=do.warn,...)
      if (length(peptides) == 0)
        return(FALSE)
      sel <- spectrumSel(x,peptide=peptides,spectrum.titles=spectrum.titles,
                         modif=modif,use.for.quant.only=use.for.quant.only,do.warn=do.warn)
      if (do.warn && (isTRUE(spectrum.titles) && length(sel) == 0 || !isTRUE(spectrum.titles) && !any(sel)))
        warning("No spectra for protein ",protein,
                " with specificity ",paste(specificity,collapse=","))
      return(sel)
    }
)

setGeneric("spectrumTitles", function(x,...) standardGeneric("spectrumTitles"))
setMethod("spectrumTitles","IBSpectra",function(x,...) 
    function(x) fData(x)[spectrumSel(x,...),.SPECTRUM.COLS['SPECTRUM']])

#setMethod("reporterTagNames","IBSpectra",function(x) x@reporterTagNames)

setGeneric("classLabels", function(x) standardGeneric("classLabels"))
setGeneric("classLabels<-", function(x,value) standardGeneric("classLabels<-"))

setMethod("classLabels",signature(x="IBSpectra"), 
    function(x) {
      class.labels <- phenoData(x)[["class.labels"]]
      names(class.labels) <- phenoData(x)[["class.label.description"]]
      return(class.labels)
    }
)
setReplaceMethod("classLabels","IBSpectra",
    function(x,value) {
      
      if (length(value) != length(reporterTagNames(x)))
        stop("Class labels [",length(value),"] need to habe the same length as reporterTagNames [",length(reporterTagNames(x)),"]")

      if (is.null(names(value))) {
        phenoData(x)[["class.labels",labelDescription="class labels"]] <- value
      } else {
        phenoData(x)[["class.labels",labelDescription="class labels"]] <- as.character(value)
        phenoData(x)[["class.label.description",labelDescription="class label descriptions"]] <- names(value)
      }

      validObject(x)

      x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### correctIsotopeImpurities, normalize, 
### substractAdditiveNoise and exclude methods.
###

setGeneric("correctIsotopeImpurities",
           function(x) standardGeneric("correctIsotopeImpurities") )
setGeneric("subtractAdditiveNoise",
           function(x,...) standardGeneric("subtractAdditiveNoise") )
setGeneric("exclude",
           function(x,protein,peptide=NULL,...) standardGeneric("exclude"))

setMethod("correctIsotopeImpurities",signature(x="IBSpectra"),
    function(x) {
      
      if (is.logged(x,"isotopeImpurities.corrected")) {
        warning("isotope impurities have been corrected already.");
      }
      
      ri <- reporterIntensities(x)
      AA <- isotopeImpurities(x)
      ri.corrected <- t(apply(ri,1,function(b) {
            ok <- !is.na(b)
            if (sum(ok) > 1){
              A <- AA[ok,ok]
              x <- base::solve(A,b[ok])
              b[ok] <-x
            }
            return(b)
          }
      ))
      colnames(ri.corrected) <- reporterTagNames(x)
      ri.corrected[ri.corrected<1] <- NA
      reporterIntensities(x) <- ri.corrected

      x <- do.log(x,"isotopeImpurities.corrected",TRUE)

      return(x)
    }
)

normalize <- function(x,f=median,target="intensity",
                      exclude.protein=NULL, use.protein=NULL, peptide.specificity=NULL,
                      f.doapply=TRUE,log=TRUE,
                      channels=NULL,na.rm=FALSE,per.file=TRUE,
                      normalize.factors=NULL,...){
  
  ## NOTE: median normalizes might normalize too much when a lot of NA
  ##         values are present - mean or colSums better?
  
  if (!is.logged(x,"isotopeImpurities.corrected")) 
    warning("Isotope impurity correction has not been logged",
            " - data might be uncorrected. See ?correctIsotopeImpurities.")

  x <- do.log(x,"is.normalized",TRUE)
  ## save original reporter intensities for noise estimation
  if (is.null(assayDataElement(x,"ions_not_normalized")))
    assayDataElement(x,"ions_not_normalized") <- reporterIntensities(x)


  if (!is.null(normalize.factors)) {
    reporterIntensities(x) <- reporterIntensities(x)/
      rep(normalize.factors,each=nrow(reporterIntensities(x)))
    for (i in seq_along(normalize.factors)) {
      x <- do.log(x,paste("normalization.multiplicative.factor channel",
                          colnames(reporterIntensities(x))[i]),
                  round(normalize.factors[i],4))
    }
    return(x)
  }

  if (per.file) {
    if (!.SPECTRUM.COLS['FILE'] %in% colnames(fData(x))) {
      warning("No 'file' column is specified in IBSpectra, setting normalize per.file to FALSE")
      per.file <- FALSE
    } else if (length(unique(fData(x))[,.SPECTRUM.COLS['FILE']])==1) {
      warning("Only one file - setting normalize per.file to FALSE.")
      per.file <- FALSE
    }
  }

  ##if (is.logged(x,"is.normalized"))
  ##  warning("Normalization is logged already.")
  
  if (!is.null(exclude.protein) & !is.null(use.protein))
    stop("Provide either exclude.protein or use.protein, not both.")

  if (is.null(channels) && length(classLabels(x)) > 0)
    channels <- reporterTagNames(x)[!is.na(classLabels(x))]

  if (!is.null(channels) ) {
    if (is.list(channels)) {
      for (channels.set in channels)
        x <- normalize(x,f=f,target=target,exclude.protein=exclude.protein,
                       use.protein=use.protein,peptide.specificity=peptide.specificity,
                       f.doapply=f.doapply,
                       log=log,channels=channels.set,na.rm=na.rm,...)
      return(x)
    } else {
      if (!all(channels %in% reporterTagNames(x)))
        stop("channels must be reporterTagNames.")
                                                         
      message("normalizing channels ",paste(channels,collapse=", "))
      ri <- reporterIntensities(x)[,channels,drop=FALSE]
    }
  } else {
    ri <- reporterIntensities(x)
  } ## else is.null(channels)
  
  if (!is.null(exclude.protein))
    ri <- ri[!spectrumSel(x,protein=exclude.protein,
                          specificity=if(is.null(peptide.specificity)) c(REPORTERSPECIFIC,GROUPSPECIFIC,UNSPECIFIC) else peptide.specificity),]
  
  if (!is.null(use.protein)) 
    ri <- ri[spectrumSel(x,protein=use.protein,
                         specificity=if(is.null(peptide.specificity)) REPORTERSPECIFIC else peptide.specificity),]

  if (na.rm) 
    sel.na <- apply(!is.na(ri),1,all)
  else
    sel.na <- rep(TRUE,nrow(ri))
  
  ## TODO: warning when ri is empty

  if (per.file && .SPECTRUM.COLS['FILE'] %in% colnames(fData(x))) {
    fd <- fData(x)
    for (n.file in sort(unique(fd[,.SPECTRUM.COLS['FILE']]))) {
      sel <- sel.na & fd[,.SPECTRUM.COLS['FILE']] == n.file
      message("\tnormalizing ",n.file," [",sum(sel)," spectra]")
      ri.sel <- ri[sel,,drop=FALSE]
      normalize.factors <- .get.normalization.factors(ri.sel,f,target,f.doapply,...)
      reporterIntensities(x)[sel,colnames(ri)] <- 
        reporterIntensities(x)[sel,colnames(ri)]/rep(normalize.factors,each=sum(sel))
      for (i in seq_along(normalize.factors)) {
        x <- do.log(x,paste("normalization.multiplicative.factor file",n.file,"channel",colnames(ri)[i]),
                    round(normalize.factors[i],4))
      } 
    }
  } else {
    normalize.factors <- .get.normalization.factors(ri[sel.na,,drop=FALSE],f,target,f.doapply,...)
    reporterIntensities(x)[,colnames(ri)] <- 
      reporterIntensities(x)[,colnames(ri)]/
      rep(normalize.factors,each=nrow(reporterIntensities(x)))
    for (i in seq_along(normalize.factors)) {
      x <- do.log(x,paste("normalization.multiplicative.factor channel",
                          colnames(ri)[i]),
                  round(normalize.factors[i],4))
    }
  }
  ## FIXME: logging a function for f does not work
  ##       x <- do.log(x,"normalization.method",
  ##              sprintf("%s of %s",as.character(substitute(f)),target))
 
  x
}

.get.normalization.factors <- function(ri,f,target,f.doapply,...) {
  
  if (target=="ratio") {
    rs <- rowSums(ri,na.rm=TRUE)
    if (log)
      ri <- apply(ri,2,function(x) x/rs)
    else
      ri <- apply(ri,2,function(x) x-rs)
  } else if (target=="intensity") {
    ## take ri as specified
  } else {
    stop(paste("target",target,"not known"))
  }
  
  if (!f.doapply)
    res <- f(ri,na.rm=TRUE,...)
  else
    res <- apply(ri,2,f,na.rm=TRUE,...)

  return(res/max(res,na.rm=T))
 
}

setGeneric("normalize")


## Estimates the noise level througth the 1st percentile of the >0
## signals (in the linear space) Method quantile: It take's the prob
## (0.01) quantile to estimate the noise level.  This value is
## subtracted from all intensities, and all remaining intensities have
## to be at least that value.

## If channels are assumed similar in intensity and hence a shared
## noise level is reasonable If not, then one level per channel is
## necessary
setMethod("subtractAdditiveNoise",signature(x="IBSpectra"),
    function(x,method="quantile",shared=TRUE,prob=0.01) {
      if (method=="quantile") {
        if (shared) {
          q <- quantile(reporterIntensities(x),probs=prob,na.rm=T)
          ri.sub <- apply(reporterIntensities(x),2,function(x) {
                res <- x - q
                res[res<=q] <- NA
                res
              }
          )
          
        } else {
          ri.sub <- apply(reporterIntensities(x),2,function(x) {
                q <- quantile(x,probs=prob,na.rm=T)
                res <- x - q
                res[res<=q] <- NA
                res
              }
          )
        }
        reporterIntensities(x) <- ri.sub
        x
      } else {
        stop(paste("Method",method,"not implemented"))
      }
    }
)

setMethod("exclude",
    signature(x="IBSpectra",protein="character",peptide="ANY"),
    function(x,protein,peptide=NULL,specificity=c(UNSPECIFIC,GROUPSPECIFIC,REPORTERSPECIFIC)) {
      
      # select peptides for removal
      sel.peptides <- c(peptides(proteinGroup(x),protein=protein,
                                 specificity=specificity),peptide)
      
      # select spectra to keep
      sel.spectra <- !spectrumSel(x,peptide=sel.peptides)
      #sel.spectra <- !spectrumSel(x,protein=proteins.to.exclude,
      #                            specificity=c(REPORTERSPECIFIC,GROUPSPECIFIC,UNSPECIFIC))
      
      # remove spectra from assayData
      for (aden in assayDataElementNames(x)) {
        assayDataElement(x,aden) <- assayDataElement(x,aden)[sel.spectra,]
      }
      
      # remove from featureData
      featureData(x) <- as(fData(x)[sel.spectra,],"AnnotatedDataFrame")
      
      pg.df <- as.data.frame(proteinGroup(x))
      # remove peptides and proteins from proteinGroup
      proteinGroup(x) <- ProteinGroup(pg.df[!pg.df$peptide %in% sel.peptides,] )

      x
    }
)



subsetIBSpectra <- 
  function(x, protein=NULL, peptide=NULL, direction="exclude",
           specificity=c(REPORTERSPECIFIC,GROUPSPECIFIC,UNSPECIFIC),...) {

  if ((is.null(protein) && is.null(peptide)) || (!is.null(protein) && !is.null(peptide)))
    stop("define either protein or peptide to include or exclude")

  if (!is.null(peptide))
    sel.spectra <- spectrumSel(x,peptide=peptide,...)
  else
    sel.spectra <- spectrumSel(x,protein=protein,specificity=specificity,...)

  sel.spectra <- switch (direction,
                         exclude = !sel.spectra,
                         include = sel.spectra,
                         stop("direction must be either 'exclude' or 'include'."))

  for (aden in assayDataElementNames(x)) {
    assayDataElement(x,aden) <- assayDataElement(x,aden)[sel.spectra,]
  }

  pg.df <- as.data.frame(proteinGroup(x))
  # remove peptides and proteins from proteinGroup
  proteinGroup(x) <- ProteinGroup(pg.df[pg.df[,"spectrum"] %in% 
                                        fData(x)[sel.spectra,"spectrum"],] )


  # remove from featureData
  featureData(x) <- as(fData(x)[sel.spectra,],"AnnotatedDataFrame")

  x
}






### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show",signature(object="IBSpectra"),
    function(object){
      p <- proteinGroup(object)
      callNextMethod(object)
      cat(nrow(fData(object))," spectra\n")
      show(p)
    }
)


###
### Logging
###

setGeneric("is.logged",function(x,name) standardGeneric("is.logged"))
setMethod("is.logged",signature("IBSpectra","character"), function(x,name) {
    name %in% rownames(x@log)
})

setGeneric("get.log",function(x,name) standardGeneric("get.log"))
setMethod("get.log",signature("IBSpectra","character"), function(x,name) {
  paste("\t",paste(apply(x@log[rownames(x@log) %in% name,,drop=FALSE],
                   1,paste,collapse="\t"),
        collapse="\n\t",sep=""),"\n",sep="")
})

setGeneric("do.log", function(x,name,msg) standardGeneric("do.log"))
setMethod("do.log",signature("IBSpectra","character","ANY"),function(x,name,msg) {
    if (is.logged(x,name)) {
        warning("operation ",name," already logged ",
                sum(rownames(x@log) == name)," times:\n",get.log(x,name))
    }
    if (is.logical(msg)) {
        if (msg) { msg = "TRUE"
        } else   { msg = "FALSE" }
    }
    tmp.matrix <- matrix(c(format(Sys.time(), "%F %T %Z"),msg),
                         ncol=2,dimnames=list(name))
                         
    message("LOG: ",sprintf("%10s: %s",name,msg))
    x@log <- rbind(x@log,tmp.matrix)
    return(x)
})

