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
             " has different length then reporterTagNames [",n,"]!")
    if (!all(dim(object@isotopeImpurities) == n))
        stop("[IBSpectra: validation] isotopeImpurities has dimensions ",
             paste(dim(object@isotopeImpurities),collapse="x"),"!",
             " Expected is ",n,"x",n,".")
    
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

.ID.COLS <- c(SEARCHENGINE="search.engine",SCORE="score",SCORE.MASCOT="score.mascot",SCORE.PHENYX="score.phenyx")

.PTM.COLS <- c(SCORE.PHOSPHORS='pepscore',PROB.PHOSPHORS='pepprob',SEQPOS='seqpos',SITEPROBS='site.probs',PHOSPHO.SITES='phospho.sites')

.SPECTRUM.COLS <- c(PEPTIDE="peptide",MODIFSTRING="modif",CHARGE="charge",
                   THEOMASS="theo.mass",EXPMASS="exp.mass",
                   PARENTINTENS="parent.intens",RT="retention.time",
                   SPECTRUM="spectrum",SPECTRUM.QUANT="spectrum.quant",
                   .ID.COLS,USEFORQUANT="use.for.quant",
                   .PTM.COLS,
                   DISSOCMETHOD="dissoc.method",
                   PRECURSOR.PURITY="precursor.purity",
                   SCANS="scans",SCANS.FROM="scans.from",SCANS.TO="scans.to",
                   RAWFILE="raw.file",NMC="nmc",DELTASCORE="deltascore",
                   MASSDELTA.ABS="massdelta.abs",MASSDELTA.PPM="massdelta.ppm",
                   SAMPLE="sample",FILE="file",NOTES="notes")

.PEPTIDE.COLS <- c(PROTEINAC="accession",STARTPOS="start.pos",
                  REALPEPTIDE="real.peptide")

.PROTEIN.COLS <- c(PROTEINAC="accession",PROTEINAC_CONCISE="accessions",
                   NAME="name",PROTEIN_NAME="protein_name",
                   GENE_NAME="gene_name",ORGANISM="organism")

.ROCKERBOX.COLS <- c(PROTEINAC="accession",STARTPOS="start.pos",MODIFSTRING="modif",
                     QN="mascot.query.number",RANK="rank",SCANS="scan.number.s.",RT="retention.time",
                     RAWFILE="raw.file",PEPTIDE="sequence", "phosphosequence",NMC="miscleavages",
                     SCORE="score",DELTASCORE="deltascore",PHOSPHO.DELTASCORE="phospho_deltascore",
                     EXPMASS="peptide.mr",MASSDELTA.ABS="mass.delta..abs.",MASSDELTA.PPM="mass.delta..ppm.",
                     ROCKERBOX_MOD="modifications",PROTEINACS="all.protein.matches",SPECTRUM="scan.title")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor and readIBSpectra.
###

setMethod("initialize","IBSpectra",
    function(.Object,data=NULL,data.ions=NULL,data.mass=NULL,
             proteinGroupTemplate=NULL,fragment.precision=NULL,
             assayDataElements=list(),allow.missing.columns=FALSE,
             annotate.spectra.f=NULL,merge.ids=TRUE,
             write.excluded.to=NULL,...) { 
      if (is.null(data)){
        callNextMethod(.Object,...)
      } else {

        reporterTagNames <- reporterTagNames(.Object)
        data <- .factor.to.chr(data)
        if (is.function(annotate.spectra.f)) {
          data <- annotate.spectra.f(data)
        } else {
          if (!is.null (annotate.spectra.f))
            stop("annotate.spectra.f should be a function!")
        }

        SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(data)]
        PC <- .PROTEIN.COLS[.PROTEIN.COLS %in% colnames(data)]

        # check that obligatory columns are present
        missing.cols = c()
        for (col in c('SPECTRUM','PEPTIDE','MODIFSTRING'))
          if (!is.element(.SPECTRUM.COLS[col],SC))
            missing.cols <- c(missing.cols,.SPECTRUM.COLS[col])
        for (col in c('PROTEINAC'))
          if (!is.element(.PROTEIN.COLS[col],PC))
            missing.cols <- c(missing.cols,.PROTEIN.COLS[col])
        if (is.null(data.ions) || is.null(data.mass)) {
          reagentfields <- c(sprintf(.PEAKS.COLS['MASSFIELD'],reporterTagNames),
                             sprintf(.PEAKS.COLS['IONSFIELD'],reporterTagNames))
          if (!all(reagentfields %in% colnames(data))) {
            missing.cols <- c(missing.cols,
                              reagentfields[!reagentfields %in% colnames(data)])
          }
        }

        # handle missing columns
        if (length(missing.cols) > 0) {
            msg <- paste("not all required columns in data, the following are missing: \n\t",
                     paste(missing.cols,collapse="\n\t"))
            if (allow.missing.columns) {
                warning(msg)
                data[,missing.cols] <- NA
            } else {
                stop(msg)
            }
        }
        
        message("merge ids")
        data <- .merge.identifications(data)
        SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(data)]
        
        data.sc <- unique(data[,SC])
        # Check for divergent identifications of spectra
        if (anyDuplicated(data.sc[,SC['SPECTRUM']])) {
          dupl.spectra <- .get.dupl.n.warn(data.sc,SC['SPECTRUM'],write.to=write.excluded.to)
          #dupl.spectra <- .get.dupl.n.warn(data.sc[,SC[c("SPECTRUM","PEPTIDE","MODIFSTRING")]],SC['SPECTRUM'])
          # problem: different values in any column - what info should you choose??
          data <- data[!data[,SC['SPECTRUM']] %in% dupl.spectra,]
        }
        
        # create ProteinGroup
        PROTEINGROUP.COLS <- SC[c('SPECTRUM','PEPTIDE','MODIFSTRING')]
        if (.PEPTIDE.COLS['STARTPOS'] %in% colnames(data)) 
          PROTEINGROUP.COLS <- c(PROTEINGROUP.COLS,.PEPTIDE.COLS['STARTPOS'])
        PROTEINGROUP.COLS <- c(PROTEINGROUP.COLS,.PROTEIN.COLS['PROTEINAC'])

        proteinGroup <-
          ProteinGroup(unique(data[,PROTEINGROUP.COLS]),
                       template=proteinGroupTemplate)
        
        if (is.null(assayDataElements))
          assayDataElements <- list()
          
        # featureData: only keep 'spectrum columns'
        # assayData: create mass and ions matrices
        if (is.null(data.ions) || is.null(data.mass)) {
          reagentfields.mass <- sprintf(.PEAKS.COLS['MASSFIELD'],reporterTagNames)
          reagentfields.ions <- sprintf(.PEAKS.COLS['IONSFIELD'],reporterTagNames)
          data <- unique(data[,c(SC,reagentfields.mass,reagentfields.ions)])
          rownames(data) <- data[,SC['SPECTRUM']]
          assayDataElements$mass <- as.matrix(data[,reagentfields.mass])
          assayDataElements$ions <- as.matrix(data[,reagentfields.ions])
          dimnames(assayDataElements$mass) <- list(data$spectrum,reporterTagNames)
          dimnames(assayDataElements$ions) <- list(data$spectrum,reporterTagNames)
        } else {
          data <- unique(data[,SC])
          rownames(data) <- data[,SC['SPECTRUM']]
          
          if (is.null(rownames(data.ions)) || is.null(rownames(data.mass)))
            stop("provide spectrum ids as rownames for data.ions and data.mass")
          if (!identical(rownames(data.ions),rownames(data.mass)))
            stop("data.ions and data.mass do not have identical rownames!")
          for (elem in assayDataElements) {
            if (is.null(rownames(elem)))
              stop("provide rownames for all assayDataElements")
            if (!identical(rownames(elem),rownames(data.ions)))
              stop("all assayDataElements must have the same rownames")
            if (!identical(dim(data.ions),dim(elem)))
              stop("all assayDataElemetns must have the same dim")
          }
          colnames(data.ions) <- reporterTagNames
          colnames(data.mass) <- reporterTagNames
          assayDataElements$ions <- data.ions
          assayDataElements$mass <- data.mass
          
          if (!identical(rownames(data),rownames(data.ions))) {
            id.n.quant <- intersect(rownames(data),rownames(data.mass))
            id.not.quant <- setdiff(rownames(data),rownames(data.mass))
            quant.not.id <- setdiff(rownames(data.mass),rownames(data))

            if (length(id.n.quant)==0) stop("No spectra could be matched between identification and quantifications")
            
            if (length(quant.not.id) > 0)
              warning(length(quant.not.id)," spectra could not be mapped from quant to id: ",
                      paste(quant.not.id[1:min(length(quant.not.id),10)],collapse=";"))
            if (length(id.not.quant) > 0)
              warning(length(id.not.quant)," spectra could not be mapped from id to quant: ",
                      paste(id.not.quant[1:min(length(id.not.quant),10)],collapse=";"))

            na.matrix <- matrix(NA,nrow=nrow(data),ncol=ncol(data.mass),
                                dimnames=list(rownames(data),colnames(data.mass)))

            for (elem in names(assayDataElements)) {
              tmp <- na.matrix
              tmp[id.n.quant,] <- assayDataElements[[elem]][id.n.quant,]
              assayDataElements[[elem]] <- tmp
            }
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

        if (!.SPECTRUM.COLS['USEFORQUANT'] %in% colnames(data)) {
          # not perfect yet - better: set spectra of peptides shared between groups to FALSE
          #                            and spectra with no values
          data[,.SPECTRUM.COLS['USEFORQUANT']] <- TRUE
          #data[,.SPECTRUM.COLS['USEFORQUANT']] <- !apply(is.na(assayDataElements$ions),1,all)
        }

        nn <- .SPECTRUM.COLS %in% colnames(data)
        fdata <- data[,colnames(data) %in% .SPECTRUM.COLS]

        label.desc <- c(PEPTIDE='peptide sequence',
            MODIFSTRING='modifications of peptide',
            CHARGE='peptide charge state',
            THEOMASS='theoretical peptide mass',
            EXPMASS='experimental peptide mass',
            PARENTINTENS='parent ion intensity',
            RT='retention time',
            DISSOCMETHOD='dissociation METHOD',
            SPECTRUM='spectrum title',
            SPECTRUM.QUANT='title of spectrum used for quantiation',
            PRECURSOR.PURITY="precursor purity",
            SCANS.FROM="scans from",SCANS.TO="scans to",
            RAWFILE="raw file",NMC="nmc",DELTASCORE="deltascore",
            SCANS="scans",MASSDELTA.ABS="massdelta (abs)",MASSDELTA.PPM="massdelta (ppm)",
            SEARCHENGINE='protein search engine',
            SCORE='protein search engine score',
            SCORE.MASCOT="Mascot search score",
            SCORE.PHENYX="Phenyx search score",
            USEFORQUANT='use spectrum for quantification',
            SCORE.PHOSPHORS='PhosphoRS pepscore',
            PROB.PHOSPHORS="PhosphoRS probability",
            PHOSPHO.SITES="phosphorylation sites",
            PEPPROB='PhosphoRS pepprob',
            SEQPOS='PTM seqpos',
            SITEPROBS='PhosphoRS site.probs',
            FILE='file',SAMPLE='sample',NOTES='notes'
            )
        if (!all(names(.SPECTRUM.COLS) %in% names(label.desc)))
          stop("Not all SPECTRUM COLS have a label description:\n\t",
               paste(names(.SPECTRUM.COLS)[!names(.SPECTRUM.COLS) %in% names(label.desc)],collapse="\n\t"))

        VARMETADATA=data.frame(labelDescription=label.desc[names(.SPECTRUM.COLS)],
                               row.names=.SPECTRUM.COLS)
	    
        featureData <- new("AnnotatedDataFrame",data=fdata,
            varMetadata=VARMETADATA[nn,,drop=FALSE])
	   
        assayData=do.call(assayDataNew, assayDataElements, envir=parent.frame())
        
        ## create ibspectra object
        callNextMethod(.Object,
            assayData=assayData,
            featureData=featureData,
            proteinGroup=proteinGroup,
            reporterTagNames=reporterTagNames,...)
      }
    }
)

setGeneric("readIBSpectra", function(type,id.file,peaklist.file,...)
           standardGeneric("readIBSpectra"))
setMethod("readIBSpectra",
          signature(type="character",id.file="character",peaklist.file="missing"),
    function(type,id.file,...) {
      ll <- lapply(seq_along(id.file),function(i) {
                   df <- read.table(id.file[i],header=T,sep="\t")
                   df[,.SPECTRUM.COLS['FILE']]  <- id.file[i]
                   if (!is.null(names(id.file)))
                     df[,.SPECTRUM.COLS['SAMPLE']]  <- names(id.file)[i]
                   df
            })
      new(type,data=do.call(rbind,ll),...)
    }
)

## TODO: log is not returned
.read.idfile <- function(id.file,id.format=NULL,decode.titles=TRUE,log=NULL) {
  data <- unique(do.call("rbind",lapply(id.file,function(f) {
    
    if (is.null(id.format)) {
      if (grepl(".mzid$",f,ignore.case=TRUE)) 
        id.format.f <- "mzid"
      else if (grepl(".ibspectra.csv$",f,ignore.case=TRUE) ||
               grepl(".id.csv$",f,ignore.case=TRUE))
        id.format.f <- "ibspectra.csv"
      else if (grepl(".peptides.csv$",f) || grepl(".peptides.txt$",f)) 
        id.format.f <- "rockerbox"
      else
        stop(paste("cannot parse file ",f," - cannot deduce format based on extenstion (it is not ibspectra.csv, id.csv, peptides.txt or mzid). Please provide id.format to readIBSpectra",sep=""))
    } else {
      id.format.f <- id.format
    }
    
    if (id.format.f == "ibspectra.csv") {
      data <- read.table(f,header=T,stringsAsFactors=F,sep="\t")
      log <- rbind(log,c("identification file [id.csv]",f))
    } else if (id.format.f == "mzid") {
      data <- read.mzid(f)
      log <- rbind(log,c("identification file [mzid]",f))
    } else if (id.format.f == "rockerbox") {
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

.merge.identifications <- function(identifications) {
  SC <- .SPECTRUM.COLS[.SPECTRUM.COLS %in% colnames(identifications)]
  PC <- .PROTEIN.COLS[.PROTEIN.COLS %in% colnames(identifications)]
  ## Substitute Isoleucins with Leucins (indistinguishable by Masspec)
  identifications[,.PEPTIDE.COLS['REALPEPTIDE']] <- identifications[,SC['PEPTIDE']]
  identifications[,SC['PEPTIDE']] <- gsub("I","L",identifications[,SC['PEPTIDE']])
  protein.colnames <- which(colnames(identifications) %in% c(.PROTEIN.COLS['PROTEINAC'],.PEPTIDE.COLS['STARTPOS'],SC['PEPTIDE']))
  pept.n.prot <- unique(identifications[,protein.colnames])
  identifications <- unique(identifications[,-which(colnames(identifications) %in% c(.PROTEIN.COLS['PROTEINAC'],.PEPTIDE.COLS['STARTPOS']))])

  ## Merge results of different search engines / on different spectra
  if (.SPECTRUM.COLS['SEARCHENGINE'] %in% colnames(identifications)) {
    if (any(grepl("|",identifications[,SC['SEARCHENGINE']],fixed=TRUE))) {
      ## scores are merged together
      engines <- strsplit(identifications[,SC['SEARCHENGINE']],"|",fixed=TRUE)
      scores <- strsplit(identifications[,SC['SCORE']],"|",fixed=TRUE)
      for (engine in unique(unlist(engines))) {
        name <- paste0('score.',tolower(engine))
        e.scores <- mapply(function(e,s) if(any(e==engine)) as.numeric(s[e==engine]) else NA,engines,scores)
        identifications[,name] <- e.scores
      }
      identifications[,SC['SEARCHENGINE']] <- NULL
      identifications[,SC['SCORE']] <- NULL
    }
    ## keep only spectra which have same identification with both engines
  }

  score.colname <- .SPECTRUM.COLS[c('SCORE.PHOSPHORS','SCORE.MASCOT','SCORE.PHENYX')]
  score.colname <- score.colname[score.colname %in% colnames(identifications)]
  if (.SPECTRUM.COLS['SPECTRUM.QUANT'] %in% colnames(identifications)) {
    tt <- table(identifications[,SC['SPECTRUM.QUANT']])
    if (any(tt>1)) {
      identifications <- ddply(identifications,.SPECTRUM.COLS['SPECTRUM.QUANT'],function(x) {
                               if (nrow(x) == 1) return(x)
                               my.args <- as.list(x[,score.colname])
                               my.args$decreasing=TRUE
                               max.hit <- do.call(order,my.args)[1]
                               if (!all(x[,SC['PEPTIDE']] == x[max.hit,SC['PEPTIDE']]))
                                 return(NULL)
                               if (all(x[,SC['MODIFSTRING']] == x[max.hit,SC['MODIFSTRING']]) && 'DISSOCMETHOD' %in% names(SC))
                                 x[max.hit,SC['DISSOCMETHOD']] <- paste(sort(unique(x[,SC['DISSOCMETHOD']])),collapse="&")
                               return(x[max.hit,])
      })
    }
  }

  tt <- table(identifications[,SC['SPECTRUM']])
  if (any(tt>1)) {
    identifications <- ddply(identifications,.SPECTRUM.COLS['SPECTRUM'],function(x) {
          if (nrow(x) == 1) return(x)
          my.args <- as.list(x[,score.colname])
          my.args$decreasing=TRUE
          max.hit <- do.call(order,my.args)[1]
          if (!all(x[,SC['PEPTIDE']] == x[max.hit,SC['PEPTIDE']]))
            return(NULL)
          if (all(x[,SC['MODIFSTRING']] == x[max.hit,SC['MODIFSTRING']])  && 'DISSOCMETHOD' %in% names(SC))
            x[max.hit,SC['DISSOCMETHOD']] <- paste(sort(unique(x[,SC['DISSOCMETHOD']])),collapse="&")
          return(x[max.hit,])
    })
  }
  return(merge(pept.n.prot,identifications,by="peptide",all.y=TRUE))
}

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
  return(.merge.identifications(identifications))
}



writeIBSpectra <- function(ibspectra,file,sep="\t",row.names=FALSE,...) {
  write.table(as.data.frame(ibspectra),file=file,sep=sep,row.names=row.names,...)
  message("finished writing ibspectra to ",file)
}

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
          signature(type="character",id.file="character",peaklist.file="character"),
    function(type,id.file,peaklist.file,
             proteinGroupTemplate=NULL,
             mapping.file=NULL,mapping=c(peaklist="even",id="odd"),
#             mapping.file.readopts=list(header=TRUE,stringsAsFactors=FALSE,sep=","),
             id.file.domap=NULL,annotate.spectra.f=NULL,
             peaklist.format=NULL,id.format=NULL,fragment.precision=NULL,fragment.outlier.prob=NULL,
             decode.titles=TRUE,scan.lines=0,...) {
      
      log <- data.frame(key=c(),message=c())

      ## get identified spectra
      #data <- .read.idfile(id.file,id.format,decode.titles=decode.titles,log)
      data <- .read.identifications(id.file,mapping=mapping.file,mapping.names=mapping,
                                    identifications.quant=id.file.domap,
                                    identifications.format=id.format,decode.titles=decode.titles)

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
.read.mapping <- function(mapping.file,readopts,mapping) {
  ## STUB
  if (is.null(mapping.file)) return(NULL)

  if (!all(c("peaklist","id") %in% names(mapping))) {
    stop("readIBSpectra/mapping must be a named vector with the",
         " names peaklist and id.")
  }
  ## read mapping file(s)
  mapping.quant2id <- do.call(rbind,lapply(mapping.file,function(f) {
                                           readopts$file <- f
                                           do.call(read.table,readopts)
                }))

  cn <-  colnames(mapping.quant2id)
  colnames(mapping.quant2id)[c(mapping['id'],mapping['peaklist'])] <- c('id','peaklist')
  colnames(mapping.quant2id)[cn == mapping['peaklist']] <- 'peaklist'
  log <- rbind(log,data.frame(rep("mapping file",length(mapping.file)),mapping.file,stringsAsFactors=FALSE))
  if (!is.null(id.file.domap)) {
    data.m <- .read.idfile(id.file.domap,id.format,decode.titles=decode.titles,log)
    map.spectrum <- mapping.quant2id[,"id"]
    names(map.spectrum) <- mapping.quant2id[,"peaklist"]
    data.m[,"spectrum"] <- map.spectrum[data.m[,"spectrum"]]
    data <- rbind(data,data.m)
  }
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

.get.dupl.n.warn <- function(df,col,msg="ibspectra",write.to=NULL) {
  dupl <- .all.duplicate.rows(df,col)
  if (!is.null(write.to))
    write.table(dupl,file=write.to,row.names=FALSE,sep="\t")
  dupl.msg <- paste(apply(dupl,1,paste,collapse="; "),collapse="\n\t")
  warning(sprintf("%s> divergent identifications in %s spectra [%s ids]:\n\t%s",
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

.IBSpectraAsConciseDataFrame  <- function(from) {

      # prepare ProteinGroup data.frame
      pg.df <- as(proteinGroup(from),"data.frame.concise")
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
  pg.df <- as(proteinGroup(from),"data.frame.concise")
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
    function(x,element="ions",na.rm=FALSE,...) {
      sel <- spectrumSel(x,...)
      data <- assayDataElement(x,element)[sel,,drop=FALSE]

      if (na.rm & length(data) > 0) 
        return(data[apply(!is.na(data),1,all),,drop=FALSE])
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

setMethod("spectrumSel",signature(x="IBSpectra",peptide="matrix",protein="missing"),
    function(x,peptide,modif=NULL,spectrum.titles=FALSE,use.for.quant.only=TRUE,do.warn=TRUE) {
        if (length(peptide) == 0) {
          warning("0L peptide provided")
          return(FALSE)
        }
        if (ncol(peptide) != 2)
          stop("don't know how to handle matrix with ",ncol(peptide)," columns!")
        
        sel <- fData(x)[,.SPECTRUM.COLS['PEPTIDE']]  %in% peptide[,1] & 
               fData(x)[,.SPECTRUM.COLS['MODIFSTRING']]  %in% peptide[,2]

        if (use.for.quant.only && .SPECTRUM.COLS['USEFORQUANT'] %in% colnames(fData(x)))
          sel <- sel & fData(x)[,.SPECTRUM.COLS['USEFORQUANT']] 
        
        for (m in modif)
          sel <- sel & grepl(m,fData(x)[,.SPECTRUM.COLS['MODIFSTRING']])
        if (!any(sel) && do.warn) warning("No spectra for peptide ",peptide)
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
      if ((spectrum.titles & any(sel)) || (!spectrum.titles & !any(sel)))
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
      if (is.null(names(value))) {
        phenoData(x)[["class.labels",labelDescription="class labels"]] <- value
      } else {
        phenoData(x)[["class.labels",labelDescription="class labels"]] <- as.character(value)
        phenoData(x)[["class.label.description",labelDescription="class label descriptions"]] <- names(value)
      }
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
      ri.corrected[ri.corrected<0] <- NA
      reporterIntensities(x) <- ri.corrected

      x <- do.log(x,"isotopeImpurities.corrected",TRUE)

      return(x)
    }
)

normalize <- function(x,f=median,target="intensity",exclude.protein=NULL,
                      use.protein=NULL,f.doapply=TRUE,log=TRUE,
                      channels=NULL,na.rm=FALSE,per.file=TRUE,...){
  
  ## NOTE: median normalizes might normalize too much when a lot of NA
  ##         values are present - mean or colSums better?
  
  if (!is.logged(x,"isotopeImpurities.corrected")) 
    warning("Isotope impurity correction has not been logged",
            " - data might be uncorrected. See ?correctIsotopeImpurities.")

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


  if (!is.null(channels) ) {
    if (is.list(channels)) {
      for (channels.set in channels)
        x <- normalize(x,f=f,target=target,exclude.protein=exclude.protein,
                       use.protein=use.protein,f.doapply=f.doapply,
                       log=log,channels=channels.set,na.rm=na.rm,...)
      return(x)
    } else {
      if (!all(channels %in% colnames(ri)))
        stop("channels must be reporterTagNames.")
                                                         
      ri <- reporterIntensities(x,na.rm=na.rm)[,channels,drop=FALSE]
    }
  } else {
    ri <- reporterIntensities(x,na.rm=na.rm)
  } ## else is.null(channels)
  
  if (!is.null(exclude.protein))
    ri <- ri[!spectrumSel(x,protein=exclude.protein,
                          specificity=c(REPORTERSPECIFIC,GROUPSPECIFIC,UNSPECIFIC)),]
  
  if (!is.null(use.protein))
    ri <- ri[spectrumSel(x,protein=use.protein,specificity=REPORTERSPECIFIC),]
  
  ## TODO: warning when ri is empty

  x <- do.log(x,"is.normalized",TRUE)
  ## save original reporter intensities for noise estimation
  if (is.null(assayDataElement(x,"ions_not_normalized")))
    assayDataElement(x,"ions_not_normalized") <- reporterIntensities(x)

  if (per.file && .SPECTRUM.COLS['FILE'] %in% colnames(fData(x))) {
    fd <- fData(x)
    for (n.file in sort(unique(fd[,.SPECTRUM.COLS['FILE']]))) {
      message(n.file)
      sel <- fd[,.SPECTRUM.COLS['FILE']] == n.file
      message(sum(sel))
      factor <- .get.normalization.factors(ri[sel,,drop=FALSE],f,target,f.doapply,...)
      reporterIntensities(x)[sel,colnames(ri)] <- 
        reporterIntensities(x)[sel,colnames(ri)]*rep(factor,each=sum(sel))
      for (i in seq_along(factor)) {
        x <- do.log(x,paste("normalization.multiplicative.factor file",n.file,"channel",colnames(ri)[i]),
                    round(factor[i],4))
      } 
    }
  } else {
    factor <- .get.normalization.factors(ri,f,target,f.doapply,...)
    reporterIntensities(x)[,colnames(ri)] <- 
      reporterIntensities(x)[,colnames(ri)]*
      rep(factor,each=nrow(reporterIntensities(x)))
    for (i in seq_along(factor)) {
      x <- do.log(x,paste("normalization.multiplicative.factor channel",
                          colnames(ri)[i]),
                  round(factor[i],4))
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

  return(max(res,na.rm=T)/res)
 
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



subsetIBSpectra <- function(x, protein=NULL, peptide=NULL, direction="exclude",
                            specificity=c(REPORTERSPECIFIC,GROUPSPECIFIC,
                                             UNSPECIFIC),...) {

  if (is.null(protein))
    sel.spectra <- spectrumSel(x,peptide=peptide,...)
  else
    sel.spectra <- spectrumSel(x,protein=protein,specificity=specificity,...)
  if (direction=="exclude")
    sel.spectra <- !sel.spectra

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

