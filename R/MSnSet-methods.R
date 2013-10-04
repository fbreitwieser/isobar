setAs("MSnSet","IBSpectra",function(from) {
  ions <- exprs(from)
  q <- qual(from)
  lowerMz <- matrix(q$lowerMz,byrow=TRUE,ncol=ncol(ions))
  upperMz <- matrix(q$upperMz,byrow=TRUE,ncol=ncol(ions))
  maxInt <- matrix(q$maxInt,byrow=TRUE,ncol=ncol(ions))
  
  data <- fData(from)[,c("ProteinAccession","PeptideSequence","charge","precursor.mz","retention.time","spectrum")]
  colnames(data) <- c("accession","peptide","charge","exp.mass","retention.time","spectrum")

  my.class <- switch(colnames(ions)[1],
                     iTRAQ4.114 = "iTRAQ4plexSpectra",
                     iTRAQ8.113 = "iTRAQ8plexSpectra",
                     TMT2.126 = "TMT2plexSpectra",
                     TMT6.126 = "TMT6plexSpectra",
                     TMT10.126 = "TMT10plexSpectra",
                     "unknown")
  if (identical(my.class,"unknown"))
    stop("I do not know how to map MSnSet w/ columns [",paste(colnames(ions),collapse=","),"] to an IBSpectra object.")
         
  o <- new(my.class)
  rownames(ions) <- data$spectrum
  colnames(ions) <- o@reporterTagNames
  mass <- 0.5*(lowerMz+upperMz)
  rownames(mass) <- data$spectrum
  colnames(mass) <- o@reporterTagNames
  assayDataElements <- list(lowerMz=lowerMz,upperMz=upperMz,maxInt=maxInt)
  for (elem in names(assayDataElements)) {
    rownames(assayDataElements[[elem]]) <- data$spectrum
    colnames(assayDataElements[[elem]]) <- o@reporterTagNames
  }

  ib <- new(my.class,data,ions,mass,
            assayDataElements=assayDataElements)
  if (validObject(ib))
    return(ib)
})

setAs("IBSpectra","MSnSet",function(from) {
  ## based on quantify.MSnExp from MSnbase
  from <- ibspiked_set1
  
  elems <- assayDataElementNames(from)
  exprs <- reporterIntensities(from)

  get.elem <- function(x,y) {
    if (is.element(x,elems)) {
      return(assayDataElement(from,x))
    } else {
      if (is.matrix(y)) return(y)
      else return(assayDataElement(from,y))
    }
  }

  if (is(from,"iTRAQ4plexSpectra")) {
    reporters <- iTRAQ4
  } else if (is(from,"iTRAQ8plexSpectra")) {
    reporters <- iTRAQ8
  } else if (is(from,"TMT6plexSpectra")) {
    reporters <- TMT6
  } else if (is(from,"TMT10plexSpectra")) {
    reporters <- TMT10
  } else {
    stop("Cannot convert object")
  }
  
  lowerMz <- get.elem("lowerMz","mass")
  upperMz <- get.elem("upperMz","mass")
  maxInt  <- get.elem("maxInt",matrix(NA,nrow=nrow(lowerMz),ncol=ncol(lowerMz)))
  nMaxInt  <- get.elem("nMaxInt",matrix(1,nrow=nrow(lowerMz),ncol=ncol(lowerMz)))
  baseLength  <- get.elem("baseLength",matrix(NA,nrow=nrow(lowerMz),ncol=ncol(lowerMz)))

  .qual <- data.frame(maxInt=as.numeric(t(maxInt)),
                  nMaxInt=as.numeric(t(nMaxInt)),
                  baseLength=as.numeric(t(baseLength)),
                  lowerMz=as.numeric(t(lowerMz)),
                  upperMz=as.numeric(t(upperMz)),
                  reporter=mz(o),
                  precursor=rep(fData(from)$exp.mass,each=length(mz(o))))

  .exprs <- reporterIntensities(from)
  colnames(.exprs) <- reporters@reporterTagNames
  mapping <- c(spectrum="spectrum",
               ProteinAccession="accessions",
               ProteinDescription="NA",
               PeptideSequence="peptide",
               index="NA",
               file="NA",
               retention.time="retention.time",
               precursior.mz="exp.mass",
               peaks.count="NA",
               tic="NA",
               ms.level="NA",
               charge="charge",
               collision.energy="NA")
  df <- ibSpectra.as.concise.data.frame(from)
  df[,"NA"] <- NA
  
  fd <- df[,mapping]
  rownames(fd) <- fd[,"spectrum"]
  colnames(fd) <- names(mapping)
  .featureData <- new("AnnotatedDataFrame",data=fd[rownames(.exprs),])
  .phenoData <- new("AnnotatedDataFrame",
                    data=data.frame(mz=reporters@mz,
                      reporters=reporters@name,
                      row.names=reporters@reporterTagNames))
  
  msnset <- new("MSnSet",
                qual=.qual,
                exprs=.exprs,
                phenoData=.phenoData,
                featureData=.featureData,
                processingData=new("MSnProcess",processing=paste("created from IBSpectra object:",date())),
                annotation="No annotation")
  
  if (validObject(msnset))
    return(msnset)
})


