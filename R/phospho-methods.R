
getPhosphoRSProbabilities <- function(
  id.file,mgf.file,massTolerance,activationType,simplify=FALSE,
  mapping.file=NULL,mapping=c(peaklist="even",id="odd"),pepmodif.sep="##.##",besthit.only=TRUE,
  phosphors.jar=system.file("phosphors","phosphoRS.jar",package="isobar",mustWork=TRUE)) {
  
  tmpfile <- tempfile("phosphors.")
  infile <- paste(tmpfile,".in.xml",sep="")
  outfile <- paste(tmpfile,".out.xml",sep="")

  writePhosphoRSInput(infile,
                      id.file,mgf.file,massTolerance,activationType,
                      mapping.file,mapping,pepmodif.sep)
  
  system(paste("java -jar",phosphors.jar,infile,outfile))
  readPhosphoRSOutput(outfile,simplify=simplify,pepmodif.sep,besthit.only=besthit.only)
}
  
writePhosphoRSInput <- function(phosphoRS.infile,id.file,mgf.file,massTolerance,activationType,
                                mapping.file=NULL,mapping=c(peaklist="even",id="odd"),pepmodif.sep="##.##") {

  if (is.data.frame(id.file)) 
    ids <- id.file
  else
    ids <- isobar:::.read.idfile(id.file,id.format="ibspectra.csv",log=NULL)

  ids <- unique(ids[,c("peptide","modif","spectrum")])
  # data[,SC['PEPTIDE']] <- gsub("I","L",data[,SC['PEPTIDE']])

  
  if (!is.null(mapping.file)) {
    mapping.quant2id <- do.call(rbind,lapply(mapping.file,function(f) {
      read.table(f,header=TRUE,sep=',',stringsAsFactors=FALSE)
    }))
    cn <-  colnames(mapping.quant2id)
    if (!all(mapping %in% cn))
      stop("mapping not correct")
    
    colnames(mapping.quant2id)[cn == mapping['id']] <- 'id'
    colnames(mapping.quant2id)[cn == mapping['peaklist']] <- 'peaklist'  
  }
  
  con.out <- file(phosphoRS.infile,'w')
  cat.f <- function(...) cat(...,"\n",file=con.out,sep="",append=TRUE)
  cat.f("<phosphoRSInput>")
  cat.f("  <MassTolerance Value='",0.02,"' />")
  cat.f("  <Phosphorylation Symbol='1' />")
  cat.f("  <Spectra>")

  input <- c()
  for (f in mgf.file) {
    con <- file(f,'r')
    input <- c(input,readLines(con))
    close(con)
  }

  modifs <- 
    rbind(c("PHOS",        "1","1:Phospho:Phospho:79.966331:PhosphoLoss:97.976896:STY"),
          c("Oxidation_M", "2","2:Oxidation:Oxidation:15.994919:null:0:M"),
          c("Cys_CAM",     "3","3:Carbamidomethylation:Carbamidomethylation:57.021464:null:0:C"), # check Cys_CAM mass!
          c("iTRAQ4plex","4","4:iTRAQ:iTRAQ:144.1544:null:0:KX")) #average mass? or monoisotopic?

  begin_ions <- which(input=="BEGIN IONS")+1
  end_ions <- which(input=="END IONS")-1
  titles <- gsub("TITLE=","",grep("TITLE",input,value=TRUE),fixed=TRUE)
  if (!all(ids$spectrum %in% titles))
    stop("Not all id spectrum titles are in MGF titles!\n",
         sum.bool.c(ids$spectrum %in% titles))
  
  if (length(begin_ions) != length(end_ions))
    stop("mgf file is errorneous, non-matching number",
         " of BEGIN IONS and END IONS tags");

  modif <- "PHOS"
  ids$modifrs <- .convertModifToPhosphoRS(ids$modif,modifs)

  pepid <- 0
  for (title in unique(ids[grep(modif,ids$modif),"spectrum"])) {

    spectrum_i <- which(titles==title)
    spectrum <- input[begin_ions[spectrum_i]:end_ions[spectrum_i]]
    ## read header
    header <- .strsplit_vector(spectrum[grep("^[A-Z]",spectrum)],"=")
    
    cat.f("    <Spectrum ID='",URLencode(header["TITLE"],reserved=TRUE),"'",
          " PrecursorCharge='",sub("+","",header["CHARGE"],fixed=TRUE),"'",
          " ActivationTypes='",activationType,"'>")
    
    peaks <- gsub(" ?","",spectrum[grep("^[0-9]",spectrum)],fixed=TRUE)
    cat.f("    <Peaks>",paste(gsub("\\s+",":",peaks),collapse=","),"</Peaks>")

    for (id_i in which(ids$spectrum==title)) {
      pepid <- pepid + 1
      cat.f("      <IdentifiedPhosphorPeptides>")
      cat.f("        <Peptide ID='",ids[id_i,"peptide"],pepmodif.sep,ids[id_i,"modif"],"'",
            " Sequence='",ids[id_i,"peptide"],"'",
            " ModificationInfo='",ids[id_i,"modifrs"],"' />")
      cat.f("      </IdentifiedPhosphorPeptides>")
    }
    cat.f("    </Spectrum>")
  }

  cat.f("  </Spectra>")
  cat.f("  <ModificationInfos>")
  for (i in seq_len(nrow(modifs))) {
    cat.f("    <ModificationInfo Symbol='",modifs[i,2],"' Value='",modifs[i,3],"' />")
  }
  cat.f("  </ModificationInfos>")
  cat.f("</phosphoRSInput>")
  close(con.out)
}

.convertModifToPhosphoRS <- function(modifstring,modifs) {
  sapply(strsplit(paste(modifstring," ",sep=""),":"),function(x) {
    x[length(x)] <- sub(" $","",x[length(x)])
    x[x==""] <- 0
    for (i in seq_len(nrow(modifs))) 
      x[grep(modifs[i,1],x,fixed=TRUE)] <- modifs[i,2]

    y <- c(x[1],".",x[2:(length(x)-1)],".",x[length(x)]);
    paste(y,collapse="") })
}

.convertModifToPos <- function(modifstring,modif="PHOS") {
  sapply(strsplit(paste(modifstring," ",sep=""),":"),function(x) {
    x[length(x)] <- sub(" $","",x[length(x)])
    paste(which(x==modif)-1,collapse="&")
  })
}

.convertPeptideModif <- function(peptide,modifstring,modifs=c("PHOS","Oxidation_M","Cys_CAM")) {
  names(letters) <- LETTERS
  mapply(function(pep,m) {
           m <- m[-c(1,length(m))]
           for (mm in modifs)
             pep[m==mm] <- letters[pep[m==mm]]
           paste(pep,collapse="")
         },
         strsplit(peptide,""),
         strsplit(paste(modifstring," ",sep=""),":")
)
}

readPhosphoRSOutput <- function(phosphoRS.outfile,simplify=FALSE,pepmodif.sep="##.##",
                                besthit.only=TRUE) {
  require(XML)
  doc <- xmlTreeParse(phosphoRS.outfile,useInternalNodes=TRUE)  
  spectra <- xmlRoot(doc)[["Spectra"]]
  res <- xmlApply(spectra,function(spectrum) {
    spectrum.id <- URLdecode(xmlAttrs(spectrum)["ID"])
    #message(spectrum.id)
    res.s <- xmlApply(spectrum[["Peptides"]],function(peptide) {
      pep.id <- strsplit(xmlAttrs(peptide)["ID"],pepmodif.sep,fixed=TRUE)[[1]]    
      #message(pep.id[1])
      site.probs <- t(xmlSApply(peptide[["SitePrediction"]],xmlAttrs))
      isoforms <- t(xmlSApply(peptide[["Isoforms"]],function(isoform) {
        seqpos <- xmlSApply(isoform[["PhosphoSites"]],xmlGetAttr,"SeqPos")

        # get right modif string
        modifstring <- strsplit(paste(pep.id[2]," ",sep=""),":")[[1]]
        modifstring <- gsub(" $","",modifstring)
        modifstring[modifstring=='PHOS'] <- ''
        modifstring[as.numeric(seqpos)+1] <- 'PHOS'
        modifstring <- paste(modifstring,collapse=":")
 
        if (length(seqpos > 1)) seqpos <- paste(seqpos,collapse="&")
        c(modif=modifstring,
          pepscore=xmlAttrs(isoform)[['PepScore']],
          pepprob=xmlAttrs(isoform)[['PepProb']],
          seqpos=seqpos)
      }))   

      #rownames(isoforms) <- NULL
      if (isTRUE(simplify))
        data.frame(peptide=pep.id[1],isoforms,
                   site.probs=paste(apply(round(site.probs,4),1,paste,collapse=":"),collapse=";"),
                   stringsAsFactors=FALSE,row.names=NULL)
      else
        list(peptide=pep.id,
             site.probs=site.probs,
             isoforms=isoforms)
    })
    if (simplify)
      data.frame(spectrum=spectrum.id,do.call(rbind,res.s),
                 stringsAsFactors=FALSE,row.names=NULL)
    else
      res.s
  })
  if (simplify) {
    res <- do.call(rbind,res)
    res$pepscore <- as.numeric(res$pepscore)
    res$pepprob <- as.numeric(res$pepprob)
    rownames(res) <- NULL
  } else {
    names(res) <- sapply(xmlChildren(spectra),xmlGetAttr,"ID")
  }
  if(besthit.only) {
    res <- ddply(res,'spectrum',function(d) d[which.max(d$pepprob),])
    rownames(res) <- res$spectrum
  }
  res
}

annotateSpectraPhosphoRS <- function(data,peaklist.file,min.prob=NULL,...) {
  if (is(data,"character"))
    data <- .read.idfile(data)
  probs <- getPhosphoRSProbabilities(data,peaklist.file,...,simplify=TRUE)
  ## probs excludes non-PHOS peptides - we do filter them for now? (about 8-10%)
  data$peptide <- NULL
  data$modif <- NULL
  data <- merge(data,probs,by="spectrum")
  if (!is.null(min.prob)) {
    if (!'use.for.quant' %in% colnames(data)) data$use.for.quant <- TRUE
    data[,"use.for.quant"] <-
      data[,"use.for.quant"] & data[,"pepprob"] >= min.prob
  }
  return(data)
}
