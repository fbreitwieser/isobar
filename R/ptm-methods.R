
getPhosphoRSProbabilities <- function(
  id.file,mgf.file,massTolerance,activationType,simplify=FALSE,
  mapping.file=NULL,mapping=c(peaklist="even",id="odd"),pepmodif.sep="##.##",besthit.only=TRUE,
  phosphors.cmd=paste("java -jar",system.file("phosphors","phosphoRS.jar",package="isobar",mustWork=TRUE)),
  file.basename=tempfile("phosphors.")) {

  infile <- paste0(file.basename,".in.xml")
  outfile <- paste0(file.basename,".out.xml")

  writePhosphoRSInput(infile,
                      id.file,mgf.file,massTolerance,activationType,
                      mapping.file,mapping,pepmodif.sep)
  
  system(paste(phosphors.cmd,infile,outfile))
  readPhosphoRSOutput(outfile,simplify=simplify,pepmodif.sep=pepmodif.sep,besthit.only=besthit.only)
}

.iTRAQ.mass = c(monoisotopic = 144.102063, average = 144.1544)
.iTRAQ8.mass = c(monoisotopic = 304.2, average = 304.308)
.CysCAM.mass = c(monoisotopic = 57.021464, average =  57.0513)
.OxidationM.mass = c(monoisotopic = 15.994915, average =  15.9994)
  
writePhosphoRSInput <- 
  function(phosphoRS.infile,id.file,mgf.file,massTolerance,activationType,
           mapping.file=NULL,mapping=c(peaklist="even",id="odd"),pepmodif.sep="##.##",
           modif.masses=
           rbind(c("PHOS",       "1","1:Phospho:Phospho:79.966331:PhosphoLoss:97.976896:STY"),
                 c("Oxidation_M","2","2:Oxidation:Oxidation:15.994919:null:0:M"),
                 c("Cys_CAM",    "3","3:Carbamidomethylation:Carbamidomethylation:57.021464:null:0:C"),
                 c("iTRAQ4plex", "4","4:iTRAQ4:iTRAQ4:144.1544:null:0:KX"),
                 c("iTRAQ8plex", "5","5:iTRAQ8:iTRAQ8:304.308:null:0:KX"))
) {

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
  ids$modifrs <- .convertModifToPhosphoRS(ids$modif,modif.masses)

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
  for (i in seq_len(nrow(modif.masses))) {
    cat.f("    <ModificationInfo Symbol='",modif.masses[i,2],"' Value='",modif.masses[i,3],"' />")
  }
  cat.f("  </ModificationInfos>")
  cat.f("</phosphoRSInput>")
  close(con.out)
}

calc.delta.score <- function(data) {
  pep.n.prot <- unique(data[,c("accession","peptide","start.pos")])
  data$accession <- NULL
  data$start.pos <- NULL
  data <- unique(data)
  if (!any(by(data$score,data$spectrum, length)>1)) {
    stop("Cannot calculate delta score: Only one hit per spectrum available")
  }

  data$delta.score <- data$score
  data$n.pep <- 1
  data$n.loc <- 1

  res <- ddply(data,"spectrum",function(x) {
    if (nrow(x) == 1) return(x);
    res <- x[which.max(x$score),,drop=FALSE]
    res$n.pep <- length(unique(x$peptide))
    x <- x[-which.max(x$score),] # remove best hit from x
    res$delta.score <- res[,'score'] - x[which.max(x$score),'score'] # calc delta score w/ max
    res$delta.score.pep <- res$delta.score
    x <- x[x$peptide == res[,'peptide',],] # only keep same peptide hits in x
    if (nrow(x) == 0) return(res);
    res$delta.score.pep <- res[,'score'] - x[which.max(x$score),'score'] # calc delta score w/ max (same pep)
    res$n.loc <- nrow(x) + 1
    return(res);
  })

  data <- merge(pep.n.prot,res,by="peptide",all.y=TRUE)

  return(data[order(data[,"accession"],data[,"peptide"]),])
}

filterSpectraDeltaScore <- function(data, min.delta.score=10, do.remove=FALSE) {
  if (!"delta.score" %in% colnames(data))
    data <- calc.delta.score(data)
  
  if (!is.null(min.delta.score)) {
    sel.mindeltascore <- data[,"delta.score"] >= min.delta.score
    data[,"use.for.quant"] <- data[,"use.for.quant"] & sel.mindeltascore
    if (isTRUE(do.remove))
      data <- data[,sel.mindeltascore]
  }
  return(data)
}


.convertModifToPhosphoRS <- function(modifstring,modifs) {
  sapply(strsplit(paste0(modifstring," "),":"),function(x) {
    x[length(x)] <- sub(" $","",x[length(x)])
    x[x==""] <- 0
    for (i in seq_len(nrow(modifs))) 
      x[grep(modifs[i,1],x,fixed=TRUE)] <- modifs[i,2]

    y <- c(x[1],".",x[2:(length(x)-1)],".",x[length(x)]);
    paste(y,collapse="") })
}

## TODO:
#.convertModifToPosAndName <- function(modifstring,modif=c(p="PHOS",o="Oxidation_M",c="Cys_CAM",
#                                          me="METH_KR",me="METH_K",me="METH_R",
#                                          me2="BIMETH_KR",me2="BIMETH_K",me2="BIMETH_R",
#                                          me3="TRIMETH_K",
#                                          ac="ACET_K"),collapse="&",simplify=TRUE) {
#  sapply(strsplit(paste0(modifstring," "),":"),function(x) {
#    x[length(x)] <- sub(" $","",x[length(x)])
#    res <- sapply(seq_along(x),function(x.i) {
#      which(modif==indiv.modif)
#    })
#    if (!is.null(collapse))
#      paste(which(x%in%modif)-1,collapse=collapse)
#    else
#      which(x%in%modif)-1
#  },simplify=simplify)
#}
#
.convertModifToPos <- function(modifstring,modif="PHOS",collapse="&",simplify=TRUE,and.name=FALSE) {
  split.modif <- strsplit(paste0(modifstring," "),":")
  if (and.name) {
    name.modif <- .names.as.vector(modif)
  }
  sapply(split.modif,function(x) {
    x[length(x)] <- sub(" $","",x[length(x)])
    modification.pos <- which(x%in%modif)-1 
    if (!is.null(collapse))
      paste(modification.pos,collapse=collapse)
    else {
      if (and.name) 
        data.frame(modif.pos=modification.pos,modif=name.modif[x[modification.pos+1]],stringsAsFactors=FALSE)
      else
        modification.pos
    }
  },simplify=simplify)
}



.convertPhosphoRSPepProb <- function(peptide,pepprob) {
  mapply(function(pep,pprob) {
           pprob <- as.numeric(pprob)
           prob <- rep(-1,length(pep))
           prob[pprob[seq(from=1,to=length(pprob),by=2)]] <- pprob[seq(from=2,to=length(pprob),by=2)]

           pep[prob>=0] <- paste0(pep[prob>=0],"(",prob[prob>=0],")")
           paste0(pep,collapse="")
         },
         strsplit(peptide,""),
         strsplit(pepprob,"[;:]"))
}

.convertPeptideModif <- function(peptide,modifstring,
                                 modifs=c(p="PHOS",o="Oxidation_M",c="Cys_CAM",
                                          me="METH_KR",me="METH_K",me="METH_R",
                                          dime="BIMETH_KR",dime="BIMETH_K",dime="BIMETH_R",
                                          trime="TRIMETH_K",
                                          ac="ACET_K")) {
  names(letters) <- LETTERS
  if (length(peptide)==0 || all(nchar(peptide)==0))
	  stop("peptide length=0")
  mapply(function(pep,m) {
           m <- m[-c(1,length(m))]
           if (is.null(names(modifs))) {
             for (mm in modifs) 
               pep[m==mm] <- letters[pep[m==mm]]
           } else {
             for (i in seq_along(modifs)) 
               pep[m==modifs[i]] <- paste0(pep[m==modifs[i]],"(",names(modifs)[i],")")
           }
           paste(pep,collapse="")
         },
         strsplit(peptide,""),
         strsplit(paste0(modifstring," "),":")
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
        modifstring <- strsplit(paste0(pep.id[2]," "),":")[[1]]
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
      storage.mode(site.probs) <- "numeric"
      site.probs[,2] <- round(site.probs[,2],2)
      if (isTRUE(simplify))
        data.frame(peptide=pep.id[1],isoforms,
                   site.probs=paste(apply(site.probs,1,paste,collapse=":"),collapse=";"),
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
  if(besthit.only & simplify) {
    res <- ddply(res,'spectrum',function(d) d[which.max(d$pepprob),])
    rownames(res) <- res$spectrum
  }
  res
}

filterSpectraPhosphoRS <- function(id.file,mgf.file,...,min.prob=NULL, do.remove=FALSE) {
  if (is(id.file,"character"))
    id.file <- .read.idfile(id.file)
  probs <- getPhosphoRSProbabilities(id.file,mgf.file,...,simplify=TRUE)
  ## probs excludes non-PHOS peptides - we do filter them for now? (about 8-10%)
  id.file$peptide <- NULL
  id.file$modif <- NULL
  id.file <- merge(id.file,probs,by="spectrum")
  if (!is.null(min.prob)) {
    if (!'use.for.quant' %in% colnames(id.file)) id.file$use.for.quant <- TRUE
    sel.minprob <- id.file[,"pepprob"] >= min.prob
    id.file[,"use.for.quant"] <- id.file[,"use.for.quant"] & sel.minprob
    if (isTRUE(do.remove))
      id.file <- id.file[,sel.minprob]
  }
  return(id.file)
}

modif.sites <- function(protein.group,protein.g=reporterProteins(protein.group),modif) {
  pi <- protein.group@peptideInfo
  ip <- indistinguishableProteins(protein.group)
  sapply(protein.g,
         function(my.protein.g) {
           isoform.acs <- names(ip)[ip==my.protein.g]
           res <- lapply(isoform.acs,function(ac) {
                         sel <- pi[,"protein"] == ac
                         pep.pos <- .convertModifToPos(pi[sel,"modif"],modif,simplify=FALSE,collapse=NULL) 
                         modif.pos <- unlist(mapply(function(start.pos,pep.posi) start.pos + pep.posi -1,
                                                    pi[sel,"start.pos"],pep.pos))

                         seen.sites <- rep(FALSE,max(modif.pos))
                         seen.sites[modif.pos] <- TRUE
                         which(seen.sites)
                 })
           names(res) <- isoform.acs
           res
         })
}

modif.site.count <- function(protein.group,protein.g=reporterProteins(protein.group),modif,take=max) {
  modif.sites <- modif.sites(protein.group,protein.g,modif)
  sapply(modif.sites,function(modifs.protein.g) {
         n.modifs.protein.acs <- sapply(modifs.protein.g,length)
         take(n.modifs.protein.acs)
  })
}

observedKnownSites <- function(protein.group,protein.g,ptm.info,modif,modification.name=NULL) {
  modif.sites <- modif.sites(protein.group,protein.g,modif)

  lapply(modif.sites,function(modifs.protein.g) {
         lapply(names(modifs.protein.g),function(isoform.ac) {
                observed.sites <- modifs.protein.g[[isoform.ac]]

                sel.known.sites <- ptm.info$isoform_ac==ifelse(grepl("-[0-9]$",isoform.ac),
                                                   isoform.ac,paste(isoform.ac,"-1",sep=""))
                if (!is.null(modification.name))
                  sel.known.sites <- sel.known.sites & grepl(modification.name,ptm.info$modification.name)

                c(n.observed.sites=length(observed.sites),
                  n.known.sites=sum(sel.known.sites),
                  n.known.sites.observed=sum(ptm.info[sel.known.sites,"first_position"] %in% observed.sites))
         })
  })
}


.proteinPtmInfo <- function(isoform.ac,protein.group,ptm.info,modif,modification.name=NULL,simplify=TRUE) {
  require(OrgMassSpecR)
  if (length(proteinInfo(protein.group)) == 0)
    stop("no protein info attached to protein.group: see ?getProteinInfoFromUniprot on how to get it.")

  protein.length <- as.numeric(proteinInfo(protein.group,protein.ac=isoform.ac,select="length") )
  if(all(is.na(protein.length))) 
    stop("no protein info for ",isoform.ac,"; need protein length an sequence")

 
  obs.peptides <- observable.peptides(proteinInfo(protein.group,protein.ac=isoform.ac,select="sequence"),nmc=2)
  possible.sites <- t(sapply(seq_len(protein.length),function(p) c(possible.nmc1=any(p>=obs.peptides$start & p<=obs.peptides$stop & obs.peptides$mc <=1),
                                                                   possible.nmc2=any(p>=obs.peptides$start & p<=obs.peptides$stop & obs.peptides$mc <=2))))
  my.ptm.info <- ptm.info[ptm.info$isoform_ac==ifelse(grepl("-[0-9]$",isoform.ac),
                                                   isoform.ac,paste(isoform.ac,"-1",sep="")),]
  if (!is.null(modification.name)) 
    my.ptm.info <- my.ptm.info[my.ptm.info$modification.name%in%modification.name,]
 
  # TO CHECK: first_position might be bigger than protein.length
  known.sites <- rep(FALSE,protein.length)
  if (nrow(my.ptm.info) > 0)
    known.sites[my.ptm.info$first_position] <- TRUE

  
  pi <- protein.group@peptideInfo
  sel.has.modif <- sapply(strsplit(pi[,"modif"],":"),function(x) any(x %in% modif))
  pi <- pi[pi[,"protein"]==isoform.ac & sel.has.modif,]

  pep.pos <- .convertModifToPos(pi[,"modif"],modif,simplify=FALSE,collapse=NULL)

  if (nrow(pi) > 0)
    modif.pos <- unlist(mapply(function(start.pos,pep.posi) start.pos + pep.posi -1,
                               pi[,"start.pos"],pep.pos))
  else
    modif.pos <- NULL

  seen.sites <- rep(FALSE,protein.length)
  seen.sites[modif.pos] <- TRUE

  if (simplify) {
  return(
         data.frame(observed.site.pos=paste(which(seen.sites),collapse=","),
           observed.sites=sum(seen.sites),
           known.sites=sum(known.sites),
           oberserved.known.sites=sum(known.sites&seen.sites),
           observable.known.sites.1mc=sum(known.sites&possible.sites[,"possible.nmc1"]),
           observable.known.sites.2mc=sum(known.sites&possible.sites[,"possible.nmc2"]),
           stringsAsFactors=FALSE)
         )
  } else {
  return(list(peptideInfo=pi,modif.pos=modif.pos,
              observable.peptides=obs.peptides,
              known.sites=my.ptm.info))

  }
}


getPeptideModifContext <- function(protein.group,modif,n.aa.up=5,n.aa.down=5) {
  
  peptide.info <- unique(peptideInfo(protein.group)[,c('modif','protein','peptide')])
  protein.sequences <- paste0(paste0(rep("_",n.aa.down),collapse=""),gsub("I","L",proteinInfo(protein.group)[,'sequence']),paste0(rep("_",n.aa.up),collapse="")) # enlarge sequence in case peptide starts in the beginning / end
  names(protein.sequences) <- proteinInfo(protein.group)[,'accession']
  
  pep.modif.context <- mapply(function(pepmodifs,protein.ac,pep) {
    my.seq <- protein.sequences[protein.ac]
    if (is.na(my.seq)) {
      warning("No sequence for ",protein.ac)
      return(NA)
    }
    
    peptide.startpos <- gregexpr(pep,my.seq)[[1]]
    if (peptide.startpos[1] == -1) {
      warning("Peptide [",pep,"] could not be matched to ",protein.ac)
      return(NA)    
    }

    pepmodifs[length(pepmodifs)] <- sub(" $","",pepmodifs[length(pepmodifs)])
    modification.positions <- which(pepmodifs%in%modif)-1
    if (length(modification.positions) == 0) 
      return(NA)
    
    paste(sapply(peptide.startpos,function(pep.pos) {
      modification.positions.in.protein <- pep.pos + modification.positions - 1
      res <- paste(sapply(modification.positions.in.protein,function(pos) 
                          substr(my.seq,pos-n.aa.down,pos+n.aa.up)),collapse=",")
    
      if (nchar(res) < n.aa.up+n.aa.down+1) 
        stop("extracted pattern does not have the length it should have")
      res
    }),collapse=",")
  },strsplit(paste0(peptide.info[,'modif']," "),":"),
    peptide.info[,'protein'], 
    peptide.info[,'peptide'])

  # get pep-modif context
  context.df <- unique(data.frame(peptide=peptide.info[,'peptide'],modif=peptide.info[,'modif'],context=pep.modif.context,stringsAsFactors=FALSE))
  print(head(context.df))
  ddply(context.df,
        c("peptide","modif"),
        function(x) data.frame(peptide=x[1,'peptide'],modif=x[1,'modif'],context=paste(x[,'context'],collapse=";"),stringsAsFactors=FALSE))

}
