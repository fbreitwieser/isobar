
#########################################
## Functions for XLS Report Generation 

testPerl <- function(perl.cmd) {
  test.cmd <- paste(perl.cmd, "-v")
  if (system(test.cmd) != 0)
    stop(perl.cmd," does not seem to work. It is required for spreadsheet reports in XLS and XLSX formats. ",
         "Either set 'write.xls.report=FALSE', or 'perl.cmd' to a different executable.") 
}

write.xls.report <- function(properties.env,report.env,file="isobar-analysis.xls") {

    testPerl(properties.env[["perl.cmd"]])

    write.t <- function(...,sep="\t",row.names=FALSE,quote=FALSE)
      write.table(...,sep=sep,row.names=row.names,quote=quote,na="")

    get.val <- function(name) { get(name,report.env,inherits=FALSE) }
    get.property <- function(name) { get(name,properties.env,inherits=FALSE) }
    protein.group <- proteinGroup(get.val('ibspectra'))
    indist.proteins <- indistinguishableProteins(protein.group)
    modificationSites <- NULL
    
    if (identical(properties.env[["report.level"]],"protein")) {
      ## TODO: add groups column to provide a link to the groups in report and quant table
      ## in principle, it works by defining the protein.group.ids attr, but does not here - 
      ## probably due to the environment it does not
      attr(proteinGroup(report.env[["ibspectra"]]),"protein.group.ids") <- .as.vect(unique(get.val('quant.tbl')[,c("ac","group")]))

      protein.id.df <- ibSpectra.as.concise.data.frame(get('ibspectra',report.env))
      ## make columns w/ multiple groups gray
      sel.1group  <- protein.id.df[["n.groups"]] == 1
      #sel.1ac  <- protein.id.df[["n.acs"]] == 1
      #sel.1variant  <- protein.id.df[["n.variants"]] == 1
      #protein.id.df[!sel.1ac & sel.1group,1] <- paste("#color silver#",protein.id.df[!sel.1ac & sel.1group,1],sep="")
      #protein.id.df[!sel.1group | !protein.id.df[["use.for.quant"]],1] <- paste("#color=gray#",protein.id.df[!sel.1group,1],sep="")

    } else {
      protein.id.df <- ibSpectra.as.concise.data.frame(get('ibspectra',report.env))
      protein.group <- proteinGroup(get('ibspectra',report.env))
      if (!is.null(properties.env[["phosphosite.dataset"]])) {
        sites <- do.call(rbind,lapply(properties.env[["phosphosite.dataset"]],
                                      read.delim,sep="\t",header=TRUE,skip=3,stringsAsFactors=FALSE))
        colnames(sites)[colnames(sites)=="ACC."]  <- "accession"
        colnames(sites) <- tolower(colnames(sites))
        modificationSites <- sites[sites$accession %in% proteins,]
      }
      proteins <- c(names(indistinguishableProteins(protein.group)),protein.ac(protein.group))
    }

    ## Analysis Properties:
    nn <- reporterTagNames(get.val('ibspectra'))
    ii <- rbind(c("@centeracross@Analysis Properties",rep("@centeracross@Analysis Properties",length(nn))),
                "",
                c("@centeracross@Isotope Impurity Correction Matrix",
                  rep("@centeracross@Isotope Impurity Correction Matrix",length(nn))),
                cbind(c("",nn),rbind(nn,isotopeImpurities(get.val('ibspectra')))))

    cl <- classLabels(get.val('ibspectra'))
    fill.up <- function(x,w="",n=length(nn)+1)
      if (length(x) <= n) c(x,rep(w,n-length(x)))
      else stop("Can't fill up - length(x) > n")
      
    if (!is.null(cl)) {
      ii <- rbind(ii,
                  "",
                  fill.up(c("@centeracross@Class Labels","@centeracross@Class Labels","")))
      
      for (i in seq_along(nn)) {
        ii <- rbind(ii,fill.up(c(nn[i],cl[i],names(cl)[i])))
      }
    }

    # TODO: add some info on number of quantified proteins
    # n.quant.sign <- ddply(quant.tbl,c("class1","class2"),function(x) sum(!is.na(x[["lratio"]]) & x[["is.significant"]]))
    # n.quant <- ddply(quant.tbl,c("class1","class2"),function(x) sum(!is.na(x[["lratio"]])))
    # n.noquant <- ddply(quant.tbl,c("class1","class2"),function(x) sum(is.na(x[["lratio"]])))

    ## write tables to tab seperated files files:
    protein.quant.f <- paste(get.property('cachedir'),"protein_quant.csv",sep="/")
    quant.notlocalized.f <- paste(get.property('cachedir'),"quant.notlocalized.csv",sep="/")
    protein.id.f <- paste(get.property('cachedir'),"protein_id.csv",sep="/")
    analysis.properties.f <- paste(get.property('cachedir'),"analysis_properties.csv",sep="/")
    log.f <- paste(get.property('cachedir'),"logged_operations.csv",sep="/")

    xls.quant.tbl <- get.val('xls.quant.tbl')
    if (identical(properties.env[["report.level"]],"protein")) {
      xls.quant.tbl <- xls.quant.tbl[order(xls.quant.tbl[,"group"]),]
      ## Create links
    } else {
      ## Create links
      protein.id.df[["peptide.sequence"]] <- protein.id.df[["peptide"]]
      if ("site.probs" %in% colnames(protein.id.df)) {
        protein.id.df[["site.probs"]] <- .convertPhosphoRSPepProb(protein.id.df[["peptide"]],protein.id.df[["site.probs"]])
      }
      protein.id.df[["peptide"]] <- .convertPeptideModif(protein.id.df[,"peptide"],protein.id.df[,"modif"])
      q.links <- sapply(protein.id.df[["peptide"]],function(p) {
                          res=which(xls.quant.tbl[["Sequence"]]==p)[1]
                            if (is.na(res)) ""
                            else paste0("@link=internal:Quantifications!A",res+1,"@q")
                  })
      protein.id.df <- cbind(q=q.links,protein.id.df)
      xls.quant.tbl <- cbind(i=paste0("@link=internal:Identifications!A",
                                      sapply(xls.quant.tbl[["Sequence"]],
                                             function(p) which(protein.id.df[["peptide"]]==p)[1]+1),
                                      "@",xls.quant.tbl[["Spectra"]]),xls.quant.tbl)
      #col_idx <- grep("Spectra", names(xls.quant.tbl))
      #colnames(xls.quant.tbl)[col_idx] <- "i"
      #xls.quant.tbl <- xls.quant.tbl[, c(col_idx, (1:ncol(df))[-col_idx])]

    }
    if (exists("xls.quant.tbl.notlocalized",envir=report.env)) {
      write.t(report.env[["xls.quant.tbl.notlocalized"]],file=quant.notlocalized.f)
    }

    if ('notes' %in% colnames(protein.id.df)) 
      protein.id.df[,'notes'] <- gsub("[\t\n]"," ",protein.id.df[,'notes'])
    write.t(xls.quant.tbl,file=protein.quant.f)
    write.t(protein.id.df,file=protein.id.f)  
    write.t(ii,file=analysis.properties.f,col.names=FALSE)
    write.t(get.val('ibspectra')@log,file=log.f,col.names=NA,row.names=TRUE)

    if (identical(properties.env[["report.level"]],"peptide") && !is.null(modificationSites)) {
      modifsites.f <- paste(get.property('cachedir'),"modification_sites.csv",sep="/")
      write.t(modificationSites,file=modifsites.f)
    }

    ## generate perl command line:
    tab2spreadsheet.cmd <- switch(properties.env[["spreadsheet.format"]],
                                  xlsx=system.file("pl","tab2xlsx.pl",package="isobar"),
                                  xls=system.file("pl","tab2xls.pl",package="isobar"),
                                  stop("spreadsheet.format property must be either 'xlsx' or 'xls'."))


    shq <- function(...) shQuote(paste0(...))
    
    .get.perlcl <- function(name,include.identifications=TRUE) {
      fname <- paste0(ifelse(properties.env[["use.name.for.report"]],
                             sprintf("%s.%s",properties.env[["name"]],name),"isobar-analysis"),
                       ".",properties.env[["spreadsheet.format"]])

      perl.cl <- paste(properties.env[["perl.cmd"]],
                       shq(tab2spreadsheet.cmd),shq(fname),
                       shq(":autofilter,freeze_col=4,name=Quantifications:",protein.quant.f))
      
      if (file.exists(quant.notlocalized.f))
        perl.cl <- paste(perl.cl,shq(":autofilter,freeze_col=4,name=Quantifications (not confidently localized sites):",
                            quant.notlocalized.f))
      
      if (identical(properties.env[["report.level"]],"peptide") && !is.null(modificationSites))
        perl.cl <- paste(perl.cl,shq(":autofilter,freeze_col=3,name=Modification Sites:",modifsites.f))

      if (include.identifications)
        perl.cl <- paste(perl.cl,shq(":autofilter,freeze_col=3,name=Identifications:",protein.id.f))

      perl.cl <- paste(perl.cl,shq(":name=Analysis Properties:",analysis.properties.f),
                                  shq(":name=Log:",log.f))
      perl.cl
    }
    
    ## generate Excel report (using Spreadsheet::WriteExcel)
    .call.cmd(.get.perlcl("quant"))
    .call.cmd(.get.perlcl("quantonly",FALSE))

}



.create.or.load.xls.quant.tbl <- function(env,properties.env,name="xls.quant.tbl",quant.tbl.name="quant.tbl") {
  .create.or.load(name,envir=properties.env,
                  msg.f=paste0("table for Excel export [",name,"]"),f=function() {
    message("XLS report format: ",properties.env[["xls.report.format"]])

    env[[quant.tbl.name]]$is.significant[is.na(env[[quant.tbl.name]]$is.significant)] <- FALSE
                    
    compare.to.quant <- .get.or.load('compare.to.quant',properties.env,class="data.frame",null.ok=TRUE)

    # Create a 'wide' XLS table (one row per protein / peptide)
    if (isTRUE(properties.env[["xls.report.format"]]=="wide")) {
      tbl.input  <- ratiosReshapeWide(env[[quant.tbl.name]],vs.class=properties.env[["vs.class"]],
                                      sep="###",cmbn=properties.env[["combn"]])

      if (!is.null(compare.to.quant))
        compare.to.quant <- lapply(compare.to.quant,ratiosReshapeWide,
                                   vs.class=properties.env[["vs.class"]],sep="###",cmbn=properties.env[["combn"]])
    } else {
      tbl.input <- env[[quant.tbl.name]]
    }

    message("   adding identification columns ",appendLF=FALSE)
    if (identical(properties.env[["report.level"]],"protein"))
      res.tbl <- .create.xls.protein.quant.tbl(tbl.input,proteinGroup(env[["ibspectra"]]))
    else
      res.tbl <- .create.xls.peptide.quant.tbl(tbl.input,env[["ibspectra"]],
                                           properties.env[["ptm"]],env[["ptm.info"]])
 
    message("   adding quantification columns")   
    tbl <- .add.quant.to.xls.tbl(env,properties.env,res.tbl[[1]],res.tbl[[2]],compare.to.quant)
    
    #order.c <- if(isTRUE(properties.env[["xls.report.format"]]=="long"),"Channels",NULL)
    if (identical(properties.env[["report.level"]],"protein"))
      tbl <- tbl[order(tbl[,"group"]),]
    else if (all(c("ID","Sequence") %in% colnames(tbl)))
      tbl <- tbl[order(tbl[,"ID"],tbl[,"Sequence"]),]

    tbl

  })
}

## XLS HELPER FUNCTIONS
.get.cols <- function(df,cc,cc.new=NULL,f=NULL,...) {
  data.cc <- df[,grep(cc,colnames(df)),drop=FALSE]
  if (!is.null(f)) data.cc <- f(data.cc,...)
  if (!is.null(cc.new)) colnames(data.cc) <- gsub(cc,cc.new,colnames(data.cc))
  data.cc
}

.add.quant.to.xls.tbl <- function(env,properties.env,tbl,input.tbl,compare.to.quant=NULL) {
  # Add quantification columns from compare.to.quant
  if (!is.null(compare.to.quant)){
    if (!"ac" %in% colnames(input.tbl)) {
      pnp <- subset(proteinGroup(env[["ibspectra"]])@peptideNProtein,
                    protein.g %in% reporterProteins(proteinGroup(env[["ibspectra"]])))
      pnp <- unlist(tapply(pnp[,"protein.g"],pnp[,"peptide"],paste,collapse=";",simplify=FALSE))
      input.tbl[["ac"]] <- pnp[input.tbl[["peptide"]]]
    }
    if (is.data.frame(compare.to.quant))
      compare.to.quant <- list(proteome=compare.to.quant)

    for (ii in seq_along(compare.to.quant)) 
      input.tbl = merge( input.tbl, compare.to.quant[[ii]], by=c("ac","r1","r2"), all.x=TRUE,
                         suffixes=c("",paste("###",names(compare.to.quant)[ii])) )
  }
  input.tbl <- input.tbl[order(input.tbl[["i"]]),]

  if (properties.env[["sum.intensities"]]) {
    message("summing intensities") 
    protein.intensities <- function(ib,proteins) {
      ri <- reporterIntensities(ib)
      t(sapply(proteins,
               function(p) {
                 sel <- spectrumSel(ib,protein=as.character(p),specificity=REPORTERSPECIFIC)
                 if(sum(sel) < 2)
                   ri[sel,]
                 else
                   colSums(ri[sel,],na.rm=TRUE)
               }
               ))
    }
    
    tbl <- cbind(tbl,protein.intensities(env[["ibspectra"]],tbl[["protein"]]))
  } else {
    if (isTRUE(properties.env[["xls.report.format"]]=="long")) {
     tbl <-cbind(tbl,"Classes"=paste(input.tbl[["class2"]],"/",input.tbl[["class1"]]))
     if (!all(input.tbl[['r2']]==input.tbl[['class2']])) {
       tbl <-cbind(tbl,"Channels"=paste(input.tbl[["r2"]],"/",input.tbl[["r1"]]))
      }
    }

    if ("zscore" %in% properties.env[["xls.report.columns"]]) {
      ## TODO: zscore is calculated across all classes - 
      ##       it is probably more appropriate to calculate it individual for each class
      #input.tbl[,'zscore'] <- .calc.zscore(input.tbl[,'lratio'])
    }
    round.digits <- 4;

.combine.n.append.xls.tbl <- function(cc1,cc2,cc.new,f) {
  A <- .get.cols(input.tbl,cc1)
  B <- .get.cols(input.tbl,cc2)
  data.cc <- sapply(1:ncol(A), function(i) f(A[,i],B[,i]) )
  data.cc <- round(data.cc,round.digits)
  colnames(data.cc) <- gsub(cc1,cc.new,colnames(A))
  data.cc
}
.append.xls.tbl <- function(...) { .get.cols(input.tbl,...) }

.round.n.appendend.xls.tbl <- function(...,digits=round.digits) {
  round(.get.cols(input.tbl,...),digits=digits)
}


  
    for (cc in properties.env[["xls.report.columns"]]) {
      message("    adding column ",cc)
      res <- switch(cc,
            log10.ratio =    .round.n.appendend.xls.tbl("lratio","@conditional_formatting=3_color_scale@log10.ratio"),
            log2.ratio =     .round.n.appendend.xls.tbl("lratio","@conditional_formatting=3_color_scale@log2.ratio",f=function(x) x/log10(2)),
            log10.variance = .append.xls.tbl("variance","log10.var"),
            log2.variance =  .append.xls.tbl("variance","log2.var",f=function(x) (sqrt(x)/log10(2)^2)),
            is.significant = .append.xls.tbl("is.significant"),
            n.na1 =          .append.xls.tbl("n.na1"),
            n.na2 =          .append.xls.tbl("n.na2"),
            p.value =        .append.xls.tbl("p.value"),
            p.value.adjusted = .append.xls.tbl("p.value.adjusted"),
            p.value.ratio =  .append.xls.tbl("p.value.rat"),
            p.value.ratio.adjusted =  .append.xls.tbl("p.value.rat.adjusted"),
            p.value.sample = .append.xls.tbl("p.value.sample"),
            z.score =        .round.n.appendend.xls.tbl("zscore"),
            ratio =          .round.n.appendend.xls.tbl("lratio","ratio",f=function(x) 10^x),
            CI95.lower =     .combine.n.append.xls.tbl("lratio","variance","CI95.lower",f=function(x,y) 10^qnorm(0.025,x,sqrt(y))),
            CI95.upper =     .combine.n.append.xls.tbl("lratio","variance","CI95.upper",f=function(x,y) 10^qnorm(0.975,x,sqrt(y))),
            ratio.minus.sd = .combine.n.append.xls.tbl("lratio","variance","ratio.minus.sd",f=function(x,y) 10^(x-sqrt(y))),
            ratio.plus.sd = .combine.n.append.xls.tbl("lratio","variance","ratio.plus.sd",f=function(x,y) 10^(x+sqrt(y))),
            warning("ignoring unknown column ",cc," in Excel report"))
      
      if (is(res,"data.frame") || is(res,"matrix"))
        tbl <- cbind(tbl,res)
      else
        warning("ignore res",cc)

    }
  }

  if (properties.env[["summarize"]]) {
    .append.xls.tbl("n.pos")
    .append.xls.tbl("n.neg")
  }

  if (length(properties.env[["preselected"]]) > 0) {
    ## tbl <- cbind(tbl,"is.preselected"=tbl[["is.preselected"]])
  }
  tbl[["i"]] <- NULL

  return(tbl)
}

.create.xls.peptide.quant.tbl <- function(input.tbl,ibspectra,ptm,ptm.info) {
  protein.group <- proteinGroup(ibspectra)
  ## PEPTIDE REPORT
  pnp  <- subset(as.data.frame(peptideNProtein(protein.group),stringsAsFactors=FALSE),
                 protein.g %in% reporterProteins(protein.group))

  my.ptm <- "PTM"
  if ("PHOS" %in% ptm) my.ptm="Phosphorylation"
  if ("METH" %in% ptm) my.ptm="Methylation"

  input.tbl[,'ac']  <- NULL
  pnp <- ddply(pnp,'peptide',function(x) c(peptide=x[1,'peptide'],ac=paste(x[,'protein.g'],collapse=";")))

  input.tbl <- merge(pnp,input.tbl,by="peptide")
  input.tbl[["i"]]  <- seq_len(nrow(input.tbl))
  input.tbl[["Spectra"]] <- apply(input.tbl[,grepl('^n.spectra',colnames(input.tbl)),drop=FALSE],1,max)
  pg.df <- .proteinGroupAsConciseDataFrame(protein.group,modif.pos=ptm,ptm.info=ptm.info)

  # show peptide modification context
  if (.proteinInfo.ok(protein.group)) {
    pep.modif.context <- getPeptideModifContext(protein.group,modif=ptm)
    pg.df <- merge(pg.df,pep.modif.context,by=c("peptide","modif"),all=TRUE)
  }
  
  tbl <- merge(pg.df,input.tbl[,intersect(c("peptide","pep.siteprobs","modif","i","Spectra"),colnames(input.tbl))],
               by=c("peptide","modif"),all.y=TRUE)
  tbl[["peptide"]] <- .convertPeptideModif(tbl[,"peptide"],tbl[,"modif"])
  colnames(tbl)[colnames(tbl)=="peptide"] <- "Sequence"
  
  colnames(tbl)[colnames(tbl)=="proteins"] <- "ACs"
  colnames(tbl)[colnames(tbl)=="modif.pos"] <- 
    paste0("@comment=Absolute modification position in protein. Modifications in ",
           "the same protein are separated by '&', in different proteins by ';'. ",
           "Stars denote positions which are annotated as phosphorylated in NextProt.@",
           my.ptm," Position")

  tbl[["start.pos"]] <- NULL
  tbl[["modif"]] <- NULL
  tbl[["pos"]] <- NULL
  tbl[["n.groups"]] <- NULL
  tbl <- tbl[order(tbl[["i"]]),]
  return(list(tbl,input.tbl))
} 

.proteinInfo.ok <- function(protein.group)
  is.data.frame(proteinInfo(protein.group)) && 
  length(proteinInfo(protein.group)) > 0

.create.xls.protein.quant.tbl <- function(input.tbl,protein.group,
                                          specificity=c(GROUPSPECIFIC,REPORTERSPECIFIC)) {

  message(".",appendLF=FALSE)
  input.tbl[["i"]]  <- seq_len(nrow(input.tbl))

  tbl.meta <- data.frame(i=input.tbl[["i"]], group=input.tbl[,"group"],protein.g=input.tbl[,'ac'],
                         stringsAsFactors=FALSE)
 
  # protein information
  protein.gs <- unique(input.tbl[,"ac"])
  tbl <- data.frame(protein.g=protein.gs,
                    AC=.protein.acc(protein.gs,protein.group),
                    stringsAsFactors=FALSE)

  message(".",appendLF=FALSE)
  indist.proteins <- .vector.as.data.frame(indistinguishableProteins(protein.group),colnames=c("protein","protein.g"))
  indist.proteins <- indist.proteins[indist.proteins[,'protein.g'] %in% protein.gs,]
  if (.proteinInfo.ok(protein.group)) {
    protein.info.tbl <- proteinInfo(protein.group,protein.g=protein.gs,select=c("name","protein_name","gene_name"))
    colnames(protein.info.tbl) <- c("ID","Description","Gene")
    tbl <- cbind(tbl,protein.info.tbl)
  }

  message(".",appendLF=FALSE)
  # spectra and peptide counts
  pnp <- as.data.frame(peptideNProtein(protein.group),stringsAsFactors=FALSE)
  ps <- peptideSpecificity(protein.group)
  protein.to.peptides <- merge(pnp[pnp[,'protein.g'] %in% protein.gs,],
                               ps[ps[,'specificity'] %in% specificity,],
                               by="peptide")
  protein.to.spectra <- merge(protein.to.peptides,
                              .vector.as.data.frame(protein.group@spectrumToPeptide,colnames=c("spectrum","peptide")),
                              by="peptide")

  tbl <- cbind(tbl, n=sapply(protein.gs,function(p) {sum(indist.proteins == p)}),
                    "@comment=Number of specific peptides@Peptide Count" = 
                       table(protein.to.peptides[,'protein.g'])[protein.gs],
                    "@comment=Number of specific spectra@Spectral Count" = 
                       table(protein.to.spectra[,'protein.g'])[protein.gs])

  message(".",appendLF=FALSE)
  # sequence coverage
  if (.proteinInfo.ok(protein.group) && "sequence" %in% colnames(proteinInfo(protein.group))) {
    peptide.info <- unique(peptideInfo(protein.group)[,c("protein","peptide","start.pos")])
    peptide.info[,'end.pos'] <- peptide.info[,'start.pos'] + nchar(peptide.info[,'peptide']) - 1
    peptide.info <- peptide.info[peptide.info[,'protein'] %in% indist.proteins[,'protein'],]
    protein.lengths <- nchar(unlist(setNames(proteinInfo(protein.group,protein.g=protein.gs,select=c("sequence"),simplify=FALSE),NULL)))

    if (!proteinInfoIsOnSpliceVariants(proteinInfo(protein.group))) {
      splice.df <- protein.group@isoformToGeneProduct
      splice.df <- splice.df[splice.df[,'proteinac.w.splicevariant'] %in% indist.proteins[,'protein'] &
                             splice.df[,'proteinac.wo.splicevariant'] %in% names(protein.lengths),]

      # only keep first AC
      splice.df.1 <- ddply(splice.df,'proteinac.wo.splicevariant',function(x) x[1,])
      peptide.info <- merge(peptide.info,splice.df.1,by.x='protein',by.y='proteinac.w.splicevariant')
      peptide.info[,'protein'] <- peptide.info[,'proteinac.wo.splicevariant']
    }
  
    seq.covs <- ddply(peptide.info,"protein",function(x) {
      protein.length <- protein.lengths[x[1,'protein']]
      if (is.na(protein.length)) return(c(seq.cov=NA))
      # remove invalid cases
      invalid.peptides <- x[,'end.pos'] > protein.length | is.na(x[,'start.pos']) | 
          is.nan(x[,'end.pos']) | is.na(x[, 'end.pos']) | is.nan(x[, 'end.pos'])
      x <- x[!invalid.peptides, ]
      
      if (nrow(x) == 0) return(c(seq.cov=NA))
      seqq <- rep(FALSE,protein.length);
      for (i in seq_along(nrow(x))) {
        seqq[x[i,'start.pos']:x[i,'end.pos']] <-  TRUE
      }
      return(c(seq.cov=sum(seqq)/length(seqq)))
    },.parallel=isTRUE(options('isobar.parallel')))
  
    if (!proteinInfoIsOnSpliceVariants(proteinInfo(protein.group))) {
      seq.covs <- merge(seq.covs,splice.df,by.x='protein',by.y='proteinac.wo.splicevariant')
      seq.covs[,'protein'] <- seq.covs[,'proteinac.w.splicevariant']
    }

    proteing.seq.covs <- merge(indist.proteins,seq.covs,by='protein')
    seq.covs <- tapply(proteing.seq.covs[,'seq.cov'],factor(proteing.seq.covs[,'protein.g']),mean)
  
    tbl <- cbind(tbl, "@comment=Coverage of the protein sequence with observed peptides@Sequence Coverage" = 
                        round(seq.covs[protein.gs],digits=4))
  }
  tbl <- merge(tbl.meta,tbl,by='protein.g')
  tbl[["protein.g"]] <- NULL
  tbl <- tbl[order(tbl[,'i']),]
  if (!all(tbl[,'i']==input.tbl[,'i'])) stop ("problem in protein xls report ordering")

  return(list(tbl,input.tbl))
}


