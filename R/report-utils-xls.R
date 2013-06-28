
#########################################
## Functions for XLS Report Generation 

write.xls.report <- function(report.type,properties.env,report.env,file="isobar-analysis.xls") {
    write.t <- function(...,sep="\t",row.names=FALSE,quote=FALSE)
      write.table(...,sep=sep,row.names=row.names,quote=quote,na="")

    get.val <- function(name) { get(name,report.env,inherits=FALSE) }
    get.property <- function(name) { get(name,properties.env,inherits=FALSE) }
    protein.group <- proteinGroup(get.val('ibspectra'))
    indist.proteins <- indistinguishableProteins(protein.group)
    modificationSites <- NULL
    
    if (identical(report.type,"protein")) {
      ## TODO: add groups column to provide a link to the groups in report and quant table
      ## in principle, it works by defining the protein.group.ids attr, but does not here - 
      ## probably due to the environment it does not
      attr(proteinGroup(report.env$ibspectra),"protein.group.ids") <- .as.vect(unique(get.val('quant.tbl')[,c("ac","group")]))

      protein.id.df <- ibSpectra.as.concise.data.frame(get('ibspectra',report.env))
      ## make columns w/ multiple groups gray
      sel.1group  <- protein.id.df$n.groups == 1
      #sel.1ac  <- protein.id.df$n.acs == 1
      #sel.1variant  <- protein.id.df$n.variants == 1
      #protein.id.df[!sel.1ac & sel.1group,1] <- paste("#color silver#",protein.id.df[!sel.1ac & sel.1group,1],sep="")
      #protein.id.df[!sel.1group | !protein.id.df$use.for.quant,1] <- paste("#color=gray#",protein.id.df[!sel.1group,1],sep="")

    } else {
      protein.id.df <- ibSpectra.as.concise.data.frame(get('ibspectra',report.env))
      protein.group <- proteinGroup(get('ibspectra',report.env))
      if (!is.null(properties.env$phosphosite.dataset)) {
        sites <- do.call(rbind,lapply(properties.env$phosphosite.dataset,
                                      read.delim,sep="\t",header=TRUE,skip=3,stringsAsFactors=FALSE))
        colnames(sites)[colnames(sites)=="ACC."]  <- "accession"
        colnames(sites) <- tolower(colnames(sites))
        modificationSites <- subset(sites,accession %in% proteins)
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
    # n.quant.sign <- ddply(quant.tbl,c("class1","class2"),function(x) sum(!is.na(x$lratio) & x$is.significant))
    # n.quant <- ddply(quant.tbl,c("class1","class2"),function(x) sum(!is.na(x$lratio)))
    # n.noquant <- ddply(quant.tbl,c("class1","class2"),function(x) sum(is.na(x$lratio)))

    ## write tables to tab seperated files files:
    protein.quant.f <- paste(get.property('cachedir'),"protein_quant.csv",sep="/")
    protein.id.f <- paste(get.property('cachedir'),"protein_id.csv",sep="/")
    analysis.properties.f <- paste(get.property('cachedir'),"analysis_properties.csv",sep="/")
    log.f <- paste(get.property('cachedir'),"logged_operations.csv",sep="/")

    xls.quant.tbl <- get.val('xls.quant.tbl')
    if (report.type == "protein") {
      xls.quant.tbl <- xls.quant.tbl[order(xls.quant.tbl[,"group"]),]
      ## Create links
    } else {
      ## Create links
      protein.id.df$peptide.sequence <- protein.id.df$peptide
      if ("site.probs" %in% colnames(protein.id.df)) {
        protein.id.df$site.probs <- .convertPhosphoRSPepProb(protein.id.df$peptide,protein.id.df$site.probs)
      }
      protein.id.df$peptide <- .convertPeptideModif(protein.id.df[,"peptide"],protein.id.df[,"modif"])
      q.links <- sapply(protein.id.df$peptide,function(p) {
                          res=which(xls.quant.tbl$Sequence==p)[1]
                            if (is.na(res)) ""
                            else paste0("@link=internal:Quantifications!A",res+1,"@q")
                  })
      protein.id.df <- cbind(q=q.links,protein.id.df)
      xls.quant.tbl <- cbind(i=paste0("@link=internal:Identifications!A",
                                      sapply(xls.quant.tbl$Sequence,
                                             function(p) which(protein.id.df$peptide==p)[1]+1),
                                      "@",xls.quant.tbl$Spectra),xls.quant.tbl)
      #col_idx <- grep("Spectra", names(xls.quant.tbl))
      #colnames(xls.quant.tbl)[col_idx] <- "i"
      #xls.quant.tbl <- xls.quant.tbl[, c(col_idx, (1:ncol(df))[-col_idx])]

    }
    if ('notes' %in% colnames(protein.id.df)) 
      protein.id.df[,'notes'] <- gsub("[\t\n]"," ",protein.id.df[,'notes'])
    write.t(xls.quant.tbl,file=protein.quant.f)
    write.t(protein.id.df,file=protein.id.f)  
    write.t(ii,file=analysis.properties.f,col.names=FALSE)
    write.t(get.val('ibspectra')@log,file=log.f,col.names=NA,row.names=TRUE)

    if (identical(report.type,"peptide") && !is.null(modificationSites)) {
      modifsites.f <- paste(get.property('cachedir'),"modification_sites.csv",sep="/")
      write.t(modificationSites,file=modifsites.f)
    }

    ## generate perl command line:
    tab2spreadsheet.cmd <- switch(properties.env$spreadsheet.format,
                                  xlsx=system.file("pl","tab2xlsx.pl",package="isobar",mustWork=TRUE),
                                  xls=system.file("pl","tab2xls.pl",package="isobar",mustWork=TRUE),
                                  stop("spreadsheet.format property must be either 'xlsx' or 'xls'."))

    perl.cl <- paste(tab2spreadsheet.cmd," ",
                     ifelse(properties.env$use.name.for.report,sprintf("%s.quant",properties.env$name),"isobar-analysis"),
                     ".",properties.env$spreadsheet.format,
                     " ':autofilter,freeze_col=6,name=Quantifications:",protein.quant.f,"'",
                     ifelse(identical(report.type,"peptide") && !is.null(modificationSites),
                            paste(" ':autofilter,freeze_col=3,name=Modification Sites:",modifsites.f,"'",sep=""),""),
                     " ':autofilter,freeze_col=3,name=Identifications:",protein.id.f,"'",
                     " ':name=Analysis Properties:",analysis.properties.f,"'",
                     " ':name=Log:",log.f,"'",sep="")
    
    ## generate Excel report (using Spreadsheet::WriteExcel)
    message(perl.cl)
    system(perl.cl)

    perl.cl <- paste(tab2spreadsheet.cmd," ",
                     ifelse(properties.env$use.name.for.report,sprintf("%s.quantonly",properties.env$name),"isobar-analysis-quantonly"),
                     ".",properties.env$spreadsheet.format,
                     " ':autofilter,freeze_col=6,name=Quantifications:",protein.quant.f,"'",
                     ifelse(identical(report.type,"peptide") && !is.null(modificationSites),
                            paste(" ':autofilter,freeze_col=3,name=Modification Sites:",modifsites.f,"'",sep=""),""),
                     " ':name=Analysis Properties:",analysis.properties.f,"'",
                     " ':name=Log:",log.f,"'",sep="")
    
    ## generate Excel report (using Spreadsheet::WriteExcel)
    message(perl.cl)
    system(perl.cl)

}



.create.or.load.xls.quant.tbl <- function(report.type,env,properties.env) {
  .create.or.load("xls.quant.tbl",envir=properties.env,
                  msg.f="protein table for Excel export",f=function() {
    message("XLS report format: ",properties.env$xls.report.format)
                    
    compare.to.quant <- .get.or.load('compare.to.quant',properties.env,class="data.frame",null.ok=TRUE)

    # Create a 'wide' XLS table (one row per protein / peptide)
    if (isTRUE(properties.env$xls.report.format=="wide")) {
      tbl.input  <- ratiosReshapeWide(env$quant.tbl,vs.class=properties.env$vs.class,
                                      sep="###",cmbn=properties.env$combn)

      if (!is.null(compare.to.quant))
        compare.to.quant <- lapply(compare.to.quant,ratiosReshapeWide,
                                   vs.class=properties.env$vs.class,sep="###",cmbn=properties.env$combn)
    } else {
      tbl.input <- env$quant.tbl
    }

    if (identical(report.type,"protein"))
      res.tbl <- .create.xls.protein.quant.tbl(tbl.input,proteinGroup(env$ibspectra))
    else
      res.tbl <- .create.xls.peptide.quant.tbl(tbl.input,env$ibspectra,
                                           properties.env$ptm,env$ptm.info)
 
    message(" adding quantification columns")   
    tbl <- .add.quant.to.xls.tbl(env,properties.env,res.tbl[[1]],res.tbl[[2]],compare.to.quant)
    
    #order.c <- if(isTRUE(properties.env$xls.report.format=="long"),"Channels",NULL)
    if (identical(report.type,"protein"))
      tbl <- tbl[order(tbl[,"group"]),]
    else 
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
      pnp <- subset(proteinGroup(env$ibspectra)@peptideNProtein,
                    protein.g %in% reporterProteins(proteinGroup(env$ibspectra)))
      pnp <- unlist(tapply(pnp[,"protein.g"],pnp[,"peptide"],paste,collapse=";",simplify=FALSE))
      input.tbl$ac <- pnp[input.tbl$peptide]
    }
    if (is.data.frame(compare.to.quant))
      compare.to.quant <- list(proteome=compare.to.quant)

    for (ii in seq_along(compare.to.quant)) 
      input.tbl = merge( input.tbl, compare.to.quant[[ii]], by=c("ac","r1","r2"), all.x=TRUE,
                         suffixes=c("",paste("###",names(compare.to.quant)[ii])) )
  }
  input.tbl <- input.tbl[order(input.tbl$i),]

  if (properties.env$sum.intensities) {
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
    
    tbl <- cbind(tbl,protein.intensities(ibspectra,protein.tbl$protein))
  } else {
    if (isTRUE(properties.env$xls.report.format=="long")) {
     tbl <-cbind(tbl,"Channels"=paste(input.tbl$r2,"/",input.tbl$r1))
    }

    if ("zscore" %in% properties.env$xls.report.columns) {
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


  
    for (cc in properties.env$xls.report.columns) {
      message("    adding column ",cc)
      res <- switch(cc,
            log10.ratio =    .round.n.appendend.xls.tbl("lratio","@conditional_formatting=3_color_scale@log10.ratio"),
            log2.ratio =     .round.n.appendend.xls.tbl("lratio","@conditional_formatting=3_color_scale@log2.ratio",f=function(x) x/log10(2)),
            log10.variance = .append.xls.tbl("variance","log10.var"),
            log2.variance =  .append.xls.tbl("variance","log2.var",f=function(x) (sqrt(x)/log10(2)^2)),
            is.significant = .append.xls.tbl("is.significant"),
            n.na1 =          .append.xls.tbl("n.na1"),
            n.na2 =          .append.xls.tbl("n.na2"),
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

  if (properties.env$summarize) {
    .append.xls.tbl("n.pos")
    .append.xls.tbl("n.neg")
  }

  if (length(properties.env$preselected) > 0) {
    ## tbl <- cbind(tbl,"is.preselected"=tbl$is.preselected)
  }
  tbl$i <- NULL

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
  #t <- table(pnp$peptide)
  #pnp <- pnp[pnp$peptide %in% names(t)[t==1],]
  #colnames(pnp)  <- c("peptide","ac")
  pnp <- ddply(pnp,'peptide',function(x) c(peptide=x[1,1],ac=paste(x[,2],collapse=";")))

  input.tbl <- merge(pnp,input.tbl,by="peptide")
  input.tbl$i  <- seq_len(nrow(input.tbl))
  input.tbl$Spectra <- apply(input.tbl, 1, 
                             function(x) nrow(subset(fData(ibspectra),peptide==x['peptide'] & modif==x['modif'])))
  pg.df <- .proteinGroupAsConciseDataFrame(protein.group,modif.pos=ptm,ptm.info=ptm.info)
  #if (isTRUE(properties.env$show.motifs)) {
    pep.modif.context <- getPeptideModifContext(protein.group,modif=ptm)
    pg.df <- merge(pg.df,pep.modif.context,by=c("peptide","modif"),all=TRUE)
  #}
  
  tbl <- merge(pg.df,input.tbl[,c("peptide","modif","i","Spectra")],by=c("peptide","modif"),all.y=TRUE)
  tbl$peptide <- .convertPeptideModif(tbl[,"peptide"],tbl[,"modif"])
  colnames(tbl)[colnames(tbl)=="peptide"] <- "Sequence"
  
  colnames(tbl)[colnames(tbl)=="proteins"] <- "ACs"
  colnames(tbl)[colnames(tbl)=="modif.pos"] <- 
    paste0("@comment=Absolute modification position in protein. Modifications in ",
           "the same protein are separated by '&', in different proteins by ';'. ",
           "Stars denote positions which are annotated as phosphorylated in NextProt.@",
           my.ptm," Position")

  tbl$start.pos <- NULL
  tbl$modif <- NULL
  tbl$pos <- NULL
  tbl$n.groups <- NULL
  tbl <- tbl[order(tbl$i),]
  return(list(tbl,input.tbl))
} 

.create.xls.protein.quant.tbl <- function(input.tbl,protein.group,
                                          specificity=c(GROUPSPECIFIC,REPORTERSPECIFIC)) {

  indist.proteins <- indistinguishableProteins(protein.group)
  input.tbl$i  <- seq_len(nrow(input.tbl))

  proteinInfo.ok <- is.data.frame(proteinInfo(protein.group)) && 
                      length(proteinInfo(protein.group)) > 0

  .getProteinInfo <- function(name) 
    proteinInfo(protein.group,protein.g=input.tbl[,"ac"],select=name,do.warn=FALSE)

  ## creating xls protein table
  tbl <- data.frame(i=input.tbl$i, group=input.tbl[,"group"],
                    AC=.protein.acc(input.tbl[,"ac"],protein.group))

  if (proteinInfo.ok)
    tbl <- cbind(tbl, ID=.getProteinInfo("name"),
                      Description=.getProteinInfo("protein_name"),
                      Gene=.getProteinInfo("gene_name"))

  tbl <- cbind(tbl, n=sapply(input.tbl[,"ac"],function(p) {length(names(indist.proteins)[indist.proteins == p])}),
                    "@comment=Number of specific peptides@Peptide Count" = 
                      peptide.count(protein.group,input.tbl[,'ac'],specificity=specificity,do.warn=FALSE),
                    "@comment=Number of specific spectra@Spectral Count" = 
                      spectra.count(protein.group,input.tbl[,'ac'],specificity=specificity,do.warn=FALSE))

  if (proteinInfo.ok && "length" %in% colnames(proteinInfo(protein.group))) {
    tbl <- cbind(tbl, "@comment=Coverage of the protein sequence with observed peptides@Sequence Coverage" = 
                        round(sequence.coverage(protein.group,input.tbl[,'ac'],do.warn=FALSE),digits=4))
  }
  return(list(tbl,input.tbl))
}


