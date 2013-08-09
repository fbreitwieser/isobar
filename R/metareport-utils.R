create.meta.reports <- function(properties.file="meta-properties.R",
                                global.properties.file=system.file("report","meta-properties.R",package="isobar"),
                                args=NULL,include.not.localized=FALSE,...) {
  require(ggplot2)
  if (!exists("properties.env")) {
    properties.env <- load.properties(properties.file,global.properties.file,
                                      args=args,...)
  }

  protein.group <- .get.or.load("protein.group",properties.env,"protein group object","ProteinGroup")

  pg.df <- .create.or.load("pg.df",properties.env,function() {
    if (properties.env[["report.level"]]=="peptide") {
      ptm.info  <- .create.or.load.ptm.info(properties.env,protein.group)
      pg.df <- isobar:::.proteinGroupAsConciseDataFrame(protein.group,
                                                        modif.pos=properties.env$ptm,
                                                        ptm.info=properties.env$ptm.info,
                                                        report.only.modified=TRUE)
    } else {
      acs <- reporterProteins(protein.group)
      pg.df <- data.frame(ac=acs,
                          AC=.protein.acc(acs,protein.group),
                          ID=proteinInfo(protein.group,acs,do.warn=FALSE),
                          Description=proteinInfo(protein.group,acs,select="protein_name",do.warn=FALSE),
                          Gene=proteinInfo(protein.group,acs,select="gene_name",do.warn=FALSE),
                          stringsAsFactors=FALSE)

    }
  })

  if (length(proteinInfo(protein.group)) ==0 && !is.null(properties.env$protein.info.f)) 
    proteinInfo(protein.group) <- 
      .create.or.load("protein.info",envir=properties.env,
                      f=properties.env$protein.info.f,
                      x=protein.group,
                      error=warning,default.value=proteinInfo(protein.group))
  
  if (is.null(proteinInfo(protein.group)) || length(proteinInfo(protein.group))==0)
    stop("No protein information available.")

  #if (properties.env$calculate.dnsaf)
  #  dnsaf <- .create.or.load("dnsaf",envir=properties.env,f=calculate.dNSAF,protein.group=protein.group)

  message("Merging samples ",paste(properties.env$samples,collapse=", ")," ...")
  ac.vars <- switch(properties.env[["report.level"]],
                    protein = "ac",
                    peptide = c("peptide","modif","pep.siteprobs"),
                    stop("properties.env$report.level ",properties.env[["report.level"]]," unknown"))

  merge.cols <- c("class1","class2")

  merged.table <- .create.or.load.merged.table(properties.env=properties.env)
  merged.table.notlocalized <- .create.or.load.merged.table("merged.table.notlocalized","quant.tbl.notlocalized.rda",properties.env=properties.env)
  ac.vars <- intersect(ac.vars,colnames(merged.table))
  
  tbl.pg <- .create.or.load("tbl.pg",envir=properties.env,f=function()
                            .get.tbl.pg(merged.table,pg.df,merge.cols,ac.vars))

  if (exists("merged.table.notlocalized")) {
    tbl.pg.nl <- .create.or.load("tbl.pg.nl",envir=properties.env,f=function()
                                  .get.tbl.pg(merged.table.notlocalized,
                                             pg.df,merge.cols,ac.vars))
  }
#  .create.figures(tbl.pg)

  #tbl.wide <- .create.or.load("tbl.wide",properties.env,function() reshape(merged.table,idvar=intersect(ac.vars,colnames(merged.table)),timevar="comp",direction="wide",drop=merge.cols))

  #.plot.data(properties.env,tbl.wide,ac.vars)

  .write.xls.report(properties.env,tbl.pg,tbl.pg.nl,include.not.localized)
}

.plot.data <- function(properties.env,tbl.wide,ac.vars) {
  message("plotting data")
  ratio.matrix <- as.matrix(tbl.wide[,grep("lratio",colnames(tbl.wide))])
  variance.matrix <- as.matrix(tbl.wide[,grep("var",colnames(tbl.wide))])
  rownames(ratio.matrix)  <- do.call(paste,tbl.wide[,ac.vars,drop=FALSE])
  rownames(variance.matrix)  <- do.call(paste,tbl.wide[,ac.vars,drop=FALSE])
  sel <- !apply(is.na(ratio.matrix),1,any)
  ratio.matrix <- ratio.matrix[sel,]
  variance.matrix <- variance.matrix[sel,]

  m.median <- apply(ratio.matrix,2,median)
  normalized.ratio.matrix <- ratio.matrix-matrix(m.median,nrow=nrow(ratio.matrix),ncol=ncol(ratio.matrix),byrow=T)

  .plot.heatmaps(ratio.matrix,normalized.ratio.matrix,variance.matrix,properties.env$name)
  .plot.pairs(ratio.matrix,variance.matrix,properties.env$name)
}


.write.xls.report <- function(properties.env,tbl.pg,tbl.pg.nl,include.not.localized=FALSE) {
  message("writing Excel report")
  xlsxname <- paste0(properties.env$name,"-combined_report.xlsx")
  tab2spreadsheet.cmd <- switch(properties.env$spreadsheet.format,
                                xlsx=system.file("pl","tab2xlsx.pl",package="isobar"),
                                xls=system.file("pl","tab2xls.pl",package="isobar"),
                                stop("spreadsheet.format property must be either 'xlsx' or 'xls'."))

  csvname <- paste0(properties.env$name,"-combined_report.csv")
  write.table(tbl.pg[order(tbl.pg[,'ID'],tbl.pg[,'start.pos']),],file=csvname,
              row.names=FALSE,sep="\t",na="")
  perl.cl <- paste(tab2spreadsheet.cmd,shQuote(xlsxname),shQuote(csvname))
  
  if (include.not.localized) {
    csvname.nl <- paste0(properties.env$name,".not.localized-combined_report.csv")
    write.table(tbl.pg.nl[order(tbl.pg.nl[,'ID'],tbl.pg[,'start.pos']),],file=csvname.nl,row.names=FALSE,sep="\t",na="")
    perl.cl <- paste(perl.cl,shQuote(csvname.nl))
  }
  ## generate Excel report (using Spreadsheet::WriteExcel)
  message(perl.cl)
  system(perl.cl)
}


.get.tbl.pg <- function(merged.table,pg.df,merge.cols,ac.vars) {
  tbl.wide <- reshape(merged.table,idvar=intersect(ac.vars,colnames(merged.table)),timevar="comp",direction="wide",drop=merge.cols)
  tbl.pg <- merge(pg.df,tbl.wide)#,by=intersect(ac.vars,colnames(tbl.wide)))

  sel <- apply(!is.na(tbl.pg[,grep("^lratio",colnames(tbl.pg))]),1,any)
  tbl.pg <- tbl.pg[sel,]

  tbl.pg <- tbl.pg[order(rowSums(tbl.pg[,grep("is.significant",colnames(tbl.pg))]),
                         rowSums(tbl.pg[,grep("lratio",colnames(tbl.pg))]),
                         decreasing=TRUE),]

  if ("peptide" %in% colnames(tbl.pg)) {
    tbl.pg$peptide <- isobar:::.convertPeptideModif(tbl.pg$peptide,tbl.pg$modif)
    tbl.pg$modif <- NULL
  }
  tbl.pg$ac <- NULL

  print(t(sapply(.grep_columns(tbl.wide,"lratio.P",value=TRUE,logical=FALSE),
                 function(i) isobar:::.sum.bool(!is.na(tbl.wide[,i])))))


  tbl.pg
}

.create.or.load.merged.table <- function(name="merged.table",fname="quant.tbl.rda",properties.env) {
  ac.vars <- switch(properties.env[["report.level"]],
                    protein = "ac",
                    peptide = c("peptide","modif","pep.siteprobs"),
                    stop("properties.env$report.level ",properties.env[["report.level"]]," unknown"))

  merge.cols <- c("class1","class2")

  res <- .create.or.load(name,properties.env,function() {
    merged.table <- .get.merged.table(properties.env$samples,
                                     cols=c(ac.vars,merge.cols,"lratio","variance","p.value.rat"),
                                     merge.by=c(ac.vars,merge.cols),format=properties.env[["xls.report.format"]],
                                     fname=fname)
  
    zscore.any2.5 <- apply(abs(merged.table[,grepl("zscore.indiv",colnames(merged.table))]) >= properties.env[["zscore.threshold"]],1,any,na.rm=TRUE)
    zscore.all1 <- apply(abs(merged.table[,grepl("zscore.indiv",colnames(merged.table))]) >= 1,1,all,na.rm=TRUE)

    merged.table$p.value.rat <- .combine.fisher.tblwide(merged.table)
    merged.table$is.significant.notadj <- (merged.table$p.value.rat < 0.05) & zscore.any2.5 & zscore.all1

    #merged.table$p.value.rat.fdr.adj <- p.adjust(merged.table$p.value.rat)
    merged.table <- ddply(merged.table,c("class1","class2"),function(x) 
                                 cbind(x,p.value.rat.fdr.adj=p.adjust(x$p.value.rat,"fdr")))
    merged.table$is.significant <- (merged.table$p.value.rat.fdr.adj < 0.05) & zscore.any2.5 & zscore.all1  
  
    merged.table$comp <- paste(merged.table[,merge.cols[2]],sep="/",merged.table[,merge.cols[1]])
    merged.table
  })
  #merged.table$p.value.rat <- .combine.fisher.tblwide(merged.table)
  return(res)
}

.calc.zscore <- function(qq) {
  s.median <- median(qq,na.rm=TRUE)
  s.mad <- mad(qq,na.rm=TRUE)
  (qq-s.median)/s.mad
}

.get.merged.table <- function(samples,quant.tbls=paste0(samples,"/cache/",fname),
                              cols=c("ac","r1","r2","lratio","variance"),
                              merge.by=c("ac","r1","r2"),format="wide",fname="quant.tbl.rda") {

  quant.tables <- lapply(seq_along(samples),
                         function(idx) {
                           message("loading ",quant.tbls[idx])
                           env <- new.env()
                           load(quant.tbls[idx],envir=env)
                           quant.tbl <- get(ls(envir=env),envir=env)
                           q <- quant.tbl[,intersect(cols,colnames(quant.tbl))]
                           q <- ddply(q,c("class1","class2"),function(x) 
                                 cbind(x,zscore.indiv=round(.calc.zscore(x[,"lratio"]),4)))
                           q[,"lratio"] <- round(q[,"lratio"],4)
                           q[,"variance"] <- round(q[,"variance"],4)
                           if (format=="wide") {
                             sel <- !colnames(q) %in% merge.by
                             colnames(q)[sel] <- paste(colnames(q)[sel],gsub("[./]","",samples[idx]),sep=".")
                           } else {
                             q$sample <- samples[idx]
                           }
                           message(paste(colnames(q),collapse=":"))
                           q
                         })

  if (format == "wide") {
    merged.table <- quant.tables[[1]]
    for (idx in  2:length(samples)) {
      merge.colnames <- intersect(merge.by,colnames(quant.tables[[idx]]))

      merged.table <- merge(merged.table,
                            quant.tables[[idx]],by=merge.colnames,all=TRUE)
    }
    merged.table$lratio <- rowMeans(merged.table[,.grep_columns(merged.table,"lratio")],na.rm=TRUE)
    merged.table$variance <- rowMeans(merged.table[,.grep_columns(merged.table,"variance")],na.rm=TRUE)
  } else {
    merged.table <- do.call(rbind,quant.tables)
  }
  return(merged.table)
}

.get.names <- function(p,protein.group) 
  apply(my.protein.info(protein.group,p)[,c("name","gene_name","protein_name")],2,
        function(s) paste(unique(sort(gsub("'","",s))),collapse=", "))

.calc.col <- function(tbl.wide,tag,add.name=TRUE) {
  idx.ratio <- grep("lratio",colnames(tbl.wide),fixed=TRUE)
  idx.var <- grep("var",colnames(tbl.wide),fixed=TRUE)
  x <- t(apply(tbl.wide,1,function(r) {
    pos <- grep(tag,colnames(tbl.wide),fixed=TRUE)
    idx.r <- intersect(pos,idx.ratio)
    idx.v <- intersect(pos,idx.var)
  
    r<-as.numeric(r)
    c(weighted.mean=round(10^weightedMean(r[idx.r],1/r[idx.v]),2),
      var=round(weightedVariance(r[idx.r],1/r[idx.v]),4),
      n.up=sum(r[idx.r]>0,na.rm=T),
      n.down=sum(r[idx.r]<0,na.rm=T))
  }))
  if (add.name)
    colnames(x) <- paste(colnames(x),tag,sep=".")
  return(x)
}

.write.t <- function(...,sep="\t",row.names=FALSE,quote=FALSE)
  write.table(...,sep=sep,row.names=row.names,quote=quote)

.write.summarized.table <- function(tbl.wide,allnames,cols) {
  summarized.table <- cbind(ac=rownames(tbl.wide),allnames)
  for (col in cols)
    summarized.table <- cbind(summarized.table,.calc.col(tbl.wide,col,add.name=length(cols)>1))
  summarized.table <- cbind(summarized.table,tbl.wide)
  .write.t(summarized.table,file=paste(name,"_summarized_table.csv",sep=""))
  return(summarized.table)
}

.plot.heatmaps.gd <- function(ratio.matrix,variance.matrix,name) {
  require(gplots)
  breaks <- seq(from=-max(abs(ratio.matrix)),to=max(abs(ratio.matrix)),length.out=51)
  pdf(sprintf("heatmap_%s.pdf",name),width=15,height=30,title=name)
  min.max <- quantile(ratio.matrix,probs=c(0.0005,0.9995))
  sel <- apply(ratio.matrix<min.max[1] | ratio.matrix > min.max[2],1,any)
  if (sum(sel) > 2) {
    ratio.matrix2 <- ratio.matrix[sel,,drop=FALSE]
    heatmap.2(ratio.matrix2,Colv=NA,col=greenred(50),dendrogram="row",
              margins=c(5,25),main=paste(name,"- above fold change of",round(10^min.max[2],1),
                                "or below",round(10^min.max[1],1)),
              key=FALSE, keysize=1.0, symkey=FALSE, density.info='none',
              trace='none', colsep=1:10,
              #labRow=allnames[sel,"gene_name"],
              sepcolor='white', sepwidth=0.025,
              scale="none",cexRow=2,cexCol=2,
              labCol = colnames(ratio.matrix),                 
              hclustfun=function(c){hclust(c, method='mcquitty')},
              lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),breaks=breaks)
  }
  dev.off()

  library(boot)

  cols <- colnames(ratio.matrix)
  correlation.matrix <- sapply(cols ,function(i) sapply(cols,function(j) boot::corr(ratio.matrix[,c(i,j)],
                                                            w=apply(1/variance.matrix,1,median))))

  pdf(sprintf("heatmap_correlation_%s.pdf",name),width=15,height=30,title=name)
  heatmap.2(correlation.matrix,scale="none",Rowv=NA,Colv=NA,col=greenred(50),
            density.info='none',trace="none",cellnote=round(correlation.matrix,1),main="not re")
  dev.off()

}


.plot.heatmaps <- function(ratio.matrix,normalized.ratio.matrix,variance.matrix,name) {
  require(gplots)
  breaks <- seq(from=-max(abs(ratio.matrix)),to=max(abs(ratio.matrix)),length.out=51)
  pdf(sprintf("heatmap_%s.pdf",name),width=15,height=30,title=name)
  heatmap.2(ratio.matrix,Colv=NA,col=greenred(50),margins=c(5,25),main=paste(name),
#            labRow=allnames[,"gene_name"],
            lmat=rbind( c(0, 3), c(2,1), c(4,0) ), lhei=c(0.5, 5, 0.5 ),
            scale="none",trace="none",labRow=NA,
            hclustfun=function(c){hclust(c, method='mcquitty')},breaks=breaks)

#  heatmap.2(normalized.ratio.matrix,Colv=NA,col=greenred(50),margins=c(5,25),main=paste(name,"- renormalized"),
 #           labRow=allnames[,"gene_name"],
#            lmat=rbind( c(0, 3), c(2,1), c(4,0) ), lhei=c(0.5, 5, 0.5 ),scale="none",trace="none",cexCol=0.7,labRow=NA,
#            hclustfun=function(c){hclust(c, method='mcquitty')},breaks=breaks)

  min.max <- quantile(ratio.matrix,probs=c(0.005,0.995))
  sel <- apply(ratio.matrix<min.max[1] | ratio.matrix > min.max[2],1,any)
  if (sum(sel) > 2) {
    ratio.matrix2 <- ratio.matrix[sel,,drop=FALSE]
    heatmap.2(ratio.matrix2,Colv=NA,col=greenred(50),dendrogram="row",
              margins=c(5,25),main=paste(name,"- above fold change of",round(10^min.max[2],1),
                                "or below",round(10^min.max[1],1)),
              key=FALSE, keysize=1.0, symkey=FALSE, density.info='none',
              trace='none', colsep=1:10,
#              labRow=allnames[sel,"gene_name"],
#              sepcolor='white', sepwidth=0.025,
              labRow=rownames(ratio.matrix2),
              scale="none",cexRow=1.5,cexCol=1,
              labCol = colnames(ratio.matrix),                 
              hclustfun=function(c){hclust(c, method='mcquitty')},
              lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),breaks=breaks)
  }
  dev.off()

  library(boot)
  cols <- colnames(ratio.matrix)
  correlation.matrix <- sapply(cols ,function(i) sapply(cols,function(j) boot::corr(ratio.matrix[,c(i,j)],
                                                            w=apply(1/variance.matrix,1,median))))

  pdf(sprintf("heatmap_correlation_%s.pdf",name),width=7,height=7,title=name)
  heatmap.2(correlation.matrix,scale="none",Rowv=NA,Colv=NA,col=greenred(50),
            density.info='none',trace="none",cellnote=round(correlation.matrix,1),main="not re")
  dev.off()

}

.scatter.plot <- function(my.table,name) {

  png(sprintf("scatter_plot%s.png",name),width=7,height=7)
  ggplot(my.table) + 
     geom_point(aes(x=lratio.x,y=lratio.y,color=comp)) + 
     facet_wrap(~comp)

}

.plot.pairs <- function(ratio.matrix,variance.matrix,name) {
  require(RColorBrewer)
  weights <- apply(1/variance.matrix,1,median)
  weights <- weights/sum(weights)
  
  lim <- max(abs(ratio.matrix))
  ##pdf("pairwise_correlation.pdf",title="Pairwise Correlation plot",width=10,height=10)
  png(sprintf("pairwise_correlation_%s.png",name),title="Pairwise Correlation plot",
      width=1000,height=1000,pointsize=14)
  pairs(ratio.matrix,xlim=c(-lim,lim),ylim=c(-lim,lim),
        text.panel=.panel.txt,
        upper.panel=.panel.smooth,weights=weights,
        lower.panel=function(...) .panel.cor(...,cex.cor=1.2),pch=1,cex=2)
  dev.off()

}

.panel.cor <- function(x, y, digits=2, prefix="", cex.cor, weights, ...)  {
  library(boot)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(boot::corr(matrix(c(x,y),ncol=2),w=weights))
  ##    r <- abs(cor(x,y,method="spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * sqrt(r) * 2)
}

.panel.smooth <- function(x, y, digits=2, prefix="", cex=1, weights, ...) {
  abline(c(0,0),1,col="red")
  points(x,y,cex=weights/max(weights)/2,...)
}

.panel.txt <- function(x,y,labels,cex,font,...) {
  ##abc<-c("1"="rep 1","2"="rep 2");
  
  s <- strsplit(labels,".",fixed=TRUE)[[1]]
  u <- par('usr')
  names(u) <- c("xleft", "xright", "ybottom", "ytop")
  pal <- brewer.pal(8,"Pastel1")
  bgcolor <- pal[as.numeric(s[[2]])]
  ##bgcolor <- pal[as.numeric(s[[3]])-114]
  do.call(rect, c(col = bgcolor, as.list(u)))
  ##text(0.5,0.5, paste(abc[s[[2]]],"\n",s[[3]],"/114",sep=""),cex=cex)
  text(0.5,0.5, paste(s[[2]],sep=""),cex=cex)
  #text(0.5,0.5, paste(s[[2]],"\n vs ",s[[3]],sep=""),cex=cex)
}

### depr plots
#  g <- ggplot(merged.table,aes(x=lratio,y=-log10(p.value.rat))) + 
#    geom_point(aes(color=factor(comp)),alpha=0.8) + 
#    facet_wrap(~comp,ncol=2) + 
#    geom_rug(alpha=0.5) + 
#    geom_hline(yintercept=-log10(0.05)) + 
#    geom_text(data=ddply(merged.table,"comp",
#                         function(x) c(x=Inf,y=-log10(0.05),isobar:::.sum.bool.na(x$p.value.rat < 0.05))),
#              aes(x=x,y=y,label=`TRUE`),vjust=-0.5,hjust=1)

#  pdf("lratio-pval.pdf")
#  print(g)
#  dev.off()

  #rownames(pg.df) <- do.call(paste,pg.df[,ac.vars,drop=FALSE])


