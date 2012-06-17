calc.zscore <- function(qq) {
  s.median <- median(qq,na.rm=TRUE)
  s.mad <- mad(qq,na.rm=TRUE)
  (qq-s.median)/s.mad
}

get.merged.table <- function(samples,cols=c("ac","r1","r2","lratio","variance"),merge.by=c("ac","r1","r2"),format="wide") {
  quant.tables <- lapply(seq(along=samples),
                         function(idx) {
                           load(paste(samples[idx],"/cache/quant.tbl.rda",sep=""))
                           q <- quant.tbl[,cols]
                           q <- ddply(q,c("class1","class2"),function(x) 
                                 cbind(x,zscore=round(calc.zscore(x[,"lratio"]),4)))
                           q[,"lratio"] <- round(q[,"lratio"],4)
                           q[,"variance"] <- round(q[,"variance"],4)
                           if (format=="wide") {
                             sel <- !colnames(q) %in% merge.by
                             colnames(q)[sel] <- paste(colnames(q)[sel],samples[idx],sep=".")
                           } else {
                             q$sample <- samples[idx]
                           }
                           message(paste(colnames(q),collapse=":"))
                           q
                         })
  if (format == "wide") {
    merged.table <- quant.tables[[1]]
    for (idx in  2:length(samples))
      merged.table <- merge(merged.table,
                            quant.tables[[idx]],by=merge.by,all=TRUE)
    merged.table$lratio <- rowMeans(merged.table[,.grep_columns(merged.table,"lratio")],na.rm=TRUE)
    merged.table$variance <- rowMeans(merged.table[,.grep_columns(merged.table,"variance")],na.rm=TRUE)
  } else {
    merged.table <- do.call(rbind,quant.tables)
  }
  return(merged.table)
}

get.names <- function(p,protein.group) 
  apply(my.protein.info(protein.group,p)[,c("name","gene_name","protein_name")],2,
        function(s) paste(unique(sort(gsub("'","",s))),collapse=", "))

calc.col <- function(tbl.wide,tag,add.name=TRUE) {
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

write.t <- function(...,sep="\t",row.names=FALSE,quote=FALSE)
  write.table(...,sep=sep,row.names=row.names,quote=quote)

write.summarized.table <- function(tbl.wide,all.names,cols) {
  summarized.table <- cbind(ac=rownames(tbl.wide),all.names)
  for (col in cols)
    summarized.table <- cbind(summarized.table,calc.col(tbl.wide,col,add.name=length(cols)>1))
  summarized.table <- cbind(summarized.table,tbl.wide)
  write.t(summarized.table,file=paste(name,"_summarized_table.csv",sep=""))
  return(summarized.table)
}

plot.heatmaps.gd <- function(ratio.matrix,name) {
  require(gplots)
  breaks <- seq(from=-max(abs(ratio.matrix)),to=max(abs(ratio.matrix)),length.out=51)
  pdf(sprintf("heatmap_%s.pdf",name),width=15,height=30,title=name)
  min.max <- quantile(ratio.matrix,probs=c(0.005,0.9995))
  sel <- apply(ratio.matrix<min.max[1] | ratio.matrix > min.max[2],1,any)
  if (sum(sel) > 2) {
    ratio.matrix2 <- ratio.matrix[sel,,drop=FALSE]
    heatmap.2(ratio.matrix2,Colv=NA,col=greenred(50),dendrogram="row",
              margins=c(5,25),main=paste(name,"- above fold change of",round(10^min.max[2],1),
                                "or below",round(10^min.max[1],1)),
              key=FALSE, keysize=1.0, symkey=FALSE, density.info='none',
              trace='none', colsep=1:10,
              #labRow=all.names[sel,"gene_name"],
              sepcolor='white', sepwidth=0.025,
              scale="none",cexRow=2,cexCol=2,
              labCol = colnames(ratio.matrix),                 
              hclustfun=function(c){hclust(c, method='mcquitty')},
              lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),breaks=breaks)
  }
  dev.off()

  cols <- colnames(ratio.matrix)
  correlation.matrix <- sapply(cols ,function(i) sapply(cols,function(j) corr(ratio.matrix[,c(i,j)],
                                                            w=apply(1/variance.matrix,1,median))))

  pdf(sprintf("heatmap_correlation_%s.pdf",name),width=15,height=30,title=name)
  heatmap.2(correlation.matrix,scale="none",Rowv=NA,Colv=NA,col=greenred(50),
            density.info='none',trace="none",cellnote=round(correlation.matrix,1),main="not re")
  dev.off()

}


plot.heatmaps <- function(ratio.matrix,name) {
  require(gplots)
  breaks <- seq(from=-max(abs(ratio.matrix)),to=max(abs(ratio.matrix)),length.out=51)
  pdf(sprintf("heatmap_%s.pdf",name),width=15,height=30,title=name)
  heatmap.2(ratio.matrix,Colv=NA,col=greenred(50),margins=c(5,25),main=paste(name),labRow=all.names[,"gene_name"],
            lmat=rbind( c(0, 3), c(2,1), c(4,0) ), lhei=c(0.5, 5, 0.5 ),scale="none",trace="none",
            hclustfun=function(c){hclust(c, method='mcquitty')},breaks=breaks)

  heatmap.2(normalized.ratio.matrix,Colv=NA,col=greenred(50),margins=c(5,25),main=paste(name,"- renormalized"),
            labRow=all.names[,"gene_name"],
            lmat=rbind( c(0, 3), c(2,1), c(4,0) ), lhei=c(0.5, 5, 0.5 ),scale="none",trace="none",
            hclustfun=function(c){hclust(c, method='mcquitty')},breaks=breaks)

  min.max <- quantile(ratio.matrix,probs=c(0.5,0.95))
  sel <- apply(ratio.matrix<min.max[1] | ratio.matrix > min.max[2],1,any)
  if (sum(sel) > 2) {
    ratio.matrix2 <- ratio.matrix[sel,,drop=FALSE]
    heatmap.2(ratio.matrix2,Colv=NA,col=greenred(50),dendrogram="row",
              margins=c(5,25),main=paste(name,"- above fold change of",round(10^min.max[2],1),
                                "or below",round(10^min.max[1],1)),
              key=FALSE, keysize=1.0, symkey=FALSE, density.info='none',
              trace='none', colsep=1:10,labRow=all.names[sel,"gene_name"],
              sepcolor='white', sepwidth=0.025,
              scale="none",cexRow=2,cexCol=2,
              labCol = colnames(ratio.matrix),                 
              hclustfun=function(c){hclust(c, method='mcquitty')},
              lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ),breaks=breaks)
  }
  dev.off()

  cols <- colnames(ratio.matrix)
  correlation.matrix <- sapply(cols ,function(i) sapply(cols,function(j) corr(ratio.matrix[,c(i,j)],
                                                            w=apply(1/variance.matrix,1,median))))

  pdf(sprintf("heatmap_correlation_%s.pdf",name),width=15,height=30,title=name)
  heatmap.2(correlation.matrix,scale="none",Rowv=NA,Colv=NA,col=greenred(50),
            density.info='none',trace="none",cellnote=round(correlation.matrix,1),main="not re")
  dev.off()

}

plot.pairs <- function(name) {
  require(RColorBrewer)
  weights <- apply(1/variance.matrix,1,median)
  weights <- weights/sum(weights)
  
  lim <- max(abs(ratio.matrix))
  ##pdf("pairwise_correlation.pdf",title="Pairwise Correlation plot",width=10,height=10)
  png(sprintf("pairwise_correlation_%s.png",name),title="Pairwise Correlation plot",
      width=1000,height=1000,pointsize=14)
  pairs(ratio.matrix,xlim=c(-lim,lim),ylim=c(-lim,lim),
        text.panel=panel.txt,
        upper.panel=panel.smooth,weights=weights,
        lower.panel=function(...) panel.cor(...,cex.cor=1.2),pch=1,cex=0.5)
  dev.off()

}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, weights, ...)  {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(corr(matrix(c(x,y),ncol=2),w=weights))
  ##    r <- abs(cor(x,y,method="spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * sqrt(r) * 2)
}

panel.smooth <- function(x, y, digits=2, prefix="", cex=1, weights, ...) {
  abline(c(0,0),1,col="red")
  points(x,y,cex=weights/max(weights)/2,...)
}

panel.txt <- function(x,y,labels,cex,font,...) {
  ##abc<-c("1"="rep 1","2"="rep 2");
  
  s <- strsplit(labels,".",fixed=TRUE)[[1]]
  u <- par('usr')
  names(u) <- c("xleft", "xright", "ybottom", "ytop")
  pal <- brewer.pal(8,"Pastel1")
  ##  bgcolor <- pal[as.numeric(s[[2]])]
  bgcolor <- pal[as.numeric(s[[3]])-114]
  do.call(rect, c(col = bgcolor, as.list(u)))
  ##text(0.5,0.5, paste(abc[s[[2]]],"\n",s[[3]],"/114",sep=""),cex=cex)
  text(0.5,0.5, paste(s[[2]],"\n vs ",s[[3]],sep=""),cex=cex)
}
