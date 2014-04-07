
# returns a list of all combinations of a group reporter and a group member
# where there are peptides which are only shared between reporter and member
.proteins.w.shared.peptides <- function(protein.group) {
  proteins.w.shared.pep = lapply(reporterProteins(protein.group), function(reporter.protein.g) {    
    gmp <- groupMemberPeptides(protein.group,reporter.protein.g,TRUE)           
    reporter.sp.sel <- gmp$peptide.info$specificity == REPORTERSPECIFIC      
    quant.sel <- gmp$peptide.info$n.shared.groups == 1 &                        
                 gmp$peptide.info$n.shared.proteins == 2                        
    gd.proteins <- apply(gmp$group.member.peptides[quant.sel,,drop=FALSE],2,any)         
                                                                                
    return(gd.proteins)                                                         
                                                                                
  })
  names(proteins.w.shared.pep) <- reporterProteins(protein.group)
  do.call(sum,proteins.w.shared.pep)

  proteins.w.shared.pep.dfl <- lapply(seq_len(length(proteins.w.shared.pep)), function(i) {
    gd.proteins <-proteins.w.shared.pep[[i]]
    if (any(gd.proteins))
      return(data.frame(reporter.protein.g=names(proteins.w.shared.pep)[i],
                        protein.g=names(gd.proteins)[gd.proteins]))
  })

  # remove entries w/o such cases
  proteins.w.shared.pep.dfl <- proteins.w.shared.pep.dfl[!sapply(proteins.w.shared.pep.dfl,is.null)]

  proteins.w.shared.pep.df <- do.call(rbind,proteins.w.shared.pep.dfl)
  proteins.w.shared.pep.df$reporter.protein.g <- as.character(proteins.w.shared.pep.df$reporter.protein.g)
  proteins.w.shared.pep.df$protein.g <- as.character(proteins.w.shared.pep.df$protein.g)

  return(proteins.w.shared.pep.df)
}

shared.ratios <- function(ibspectra,noise.model,channel1=NULL,channel2=NULL,protein=reporterProteins(proteinGroup(ibspectra)),...) {
  protein.group <- proteinGroup(ibspectra)
  l <- lapply(protein, function(reporter.protein.g) {    
    gmp <- groupMemberPeptides(protein.group,reporter.protein.g,TRUE)           
    reporter.sp.sel <- gmp$peptide.info$specificity == REPORTERSPECIFIC      
    quant.sel <- gmp$peptide.info$n.shared.groups == 1 &                        
                 gmp$peptide.info$n.shared.proteins == 2                        
    gd.proteins <- apply(gmp$group.member.peptides[quant.sel,,drop=FALSE],2,any)
    gd.proteins <- names(gd.proteins)[gd.proteins]
    gd.proteins <- gd.proteins[gd.proteins != reporter.protein.g]

    if (length(gd.proteins) > 0) {
      ratio.1 <- estimateRatio(ibspectra,noise.model=noise.model,               
                               protein=reporter.protein.g,channel1=channel1,      
                               channel2=channel2,...)
      if (is.na(ratio.1['lratio']))
        return(NULL)
      
      quant.df <- data.frame(reporter.protein=reporter.protein.g,
                             protein2=gd.proteins,
                             ratio1=ratio.1['lratio'],ratio1.var=ratio.1['variance'],
                             n.spectra.1=ratio.1['n.spectra'],
                             ratio2=NA,ratio2.var=NA,n.spectra.2=NA,                             
                             stringsAsFactors=FALSE)
      reporter.pep <- gmp$peptide.info$peptide[reporter.sp.sel]

      for (protein.i in seq_along(gd.proteins)) {
        protein.gd <- quant.df[protein.i,"protein2"]
        protein.pep <-  rownames(gmp$group.member.peptides)[gmp$group.member.peptides[,protein.gd]]
        quant.pep <- gmp$peptide.info[ quant.sel &
                                        gmp$peptide.info$peptide %in% protein.pep,
                                      "peptide"]
        ratio.2 <- estimateRatio(ibspectra,noise.model=noise.model,protein=NULL,
                                 peptide=quant.pep, channel1=channel1,channel2=channel2,combine=T,...)
        quant.df[protein.i,"ratio2"] <- ratio.2['lratio']
        quant.df[protein.i,"ratio2.var"] <- ratio.2['variance']
        quant.df[protein.i,"n.spectra.2"] <- ratio.2['n.spectra']
      }
      return(quant.df)
    } else {
      return(NULL)
    }                                                                                
  })
  l <- l[!sapply(l,is.null)]
  do.call(rbind,l)
}

shared.ratios.sign <- function(ress,z.shared,min.spectra=1,plot=TRUE) {
    ratio1.issmaller <- ress$ratio1+z.shared*sqrt(ress$ratio1.var) < ress$ratio2 - z.shared*(sqrt(ress$ratio2.var)) 
    ratio1.isbigger  <- ress$ratio1-z.shared*sqrt(ress$ratio1.var) > ress$ratio2 + z.shared*(sqrt(ress$ratio2.var))
 
    rr <- ress[(ratio1.issmaller | ratio1.isbigger) &
               !is.na(ratio1.issmaller) & !is.na(ratio1.isbigger) & 
               ress$n.spectra.1 >= min.spectra & ress$n.spectra.2 >= min.spectra,]

    rr<-rr[order(rr$ratio1,decreasing=T),]
    
    rr$proteins <- as.character(paste(rr$reporter.protein,"\nvs",rr$protein2))
    
    rr$n.spectra.txt <-"> 10"
    rr$n.spectra.txt[rr$n.spectra<=10]<-"6 - 10"
    rr$n.spectra.txt[rr$n.spectra<=5]<-"2 - 5"
    rr$n.spectra.txt[rr$n.spectra==1]<-"1"

    rr$n.spectra.txt.1 <-"> 10"
    rr$n.spectra.txt.1[rr$n.spectra.1<=10]<-"6 - 10"
    rr$n.spectra.txt.1[rr$n.spectra.1<=5]<-"2 - 5"
    rr$n.spectra.txt.1[rr$n.spectra.1==1]<-"1"
    xx=reshape(rr,list(c("ratio1","ratio2"),c("ratio1.var","ratio2.var"),
                       c("n.spectra.txt","n.spectra.txt.1")),
               direction="long",v.names=c("ratio","var","n.spectra"),timevar="g")
    xx$proteins=factor(as.character(xx$proteins),levels=rr$proteins,labels=rr$proteins)
#    xx$g <- as.character(xx$g)
#    xx$g[xx$g==1] <-"reporter"
#    xx$g[xx$g==2] <-"member"
    xx$g <- factor(xx$g,levels=c(1,2),labels=c("reporter","member"))
    if (min.spectra == 2)
      xx$n.spectra <- factor(as.character(xx$n.spectra),levels=c("2 - 5","6 - 10","> 10"))
    else
      xx$n.spectra <- factor(as.character(xx$n.spectra),levels=c("1","2 - 5","6 - 10","> 10"))

  if (plot) {
    require(ggplot2)
    breaks <- c(0.1,0.25,0.5,1,2,3,4,5,10,20,30,40,50)
    breaks <- breaks[log10(breaks) %inrange% range(xx$ratio)]

    print(ggplot(xx,aes_string(x="ratio",y="proteins")) +
          geom_vline(xintercept=0,alpha=0.5) + 
          geom_point(aes_string(colour="g",shape="g"),size=4) + 
          geom_errorbarh(aes(xmax=ratio+sqrt(var),xmin=ratio-sqrt(var),colour=factor(g),height=0.2)) + 
          scale_x_continuous("Ratio",breaks=log10(breaks),labels=breaks) +
          scale_colour_manual("group",values = c("blue","darkgreen")) + 
          scale_shape("group") + 
          scale_size(guide="none"))
  }
  xx
}

