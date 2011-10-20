
protein.acc <- function(prots,protein.info=NULL,ip=NULL) {
  if (is.null(ip)) {
    proteins <- list(prots)
  } else {
    proteins <- lapply(prots,function(p) {names(ip)[ip == p]})
  }

  sapply(proteins,function(prots) {
  pos <- grepl("-[0-9]*$",prots)
  if (sum(pos) == 0) {
    protl <- data.frame(protein=prots,accession=prots,splice=0)
  } else {
    if (sum(!pos) == 0) {
      protl <- data.frame(protein=prots[pos],
                      do.call(rbind,strsplit(prots[pos],"-")),
                      stringsAsFactors=FALSE,row.names=NULL)
    } else {
      protl <- data.frame(protein=c(prots[!pos],prots[pos]),
                      do.call(rbind,strsplit(c(paste(prots[!pos],"0",sep="-"),
                                               prots[pos]),"-")),
                      stringsAsFactors=FALSE,row.names=NULL)
    }
    colnames(protl) <- c("protein","accession","splice")
  }
  df <- protl
  df$splice <- as.numeric(df$splice)

  res <- ddply(df,"accession",function(y) {
      if(sum(y$splice>0) <= 1)
        return(data.frame(protein=unique(y$protein)))
      else 
        return(data.frame(protein=sprintf("%s-[%s]",
                            unique(y$accession),
                            paste(sort(y[y$splice>0,'splice']),collapse=","))))
    })
  return(paste(res$protein,collapse=", "))
})
}

protein.name <- function(prots,protein.info,href=function(x,y=x) {y},bold=FALSE,
                         name_col="protein_name",show.ac=TRUE,genename.first=FALSE,
                         tex=TRUE) {
  pos <- grepl("-",prots)
  if (sum(pos) == 0) {
    protl <- data.frame(protein=prots,accession=prots,splice=0)
  } else {
    if (sum(!pos) == 0) {
      protl <- data.frame(protein=prots[pos],
                      do.call(rbind,strsplit(prots[pos],"-")),
                      stringsAsFactors=FALSE,row.names=NULL)
      colnames(protl) <- c("protein","accession","splice")
    } else {
      protl <- data.frame(protein=c(prots[!pos],prots[pos]),
                      do.call(rbind,strsplit(c(paste(prots[!pos],"0",sep="-"),
                                               prots[pos]),"-")),
                      stringsAsFactors=FALSE,row.names=NULL)
      colnames(protl) <- c("protein","accession","splice")
    }
  }
  df <- protl
  df$splice <- as.numeric(df$splice)

  if (is.null(name_col) || nrow(protein.info) == 0 ||
      sum(protein.info$accession %in% df$accession)==0) {
    merge.ac.res <- ddply(df,"accession",function(y) {
      if(sum(y$splice>0) <= 1)
        return(data.frame(protein=href(y$protein)))
      else 
        return(data.frame(protein=href(unique(y$accession),
                            sprintf("%s%s-[%s]",
                                    unique(y$accession),
                                    ifelse(tex,"",""),
                                    paste(y[y$splice>0,'splice'],collapse=",")))))
    })
    merge.name.res <- data.frame(p=paste(merge.ac.res$protein,collapse=", "))

  } else {
    df <- merge(df, protein.info,all.x=TRUE,all.y=FALSE)

    merge.name.res <- ddply(df,name_col,function(x) {
      merge.ac.res <- ddply(x,"accession",function(y) {
        if(sum(y$splice>0) <= 1)
          return(data.frame(protein=href(y$protein)))
        else 
          return(data.frame(protein=href(unique(y$accession),
                              sprintf("%s%s-[%s]",
                                      unique(y$accession),
                                      ifelse(tex,"",""),
                                      paste(y[y$splice>0,'splice'],collapse=",")))))
      })
      p <- paste(merge.ac.res$protein,collapse=", ")
      
      name <- x[,name_col]
      if (name_col=="protein_name") {
        sel <- !is.na(x$gene_name) & !is.null(x$gene_name)  & x$gene_name!=""
        if (genename.first)
          name[sel] <- paste(x[sel,"gene_name"],": ",x[sel,"protein_name"],sep="")
        else
          name[sel] <- paste(x[sel,"protein_name"]," (",x[sel,"gene_name"],")",sep="")
      }
      
      name <- unique(name)
      
      if (length(name) > 0 & !is.na(name) & unique(name) != "") {
        if (bold)
          name <- sprintf("\\textbf{%s}",name)
        
        if (show.ac) 
          p <- sprintf("%s: %s",name,p)
        else 
          p <- name
      }
      
      data.frame(p=p)  
    })
  }
  return(paste(merge.name.res$p,collapse=", "))
}
 

draw.ratio <- function(ratio,sd,bnd,is.first=FALSE) {
  d <- new("Norm",ratio,sd)
  ci95.lower <- qnorm(0.025,ratio,sd)
  ci95.upper <- qnorm(0.975,ratio,sd)
  sprintf("\\boxplot{%.2f}{%.2f}{%.2f}{%.2f}{%.2f}{%.2f}{%s!%s}{%s}{%s}",
          .bnds(ratio,bnd),sd,.bnds(distr::q(d)(0.25),bnd),.bnds(distr::q(d)(0.025),bnd),
          .bnds(distr::q(d)(0.75),bnd),.bnds(distr::q(d)(0.975),bnd),
          ifelse(ratio > 0,"green","red"),
          min(floor(abs(ratio)/bnd*100),100),
          ifelse(distr::q(d)(0.25) < -bnd, "\\divcol", "black!1"),
          ifelse(distr::q(d)(0.75) > bnd, "\\divcol", "black!1"))
}



remove.comma <- function(s) {
	gsub(",","",s)
}

strsplit.n <- function(str,pattern=" ",max.length=10,ignore.href=TRUE) {
  lapply(strsplit(str,pattern),			
    function(s) {
      xx <- c(s[1])
      xx_i <- 1
      for (ss in s[2:length(s)]) {
        if (ignore.href)
          nch <- nchar(gsub("\\\\href\\{[^\\}]*\\}\\{([^\\}]*)\\}","\\1",xx[xx_i]))
        else
          nch <- nchar(xx[xx_i])
        if (nch <= max.length) {
          xx[xx_i] <- paste(xx[xx_i],ss,sep=pattern)
        } else {
          xx_i <- xx_i + 1
          xx[xx_i] <- ss
        }
      }
      return(xx)
    }
  )
}

paste.m <- function(x) {
  lapply(x,function(xx) paste("{",paste(xx,collapse=" \\\\ "),"}"))
}

draw.protein.ratio2 <- function(protein.tbl,bnd,diff=FALSE,ratiodistr=NULL) {
  cat("\\setcounter{nplot}{0}\n")
  cat("\\begin{tikzpicture}[x=2cm,line width=0.5pt,node distance=0.1]\n")
  
  for (p_i in seq_len(nrow(protein.tbl))){
    
    er <- protein.tbl[p_i,]
    if (is.na(er['lratio']) | is.na(er['variance'])) {
      er['lratio'] <- 0
      er['variance'] <- 0
      add <- "x"
      col <- "white"
      #col <- "black!60"
    } else {
      add <- " "
      col <- "black!60"
    }
    
    if (!is.null(ratiodistr) & length(er['is.significant']) == 1 &
        !is.na(er['is.significant']) & er['is.significant'] == 1)
      add <- "*"

    ratio.smaller.bnd <- er['lratio'] - sqrt(er['variance']) < -bnd
    ratio.bigger.bnd <- er['lratio'] + sqrt(er['variance']) > bnd
    if (!diff) {
      cat(sprintf("\\rsdboxplot{%.2f}{%.2f}{%s!%s}{%s}{%s}{%s}{%s}{%.0f}{%s}\n",
                  .bnds(er['lratio'],bnd),
                  .bnds(sqrt(er['variance']),bnd),
                  ifelse(er['lratio'] > 0,"green","red"),
                  min(floor(abs(er['lratio'])/bnd*100),100),
                  ifelse(ratio.smaller.bnd,"\\divcol", "black!1"),
                  ifelse(ratio.bigger.bnd, "\\divcol", "black!1"),
                  sprintf("%s %s/%s:",add,er['r1'],er['r2']),
                  ifelse(add=="x","NA",round(10^er['lratio'],2)),
                  10^(sqrt(er['variance'])),
                  col))
    }else {
      cat(sprintf("\\rsdboxplotnoname{%.2f}{%.2f}{%s!%s}{%s}{%s}{%.0f}{%s}\n",
                  .bnds(er['lratio'],bnd),
                  .bnds(sqrt(er['variance']),bnd),
                  ifelse(er['lratio'] > 0,"green","red"),
                  min(floor(abs(er['lratio'])/bnd*100),100),
                  ifelse(ratio.smaller.bnd, "\\divcol", "black!1"),
                  ifelse(ratio.bigger.bnd,  "\\divcol", "black!1"),
                  10^(sqrt(er['variance'])),
                  col))       
    }
  }
  cat("\\end{tikzpicture}\n")
}

draw.boxplot <- function(lratio,variance,bnd) {
  
  if (is.na(lratio) | is.na(variance)) {
    return("")
  }
  col <- "black!60"

  ratio.smaller.bnd <- lratio - sqrt(variance) < -bnd
  ratio.bigger.bnd  <- lratio + sqrt(variance) > bnd

  return(sprintf("\\boxplot{%.2f}{%.2f}{%s!%s}{%s}{%s}{%.0f}{%s}\n",
                  .bnds(lratio,bnd),
                  .bnds(sqrt(variance),bnd),
                  ifelse(lratio > 0,"green","red"),
                  min(floor(abs(lratio)/bnd*100),100),
                  ifelse(ratio.smaller.bnd, "\\divcol", "black!1"),
                  ifelse(ratio.bigger.bnd,  "\\divcol", "black!1"),
                  10^(sqrt(variance)),
                  col))
}

draw.protein.ratio <- function(ibspectra,protein,noise.model,relative.to,bnd,ratiodistr=NULL,...) {
  
  df <- data.frame()
  for (channel in setdiff(ibspectra@reporterNames,relative.to)) {
    er <- estimateRatio(ibspectra,noise.model=noise.model,
                        channel1=relative.to,channel2=channel,protein=protein,
                        ratiodistr=ratiodistr)
    er['r1'] <- relative.to
	 er['r2'] <- channel
	 df <- rbind(df,er)
  }
  draw.protein.ratio2(df,bnd,ratiodistr=ratiodistr,...)
}


draw.proteingroup2 <- function(pg,protein,protein.info,href=function(x,y=x) {y}) {
  ps <- peptideSpecificity(pg)
  sp <- spectrumToPeptide(pg)
  ip <- indistinguishableProteins(pg)
  pgs <- proteinGroupTable(pg)

  speci <- rep(3,length(ps$peptide))
  speci[ps$specificity=="group-specific"] <-  2
  speci[ps$specificity=="reporter-specific"] <- 1
  names(speci) <- ps$peptide

  
  n.to.s <- c("1","\\#","*")
  proteins <- pgs$protein.g[pgs$reporter.protein==protein]
  pepnprots <- peptideNProtein(pg)
  pp <- as.data.frame(pepnprots[pepnprots[,'protein.g'] %in% proteins,,drop=FALSE],stringsAsFactors=FALSE)

  proteins <- unique(pep.groups$protein.g)
  pep.groups <- get.pep.group(pg,protein)
  groups <- unique(pep.groups$group)

  cat("\\begin{longtable}{",rep("r",length(proteins)),"rrr}\n",sep="")
  prot.names <- c()
  for (protein in proteins) {
    prots <- names(ip)[ip == protein]
    prot.names <- c(prot.names,
                    paste("\\rotatebox{90}{",protein.name(prots,protein.info,href,bold=FALSE),"}"))
  }
  cat(paste(prot.names,collapse=" &"))

  cat(" &  \\underline{Peptide} & \\underline{Speficity} & \\underline{Spectra} \\endhead \n")
#cat("\\hline \\endhead \n")
  for (group_i in seq_along(groups)) {
    for (peptide in pep.groups$peptide[pep.groups$group==groups[group_i]]) {
      pep.string <- rep(n.to.s[speci[peptide]],length(proteins))
      pep.string[!proteins %in% pp$protein[pp$peptide==peptide]] <- "~"

      cat(pep.string,sep=" & ")
      cat("& {\\tt \\small ",peptide,"}","&",
          ps$specificity[ps$peptide==peptide],"&",
          sum(sp == peptide)," \\\\ \n")
    }
  }
  cat("\\end{longtable}\n")

}

draw.proteingroup <- function(protein.group,protein,protein.info,href=function(x,y=x) {y}) {
  ip <- indistinguishableProteins(protein.group)
  sp <- spectrumToPeptide(protein.group)
  ps <- peptideSpecificity(protein.group)
  group.table <- proteinGroupTable(protein.group)

  pep.groups <- get.pep.group(protein.group,protein)
  proteins <- unique(pep.groups$protein.g)

  protein.group.table <- group.table[group.table$reporter.protein==protein,]
  pepnprots <- peptideNProtein(protein.group)
  pp <- as.data.frame(pepnprots[pepnprots[,'protein.g'] %in%
                                protein.group.table$protein.g,,drop=FALSE],stringsAsFactors=FALSE)
  speci <- rep(3,length(ps$peptide))
  speci[ps$specificity=="group-specific"] <-  2
  speci[ps$specificity=="reporter-specific"] <- 1
  names(speci) <- ps$peptide

  n.to.s <- c("1","\\#","*")
 

#  cat("\\begin{longtable}[l]{llll}\n")


#  cat("\\end{longtable}\n")

  if (length(proteins)==1) {
    prots <- names(ip)[ip == proteins[1]]
	 cat(protein.name(prots,protein.info,href,bold=FALSE),"\n\n")
  } else {
	  draw.prot.groupmembers(pep.groups,proteins,reporter.g=protein,pp,ip,speci,href=href,protein.info=protein.info)
  }
  
  groups <- unique(pep.groups$group)
  one.group = length(groups) == 1
  cat("\\begin{longtable}[l]{rrrr}\n")
  cat(" & \\underline{peptide}",
      " & \\underline{speficity}",
      " & \\underline{spectra}",
      " \\\\ \n")
  cat(" \\endhead \n")
  for (group_i in seq_along(groups)) {
    for (peptide in pep.groups$peptide[pep.groups$group==groups[group_i]]) {
      cat(ifelse(one.group,"",group_i),"&","{\\tt \\small ",peptide,"}","&",
          ps$specificity[ps$peptide==peptide],"&",
          sum(sp == peptide,na.rm=TRUE)," \\\\ \n")
    }
  }
  cat("\\end{longtable}\n")
}


sanitize.ref <- function(str) {
	result <- str
	result <- gsub("\\\\", "_", result)
	result <- gsub("$", "_", result, fixed = TRUE)
	result <- gsub("-", "_", result, fixed = TRUE)
	result <- gsub(">", "_", result, fixed = TRUE)
	result <- gsub("<", "_", result, fixed = TRUE)
	result <- gsub("|", "_", result, fixed = TRUE)
	result <- gsub("{", "_", result, fixed = TRUE)
	result <- gsub("}", "_", result, fixed = TRUE)
	result <- gsub("%", "_", result, fixed = TRUE)
	result <- gsub("&", "_", result, fixed = TRUE)
	result <- gsub("_", "_", result, fixed = TRUE)
	result <- gsub("#", "_", result, fixed = TRUE)
	result <- gsub("^", "_", result, fixed = TRUE)
	result <- gsub("~", "_", result, fixed = TRUE)
	return(result)
}



draw.prot.groupmembers <- function(pep.groups,proteins,reporter.g,
                                   pp,ip,speci,href,protein.info) {  
  n.to.s <- c("1","\\#","*")
  groups <- unique(pep.groups$group)
  cat("\\begin{tabular*}{\\textwidth}{p{2in} r r ",rep("r",length(groups)),"r}\n")
  cat("\\underline{protein} & \\underline{\\#p}  &\\underline{}  & \\\\ \n")
  for (protein in proteins) {
    prots <- names(ip)[ip == protein]
	  name <- protein.name(prots,protein.info,href,bold=FALSE)
    #cat(" |[text width=7cm]| ",paste.m(strsplit.n(name,max.length=40))[[1]]," & ",
    #    "|text width=1cm|",sum(pp$protein.g==protein)," & \n")
    cat(name," & ",sum(pp$protein.g==protein)," & \n")

 	  # draw peptide coverage for each group
	  for (group_i in seq_along(groups)) {
      group <- groups[group_i]
      d <- pep.groups[pep.groups$group==group,]
      pep.string <- n.to.s[speci[d$pept]]
      pep.string[!d$pept %in% pp$peptide[pp$protein==protein]] <- "~"
		group <- unique(d$group)

      pep.string <-
        do.call(c,lapply(c(n.to.s,"~"),function(x) {
          if(sum(pep.string == x)>10) {
            res <- paste(x,x,x,"\\dots",x,sep="~")
            if (x=="~") {
              "~~~~~~~~~~~"
            } else{
              res
            }
          }else pep.string[pep.string==x] } ))
      
      cat("\\tikznode{",remove.comma(protein),"-",which(groups==group),
          "}{\\texttt{",paste(pep.string,collapse="~"),"}} &\n",sep="")
      data.frame()
    } 
    cat(" \\\\ \n")
  }
  cat("\\end{tabular*}\n")
  one.group = length(groups) == 1
  if (!one.group) {
   cat("\\begin{tikzpicture}[overlay]\n")
   for (group_i in seq_along(groups)) {
     cat("\\drawbrace{",remove.comma(protein),"-",group_i,"}{",group_i,"};\n",sep="")
   }
   cat("\\end{tikzpicture}\n")
  }
  cat("\n") 
}

#multirow function
mr <- function(text,width="*") {
  if (ncol(cmbn) > 1) {
    return(sprintf("\\multirow{%s}{%s}{%s}",
                   ncol(cmbn),width,
                   text))
  } else {
    return(text)
  }
}

print.proteinrow <- function(protein.idx,peptides.sel) {
  if (length(peptides.sel)>0) {
      cat(sprintf("  \\proteinrow{%s}{%s}{}\n",
                  protein.idx,paste(peptides.sel,collapse=",")))
  }
}

connect.nodes <- function(protein.idx,pep.pos) {
   if (length(pep.pos) == 1) return();
#   message(protein.idx," (",length(pep.pos),"):   ",paste(pep.pos,collapse=","),"\n")
   last_pos = pep.pos[1];
   cat(sprintf("  \\draw (prot%s pep%s)",protein.idx,last_pos))
   for (i in 2:length(pep.pos)) {
      if (pep.pos[i] == (last_pos+1)) {
         cat(" --")
      }
      last_pos <- pep.pos[i]
      cat(sprintf(" (prot%s pep%s)",protein.idx,last_pos))
   }
   cat(";\n")
}

tikz.proteingroup <- function(protein.group,reporter.protein.g,show.pos) {

  gmp <- groupMemberPeptides(protein.group,reporter.protein.g,TRUE)
  reporter.sp.sel <- gmp$peptide.info$specificity == "reporter-specific"
  quant.sel <- gmp$peptide.info$n.shared.groups == 1 & 
               gmp$peptide.info$n.shared.proteins == 2
  group.sp.sel <- gmp$peptide.info$specificity == "group-specific" & !quant.sel
  unspecific.sel <- gmp$peptide.info$specificity == "unspecific"
  
  peptide.styles <- rep("us",nrow(gmp$peptide.info))
  peptide.styles[reporter.sp.sel] <- "rs"
  peptide.styles[group.sp.sel] <- "gs"

  gd.proteins <- apply(gmp$group.member.peptides[quant.sel,,drop=FALSE],2,any)
  col.i <- 1
  for (quant.prot.g in names(gd.proteins)[gd.proteins]) {
    if (quant.prot.g != reporter.protein.g) {
      quant.p.sel <- quant.sel & gmp$group.member.peptides[,quant.prot.g]
      peptide.styles[quant.p.sel] <- sprintf("quant peptide %s",colors[col.i])
      if (col.i == length(colors)) {
        col.i <- 1
      } else {
        col.i <- col.i + 1
      }
    }
  }

  n.peptides <- length(reporter.sp.sel)
  max.n.peptides <- 10
  tikz.y <- 0.5
  tikz.x <- 0.5
  tikz.minnodesize <- 0.2

  if (n.peptides > max.n.peptides) {
    tikz.x <- round(tikz.x*max.n.peptides/n.peptides,2)
    tikz.minnodesize <- round(tikz.minnodesize*max.n.peptides/n.peptides,2)
  }
 
  cat(sprintf("\\begin{tikzpicture}[x=%scm,y=%scm,every node/.style={minimum size=%scm}]\n",tikz.x,tikz.y,tikz.minnodesize))
  cat("  \\node at (1,0)[anchor=west] {peptides};\n")

  n.groupmember <- ncol(gmp$group.member.peptides)
  protein.peptides.df <- data.frame(
    i=seq_len(n.groupmember),
    rs=0,
    gs=0,
    us=0)

  for (protein.i in seq_len(n.groupmember)) {
    x = colnames(gmp$group.member.peptides)[protein.i]

    my.protein.info <- my.protein.info(protein.group,x)
    human.protein.name <- human.protein.names(my.protein.info)
    
    name.i <- 1
    table.name.nolink <- human.protein.name$ac_nolink[1]
    nch.tablename <- nchar(table.name.nolink)
    while (name.i < length(human.protein.name$ac_link)) {
      nch.tablename <- nch.tablename + nchar(human.protein.name$ac_nolink[name.i])
      if (nch.tablename < 20) {
        table.name.nolink <- paste(table.name.nolink,
                             human.protein.name$ac_nolink[name.i],sep=",")
      } else {
        table.name.nolink <- paste(table.name.nolink,"\\dots",sep=",")
        name.i <- Inf
      }
      name.i <- name.i + 1
    }

    peptide.idx.sel <- gmp$group.member.peptides[,protein.i]
    protein.peptides.df[protein.i,] <- c(protein.i,
                sum(peptide.idx.sel&reporter.sp.sel),
                sum(peptide.idx.sel&(group.sp.sel|quant.sel)),
                sum(peptide.idx.sel&unspecific.sel))
   
     print.proteinrow(
        -protein.i,paste(which(peptide.idx.sel),peptide.styles[peptide.idx.sel],sep="/")
     )
   
     if (n.peptides < 15) {
       connect.nodes(-protein.i,which((reporter.sp.sel | group.sp.sel | unspecific.sel | quant.sel)&gmp$group.member.peptides[,protein.i]))
     }
  }

  cat("  \\matrix[ampersand replacement=\\&,matrix anchor=us node.east,anchor=base] at (0,0) {\n",
      "    \\node{}; \\& \\node {rs};  \\& \\node {gs}; \\& \\node (us node){us};  \\\\\n");
  for (protein.i in seq_len(n.groupmember)) {
     if (show.pos) cat(sprintf("    \\node[draw,gray]{%s};",protein.i))
     cat(sprintf("    \\& \\node{%s};\\& \\node{%s};\\& \\node{%s}; \\\\\n",
                 protein.peptides.df[protein.i,"rs"],
                 protein.peptides.df[protein.i,"gs"],
                 protein.peptides.df[protein.i,"us"]))
  }
  cat("  };\n")
  cat("\\end{tikzpicture}\n")

}

draw.protein.group <- function(protein.group,reporter.protein.g) {


  gmp <- groupMemberPeptides(protein.group,reporter.protein.g,TRUE)
  pgt <- proteinGroupTable(protein.group)
  show.pos <- ncol(gmp$group.member.peptides) > 1
# cat("\\paragraph{Proteins}\n")
  if (show.pos) {
    cat("\\begin{tabular}{rl@{}llp{5cm}l}\n")
    cat(paste("pos","accession","","gene name","protein name",sep=" & "))
  } else {
    cat("\\begin{tabular}{l@{}llp{5.5cm}l}\n")
    cat(paste("accession","","gene name","protein name",sep=" & "))
  }
  cat(" & \\multirow{",length(unique(protein.ac(protein.group,
                                                pgt[pgt$reporter.protein==reporter.protein.g,"protein.g"])))+1,
      "}{*}{%\n",sep="")
  tikz.proteingroup(protein.group,reporter.protein.g,show.pos)
  cat("} \\\\ \n")
  for (protein.i in seq_len(ncol(gmp$group.member.peptides))) {
    x = colnames(gmp$group.member.peptides)[protein.i]

    my.protein.info <- my.protein.info(protein.group,x)
    for (ac in unique(my.protein.info$accession)) {
      sel <- my.protein.info$accession == ac
      var.string <- number.ranges(my.protein.info$splicevariant[sel])
      if (show.pos) cat(protein.i,"&")
      cat(paste(sprintf("\\uniprotlink{%s}",sanitize(ac)),
          ifelse(is.na(var.string),"",var.string),
          sanitize(unique(my.protein.info$gene_name[sel])),
          sanitize(unique(my.protein.info$protein_name[sel])),"",
          sep=" & "),"\\\\ \n")
          
    }
    human.protein.name <- human.protein.names(my.protein.info)
  } 
  cat("\\end{tabular}\n")

}

print.longtablehdr <- function(coldef,draw.channels,ncol.p,draw.signcol) {
  #cat("\n\n\\renewcommand{\\arraystretch}{0.75}\n")
  cat("\\begin{longtable}{",coldef,"}",'\n',sep="")
  cat("\t \\# ",
      " & \\textbf{protein}",
      " & \\rh{group}",
      " & \\rh{peptides}",
      " & \\rh{spectra}",sep="\n\t")
  
  # no ch1/ch2 columns when only one comparision available
  if (draw.channels) {
    cat('\n')
    cat(" & \\rh{ch1}",
        " & \\rh{ch2}",sep="\n\t")
  }
  cat("\n\t"," & \\rh{quant}")
  cat("\n\t"," & \\textbf{ratio}")
  if (draw.signcol) {
    cat("\n\t"," & ")
    #cat("\n\t"," & \\textbf{*}")
  }
  cat("\n\t"," & {\\hfill \\drawaxis{3pt}{south}}",
      "\n \\endhead \n",rep(" &",  ncol.p-1),
      " {\\hfill \\drawaxis{-3pt}{north}}",
      "\n \\endfoot \n")
}
