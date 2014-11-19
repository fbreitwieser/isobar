#################################################
## Helper Functions for LaTeX Report Generation

load.tex.properties <- function(env) {
  get.property <- function(name) { get(name,env$properties.env,inherits=FALSE) }
  if (!exists("properties.env",envir=env)) {
    warning("No properties.env - loading properties")
    assign("properties.env",load.properties(),envir=env)
  }
  needed.objects <- c('ibspectra','noise.model','quant.tbl','ratiodistr')
  if (!all(needed.objects %in% objects(envir=env))) {
    warning("Not all necessary objects present - calling initialize.env")
    initialize.env(env,env$properties.env) 
  }
  assign("get.property",get.property,envir=env)
}

write.tex.commands <- function() {
  get.property <- get("get.property")
  cat("\\newcommand{\\analysisname}{",sanitize(get.property('name'),dash=FALSE),"}\n")
  cat("\\newcommand{\\analysisauthor}{",sanitize(get.property('author'),dash=FALSE),"}\n")

  cat("\\newcommand{\\isobarthanks}{\\thanks{This report was 
    generated using the \\texttt{isobar} R package version ",
    packageDescription("isobar")$Version,
    " [built using ",packageDescription("isobar")$Built,"]",
    ". If you use it in published work, please cite ",
    "'Breitwieser FP \\textit{et~al.}: General statistical ",
    " modeling of data from protein relative expression isobaric tags, ",
    " \\textit{Journal of Proteome Research} 2011' and ",
    "'Breitwieser FP and Colinge J: isobarPTM: A software tool for the quantitative analysis of post-translationally modified proteins,",
    "\\textit{Journal of Proteomics} 2013","'}}\n",sep="")
}

print_longtablehdr <- function(level,is.single.comparision,is.quant.tbl,file="") {

  draw.channels <- !is.single.comparision
  draw.signcol <- is.quant.tbl

  if (is.single.comparision) {
    coldef <- "@{}rXfrrrr@{}r@{}r@{}"
    coldef.s <- "@{}rXfrrrr@{}r@{}"
    ncol.p <- 9
  } else {
    coldef <- "@{}rXfrrrrrrr@{}r@{}"
    coldef.s <- "@{}rXfrrrrrrr@{}"
    ncol.p <- 11
  }

  if (!is.quant.tbl) {
    coldef <- coldef.s 
    ncol.p <- ncol.p -1
  }


  mycat <- function(...,append=TRUE,sep="") 
    cat(...,file=file,append=append,sep=sep)

  mycat("\\begin{longtable}{",coldef,"}\n",
        "  \\# ",paste(" & \\textbf{",level,"}")," & \\rh{group}",
        " & \\rh{peptides}"," & \\rh{spectra}",append=FALSE)
  
  # no ch1/ch2 columns when only one comparision available
  if (draw.channels) 
    mycat("\n  & \\rh{ch1} & \\rh{ch2}")
  
  mycat("\n  & \\rh{quant} & \\textbf{ratio}")
  if (draw.signcol) 
    mycat(" & ")
  
  mycat("\n  & {\\hfill \\drawaxis{3pt}{south}}",
        "\n\\endhead \n",rep(" &",  ncol.p-1),
       " {\\hfill \\drawaxis{-3pt}{north}}",
       "\n \\endfoot \n")
}

print_longtablehdr_peptide <- function(coldef,draw.channels,ncol.p,draw.signcol) {
  #cat("\n\n\\renewcommand{\\arraystretch}{0.75}\n")
  cat("\\begin{longtable}{",coldef,"}",'\n',sep="")
  cat("\t \\# ",paste(" & \\textbf{peptide}"))  
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


draw.boxplot <- function(lratio,sd,bnd) { 
  if (is.na(lratio) || is.na(sd) || !is.finite(lratio) || !is.finite(sd)) {
    return("")
  }

  ratio.smaller.bnd <- lratio - sd < -bnd
  ratio.bigger.bnd  <- lratio + sd > bnd

  return(sprintf("\\boxplot{%.2f}{%.2f}{%s!%s}{%s}{%s}{%.0f}\n",
                  lratio,
                  sd,
                  ifelse(lratio > 0,"green","red"),
                  min(floor(abs(lratio)/bnd*100),100),
                  ifelse(ratio.smaller.bnd, "\\divcol", "black!1"),
                  ifelse(ratio.bigger.bnd,  "\\divcol", "black!1"),
                  sd))
}

transform_pepmodif <- function(pep.n.modif) {
  message(paste(pep.n.modif,collapse=" - "))
  pep <- strsplit(pep.n.modif[1],"")[[1]]
  modif <- strsplit(pep.n.modif[2],":")[[1]]
  if (length(pep)+1 != length(modif))
    stop("lengthwise i dont like them: ",pep.n.modif[1],
         " and ",pep.n.modif[2])
#  if (modif[length(modif)] != "")
#    stop("modif contains something in the end:",
#         pep.n.modif[2])
  
  modif <- modif[seq(from=2,to=length(modif))]
  pos <- as.numeric(sapply(modif,function(m) which(modifs[,1] == m)))

  pep[!is.na(pos)] <- paste("\\pdftooltip{\\textcolor{",modifs[pos[!is.na(pos)],3],
                           "}{\\underline{",pep[!is.na(pos)],"}}}{",
                            modifs[pos[!is.na(pos)],1],"}",sep="")
  return(paste(pep,collapse=""))
}

modifs <-
  matrix(c(
           "iTRAQ4plex_Nterm","i","black",
           "iTRAQ4plex_K","i","black",
           "Oxidation_M","o","black",
           "Cys_CAM","c","black",
           "PHOS","c","purple"),
         byrow=TRUE,ncol=3)

.bnds <- function(x,bnd=NULL,min.x=-bnd,max.x=bnd) {
 max(
      min(x,max.x),
      min.x)
}

# tex helper functions
draw.proteingroup.row <- function(name,protein.group,reporter.protein.g,file=file) {

  mycat <- function(...,append=TRUE,sep="") 
    cat(...,file=file,append=append,sep=sep)

  gmp <- groupMemberPeptides(protein.group,reporter.protein.g,TRUE)
  pgt <- proteinGroupTable(protein.group)

  n.multirow <- length(unique(protein.ac(protein.group,pgt[pgt$reporter.protein==reporter.protein.g,"protein.g"])))
  mycat(name,"&",sep=" ")
  n.row <- 1
  show.pos <- TRUE
  for (protein.i in seq_len(ncol(gmp$group.member.peptides))) {
    x = colnames(gmp$group.member.peptides)[protein.i]

    my.protein.info <- my.protein.info(protein.group,x)
    my.acs <- unique(my.protein.info$accession)
    for (ac.i in seq_along(my.acs)) {
      if (n.row > 1) mycat(" & ") 
      sel <- my.protein.info$accession == my.acs[ac.i]
      var.string <- number.ranges(my.protein.info$splicevariant[sel])
      #cat(protein.i,"&")
      if (show.pos && ncol(gmp$group.member.peptides)>1) {
        if (ac.i == 1)
          mycat("{\\small (",protein.i,")}",sep="")
        else
          mycat ("$\\mathord{\\cdot}$")
      }

      mycat(" & ")
      mycat(sprintf("\\uniprotlink{%s}%s",sanitize(my.acs[ac.i],dash=FALSE),ifelse(is.na(var.string),"",paste0("-",var.string)[1])),
          sanitize(unique(my.protein.info$gene_name[sel]))[1],
          sanitize(unique(my.protein.info$protein_name[sel]))[1],
          sep=", ")

      if (n.row==1) {
        # PEPTIDE GROUPING
        mycat(" & ")
        if (n.multirow > 1) mycat(" \\multirow{",n.multirow ,"}{*}{%\n",sep="")
        mycat(tikz.proteingroup(protein.group,reporter.protein.g,show.pos,show.header=FALSE))
        if (n.multirow > 1) mycat(" }")
        mycat(" \\\\ \n")
      } else {
        mycat(" & \\\\ \n")
      }
      n.row <- n.row + 1
          
    }
    #human.protein.name <- human.protein.names(my.protein.info)
  } 
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
      cat(paste(sprintf("\\uniprotlink{%s}",sanitize(ac,dash=FALSE)),
          ifelse(is.na(var.string),"",var.string),
          sanitize(unique(my.protein.info$gene_name[sel])),
          sanitize(unique(my.protein.info$protein_name[sel])),"",
          sep=" & "),"\\\\ \n")
          
    }
    #human.protein.name <- human.protein.names(my.protein.info)
  } 
  cat("\\end{tabular}\n")

}

tikz.proteingroup <- function(protein.group,reporter.protein.g,show.pos,show.header=TRUE,
                              colors=c("blue","red","green","cyan","magenta","yellow")) {

  gmp <- groupMemberPeptides(protein.group,reporter.protein.g,TRUE)
  reporter.sp.sel <- gmp$peptide.info$specificity == "reporter-specific"
  quant.sel <- gmp$peptide.info$n.shared.groups == 1 & 
               gmp$peptide.info$n.shared.proteins == 2
  group.sp.sel <- gmp$peptide.info$specificity == "group-specific" & !quant.sel
  unspecific.sel <- gmp$peptide.info$specificity == "unspecific"
  
  peptide.styles <- rep("us",nrow(gmp$peptide.info))
  peptide.styles[reporter.sp.sel] <- "rs"
  peptide.styles[group.sp.sel] <- "gs"

  res <- c()

  gd.proteins <- apply(gmp$group.member.peptides[quant.sel,,drop=FALSE],2,any)
  col.i <- 1
  for (quant.prot.g in names(gd.proteins)[gd.proteins]) {
    if (quant.prot.g != reporter.protein.g) {
      quant.p.sel <- quant.sel & gmp$group.member.peptides[,quant.prot.g]
      peptide.styles[quant.p.sel] <- sprintf("gs %s",colors[col.i])
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
  tikz.x <- 0.4
  tikz.minnodesize <- 0.2

  if (n.peptides > max.n.peptides) {
    tikz.x <- round(tikz.x*max.n.peptides/n.peptides,2)
    tikz.minnodesize <- round(tikz.minnodesize*max.n.peptides/n.peptides,2)
  }
  
  res <- paste(res,sprintf("\\begin{tikzpicture}[x=%scm,y=%scm,every node/.style={minimum size=%scm}]\n",tikz.x,tikz.y,tikz.minnodesize))
  if (show.header)
    res <- paste(res,"  \\node at (1,0)[anchor=west] {peptides};\n")

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
   
     res<- paste(res,.print.proteinrow(
        -protein.i,paste(which(peptide.idx.sel),peptide.styles[peptide.idx.sel],sep="/")
     ))
   
     if (n.peptides < 15) {
       res <- paste(res,connect.nodes(-protein.i,which((reporter.sp.sel | group.sp.sel | unspecific.sel | quant.sel)&gmp$group.member.peptides[,protein.i])))
     }
  }

  res <- paste(res,sprintf("  \\node [left=.25cm of prot-1 pep1] {%s/%s/%s};\n",
              sum(reporter.sp.sel),
              sum(group.sp.sel),
              sum(unspecific.sel)))


  if (show.header) {
    res <- paste(res,"\n","  \\matrix[ampersand replacement=\\&,matrix anchor=us node.east,anchor=base] at (0,0) {\n")
    res <- paste(res,"\n","    \\node{}; \\& \\node {rs};  \\& \\node {gs}; \\& \\node (us node){us};  \\\\\n")
    for (protein.i in seq_len(n.groupmember)) {
       if (show.pos) res <- paste(res,"\n",sprintf("    \\node[draw,gray]{%s};",protein.i))
       res <- paste(res,"\n",sprintf("    \\& \\node{%s};\\& \\node{%s};\\& \\node{%s}; \\\\\n",
                   protein.peptides.df[protein.i,"rs"],
                   protein.peptides.df[protein.i,"gs"],
                   protein.peptides.df[protein.i,"us"]))
    }
    res <- paste(res,"\n","  };\n")
  }
  res <- paste(res,"\\end{tikzpicture}\n")

  return(res)
}

.tex.combinenames <- function(protein_name,is.single.comparision,cmbn) {
  gene.names <- c()
  for (prot.gene.name in protein_name) {
    if (is.na(prot.gene.name)) next;
    prot.gene.name <- ifelse(nchar(prot.gene.name)>80,
                             paste(sanitize(substr(prot.gene.name,0,76)),"\\dots"),
                             sanitize(prot.gene.name))
    gene.names <- c(gene.names,prot.gene.name)
  }
  if (!is.single.comparision) {
    if (length(gene.names) > ncol(cmbn) - 1) {
      gene.names <- gene.names[seq_len(ncol(cmbn)-1)]
      gene.names[length(gene.names)] <- paste(gene.names[length(gene.names)],", \\dots",sep="")
    }
  }
  return(gene.names)
}

.print.proteinrow <- function(protein.idx,peptides.sel) {
  if (length(peptides.sel)>0) {
      sprintf("  \\proteinrow{%s}{%s}{}\n",
                  protein.idx,paste(peptides.sel,collapse=","))
  }
}

connect.nodes <- function(protein.idx,pep.pos) {
   if (length(pep.pos) == 1) return();
#   message(protein.idx," (",length(pep.pos),"):   ",paste(pep.pos,collapse=","),"\n")
   last_pos = pep.pos[1];
   res <- sprintf("  \\draw (prot%s pep%s)",protein.idx,last_pos)
   for (i in 2:length(pep.pos)) {
      if (pep.pos[i] == (last_pos+1)) {
         res <- paste(res,"--")
      }
      last_pos <- pep.pos[i]
      res <- paste(res,sprintf("(prot%s pep%s)",protein.idx,last_pos))
   }
   paste(res,";\n")
}


# function adapted from iquantitator
sanitize <- function(str,dash=TRUE) {
  result <- str
  result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
  result <- gsub("$", "\\$", result, fixed = TRUE)
  if (dash)
    result <- gsub("-", "\\nobreakdash-", result, fixed = TRUE)
  result <- gsub(">", "$>$", result, fixed = TRUE)
  result <- gsub("<", "$<$", result, fixed = TRUE)
  result <- gsub("|", "$|$", result, fixed = TRUE)
  result <- gsub("{", "\\{", result, fixed = TRUE)
  result <- gsub("}", "\\}", result, fixed = TRUE)
  result <- gsub("%", "\\%", result, fixed = TRUE)
  result <- gsub("&", "\\&", result, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result <- gsub("#", "\\#", result, fixed = TRUE)
  result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
  result <- gsub("~", "\\~{}", result, fixed = TRUE)
  result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", result, fixed = TRUE)
  return(result)
}


print_sign_proteins_tbl <- function(file,cmbn,protein.group,quant.tbl,my.protein.infos,bnd) {

  is.single.comparision <- ncol(cmbn) ==1
  mycat <- function(...,append=TRUE,sep="") 
    cat(...,file=file,append=append,sep=sep)

  print_longtablehdr("protein",is.single.comparision,FALSE,file=file)
  draw.hline <- FALSE
  
  for (cmb_i in seq_len(ncol(cmbn))) {
    sel <- quant.tbl[["r1"]]==cmbn[1,cmb_i] &
           quant.tbl[["r2"]]==cmbn[2,cmb_i] &
           quant.tbl[["is.significant"]] == 1
    sel[is.na(sel)] <- FALSE
    if (any(sel)) {
      if (draw.hline) {
        mycat (" \\hline\n")
      } else {
        draw.hline = TRUE
      }
      prot_i <- 1
  
      proteins <- quant.tbl[sel,"ac"][order(quant.tbl[sel,"lratio"])]
      for (protein in proteins) {
        prot.info <- my.protein.infos[[protein]]
        protein.row <- quant.tbl[quant.tbl[["ac"]]==protein & sel,,drop=FALSE]
        reporter.peptides <- peptides(protein.group,protein=protein)
        n.spectra <- length(names(spectrumToPeptide(protein.group))[spectrumToPeptide(protein.group)%in%reporter.peptides])

        protein.names <- prot.info[["collapsed.gene_name"]][["protein_name"]]
        gene.names <- .tex.combinenames(protein.names,is.single.comparision,cmbn)

        mycat(" ",sprintf("\\hyperref[protein.%s]{%s}",
                          protein.row[,"group"],prot_i),
              " & ",paste0(prot.info[["table.name"]],": ",
                           paste(gene.names,collapse=", ")),
              " & ",print_groupsize(prot.info[["n.reporter"]],prot.info[["n.groupmember"]]),
              " & ",length(reporter.peptides),
              " & ",n.spectra)
        if (!is.single.comparision) {
            mycat(" & ",cmbn[1,cmb_i]," & ",cmbn[2,cmb_i])
        }
        mycat(" & ",protein.row[1,"n.spectra"])
        mycat(" & ",sprintf("\\textbf{%.2f}",10^protein.row[1,"lratio"]),sep=" ")
        mycat(" & ",draw.boxplot(protein.row[1,'lratio'],
                               protein.row[1,'sd'],bnd))

        mycat(" \\\\  \n")
        prot_i <- prot_i + 1
      }
    }
  }
  mycat("\\end{longtable}\n")

}


print_groupsize <- function (r,t) {
  if (r>1 || t>1)
    paste0(r,"/",t)
}

get_n_proteins <- function(quant.tbl,sign=FALSE,is.na=FALSE) {
  n.quant <- ddply(quant.tbl,
                   c("class1","class2"),
                   function(x) {
                     sel <- !is.na(x[["lratio"]])
                     if (is.na) sel <- !sel
                     if (sign) sel <- sel & x[["is.significant"]]
                     sum(sel)
                   })
  n.quant <- n.quant[n.quant[["V1"]]>0,]

  paste(sprintf("%s/%s: %s",
                    sanitize(n.quant[["class2"]]),
                    sanitize(n.quant[["class1"]]),
                    n.quant[["V1"]]),collapse="; ")
}

print_protein_quant_tbl <- function(file="",
                                    cmbn,protein.group,
                                    quant.tbl,
                                    my.protein.infos,
                                    bnd) {

  message("Writing protein quantifications table ... ",append.lf=FALSE)

  is.single.comparision <- ncol(cmbn) ==1
  mycat <- function(...,append=TRUE,sep="") 
    cat(...,file=file,append=append,sep=sep)

  mr <- function(text,...) if (is.single.comparision) { text } else { sprintf("\\mr{%s}",text) }
  mrp <- function(text,...) if (is.single.comparision) { text } else {sprintf("\\mrp{%s}",text)}

  print_longtablehdr("protein",is.single.comparision,TRUE,file=file)

  proteins.n.group <- unique(quant.tbl[,c('ac','group')])
  proteins <- proteins.n.group[order(proteins.n.group[["group"]]),"ac"]
  
  for (protein in proteins) {
    protein.rows <- quant.tbl[quant.tbl[["ac"]]==protein,,drop=FALSE]
    protein.rows <- protein.rows[order(protein.rows[["r1"]],protein.rows[["r2"]]),,drop=FALSE]
    protein.groupnumber <- protein.rows[1,"group"]
    prot.info <- my.protein.infos[[protein]]
  
    reporter.peptides <- peptides(protein.group,protein=protein)
    spectra <- names(spectrumToPeptide(protein.group))[spectrumToPeptide(protein.group)%in%reporter.peptides]
  
    if (all(is.na(protein.rows[,'lratio']))) {
      next
    }
  
    protein.names <- prot.info[["collapsed.gene_name"]][["protein_name"]]
    gene.names <- .tex.combinenames(protein.names,is.single.comparision,cmbn)
  
    mycat(mr(sprintf(" \\hyperref[protein.%s]{\\textbf{%s}}",
                   protein.groupnumber,protein.groupnumber)),
        " & ",mrp(paste0(prot.info[["table.name"]],ifelse(length(gene.names)>0,": ",""),
                 paste(gene.names,collapse=", "))),
      #  " & ",mr(paste(prot.info[["n.reporter"]],prot.info[["n.groupmember"]],sep="/")),
        " & ",mr(print_groupsize(prot.info[["n.reporter"]],prot.info[["n.groupmember"]])),
        " & ",mr(length(reporter.peptides)),
        " & ",mr(length(spectra)),sep=" ")
  
    for (i in seq_len(ncol(cmbn))) {
      if (i > 1) { mycat(paste(rep("&",4),collapse=" ")) }
      if (!is.single.comparision) {
          mycat(" & ",protein.rows[i,'r1'])
          mycat(" & ",protein.rows[i,'r2'])
      }
      mycat(" & ",ifelse(is.na(protein.rows[i,'n.spectra']) | 
                    protein.rows[i,'n.spectra']==0,"",protein.rows[i,'n.spectra']))
      if (is.na(protein.rows[i,"lratio"])) {
        mycat (" & & & ")
      } else {
        mycat(" & ",sprintf("\\textbf{%.2f}",10^protein.rows[i,"lratio"]))
        mycat(" & ",ifelse(protein.rows[i,"is.significant"] == 1,"*",""))
        mycat(" & ",draw.boxplot(protein.rows[i,'lratio'],
                               protein.rows[i,'sd'],bnd))
      }
      mycat(" \\\\")
      if (i < ncol(cmbn)) mycat("*\n")
    }
    if (protein != proteins[length(proteins)] && ncol(cmbn) > 1) {
        mycat(" \\midrule[0.02em] \n\n");
    }
  }
  mycat("\\end{longtable}\n")
  message("done")
}

print_classlabels_tbl <- function(cl,reporterTagNames) {
  if (!is.null(cl)) {
    cl[is.na(cl)] <- ""
    if (!is.null(cl)) {
      if (is.null(names(cl))) {
        cat("\\begin{tabular}{rp{4cm}}\n")
        cat("Channel & Class \\\\ \\hline \\bigstrut[t]\n")
        cl <- paste("\\emph{",cl,"}",sep="")
      } else {
        cat("\\begin{tabular}{rrM}\n")
        cat("\\multicolumn{3}{c}{Class Labels} \\bigstrut \\\\ \n")
        cat("Channel & Class & Description \\\\ \\hline \\bigstrut[t]\n")
        cl <- paste(paste("\\emph{",cl,"} & "),names(cl),sep="")
      }
      for (i in seq_along(reporterTagNames)) {
        cat(reporterTagNames[i]," & ",
            cl[i],
            "\\\\ \n")
      }
      cat("\\end{tabular}\n")
    }
  }
  
}

print_protein_notquant_tbl <- function(file="",
                                    cmbn,protein.group,
                                    quant.tbl,
                                    my.protein.infos) {

  is.single.comparision <- ncol(cmbn) ==1
  mycat <- function(...,append=TRUE,sep="") 
    cat(...,file=file,append=append,sep=sep)


  tt <- table(quant.tbl[is.na(quant.tbl[["lratio"]]),"ac"])
  proteins.notquantified <- names(tt)[tt==ncol(cmbn)]
  mycat("\\begin{longtable}{rXrrr}",
  "  \\#",
  "  & \\textbf{protein}",
  "  & \\rh{group}",
  "  & \\rh{peptides}",
  "  & \\rh{spectra}",
  "\\endhead",sep="\n",append=FALSE)
  
  for (protein in proteins.notquantified) {
    protein.rows <- quant.tbl[quant.tbl[["ac"]]==protein,,drop=FALSE]
    protein.rows <- protein.rows[order(protein.rows[["r1"]],protein.rows[["r2"]]),,drop=FALSE]
    protein.groupnumber <- protein.rows[1,"group"]
    prot.info <- my.protein.infos[[protein]]
  
    reporter.peptides <- peptides(protein.group,protein=protein)
    spectra <- names(spectrumToPeptide(protein.group))[spectrumToPeptide(protein.group)%in%reporter.peptides]
  
    protein.names <- prot.info[["collapsed.gene_name"]][["protein_name"]]
    gene.names <- .tex.combinenames(protein.names,is.single.comparision,cmbn)
  
    mycat(sprintf(" \\hyperref[protein.%s]{%s}",
                   protein.groupnumber,protein.groupnumber),
        " & ",paste0(prot.info[["table.name"]],": ",
                 paste(gene.names,collapse=", ")),
        " & ",print_groupsize(prot.info[["n.reporter"]],prot.info[["n.groupmember"]]),
        " & ",length(reporter.peptides),
        " & ",length(spectra),"\\\\ \n",sep=" ")
  }
  mycat("\\end{longtable}\n")
}

print_protein_grp_tbl <- function(file="",proteins, protein.group) {
  mycat <- function(...,append=TRUE,sep="") 
    cat(...,file=file,append=append,sep=sep)

  mycat("\\begin{longtable}{rcXp{6cm}}\n",append=FALSE)
  prot_i <- 1
  for (reporter.protein.g in proteins) {
    draw.proteingroup.row(sprintf("\\label{protein.%s} %s",prot_i,prot_i),protein.group,reporter.protein.g,file=file)
    mycat("\\bigstrut\n")
    if (prot_i%%100 == 0) message(".",appendLF=FALSE);
    prot_i <- prot_i + 1
  }
  mycat("\\end{longtable}\n")
  
}

print_protein_grp_info <- function() {
cat("
\\section{Protein Group Details} \\label{sec:proteingroups}

Protein groups are created by grouping according to mass spectrometry evidence:

\\begin{tabularx}{\\textwidth}{rX}
1. & Proteins which are detected with the exact same set of peptides are clustered - they are indistinguishable based on available MS data. Very often these are splice variants. In the following text, proteins are clustered proteins (cluster size one if they have specific peptides).  \\\\ 
2. & Proteins, which have unique peptides (\\emph{i.\\,e.}~peptides not present in other detected proteins) are promoted to \\emph{reporters}. \\\\ 
3. & Remaining proteins are grouped to those reporters which have a superset of the peptides they have. These proteins are termed \\emph{group members} \\\\ 
\\end{tabularx}

Peptides are characterized based on their specificity to one or more proteins, and one or more protein groups:
\\begin{center}
\\begin{tabular}{rll}
\\emph{reporter specific} rs &  \\tikz \\node[rs,minimum size=0.2cm] {};  & Number of peptides specific to reporter. \\\\
\\emph{group specific} gs    &  \\tikz \\node[gs,minimum size=0.2cm] {};  & Number of peptides specific to the group (including rs). \\\\
\\emph{unspecific} us        &  \\tikz \\node[us,minimum size=0.2cm] {};  & Number of peptides shared with any protein. \\\\
\\end{tabular}
\\end{center}
")
}
