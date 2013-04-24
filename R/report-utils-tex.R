#################################################
## Helper Functions for LaTeX Report Generation

print_longtablehdr <- function(level,coldef,draw.channels,ncol.p,draw.signcol) {
  #cat("\n\n\\renewcommand{\\arraystretch}{0.75}\n")
  cat("\\begin{longtable}{",coldef,"}",'\n',sep="")
  cat("\t \\# ",
      paste(" & \\textbf{",level,"}"),
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
  col <- "black!60"

  ratio.smaller.bnd <- lratio - sd < -bnd
  ratio.bigger.bnd  <- lratio + sd > bnd

  return(sprintf("\\boxplot{%.2f}{%.2f}{%s!%s}{%s}{%s}{%.0f}{%s}\n",
                  lratio,
                  sd,
                  ifelse(lratio > 0,"green","red"),
                  min(floor(abs(lratio)/bnd*100),100),
                  ifelse(ratio.smaller.bnd, "\\divcol", "black!1"),
                  ifelse(ratio.bigger.bnd,  "\\divcol", "black!1"),
                  sd,
                  col))
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
draw.proteingroup.row <- function(name,protein.group,reporter.protein.g) {
  gmp <- groupMemberPeptides(protein.group,reporter.protein.g,TRUE)
  pgt <- proteinGroupTable(protein.group)

  n.multirow <- length(unique(protein.ac(protein.group,pgt[pgt$reporter.protein==reporter.protein.g,"protein.g"])))
  cat("\\multirow{",n.multirow ,"}{*}{",name,"} & \n",sep="")
  n.row <- 1
  show.pos <- FALSE
  for (protein.i in seq_len(ncol(gmp$group.member.peptides))) {
    x = colnames(gmp$group.member.peptides)[protein.i]

    my.protein.info <- my.protein.info(protein.group,x)
    for (ac in unique(my.protein.info$accession)) {
      if (n.row > 1) cat(" & ") 
      sel <- my.protein.info$accession == ac
      var.string <- number.ranges(my.protein.info$splicevariant[sel])
      #cat(protein.i,"&")
      if (show.pos) cat(protein.i,"&")
      cat(paste(sprintf("\\uniprotlink{%s}",sanitize(ac,dash=FALSE)),
          ifelse(is.na(var.string),"",var.string)[1],
          sanitize(unique(my.protein.info$gene_name[sel]))[1],
          sanitize(unique(my.protein.info$protein_name[sel]))[1],
          sep=" & "))

      if (n.row==1) {
        # PEPTIDE GROUPING
        cat(" & \\multirow{",n.multirow ,"}{*}{%\n",sep="")
        tikz.proteingroup(protein.group,reporter.protein.g,show.pos,show.header=FALSE)
        cat("} \\\\ \n")
      } else {
        cat(" & \\\\ \n")
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

tikz.proteingroup <- function(protein.group,reporter.protein.g,show.pos,show.header=TRUE) {

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
  if (show.header)
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
   
     .print.proteinrow(
        -protein.i,paste(which(peptide.idx.sel),peptide.styles[peptide.idx.sel],sep="/")
     )
   
     if (n.peptides < 15) {
       connect.nodes(-protein.i,which((reporter.sp.sel | group.sp.sel | unspecific.sel | quant.sel)&gmp$group.member.peptides[,protein.i]))
     }
  }
  cat(sprintf("\\node{%s/%s/%s};\n",protein.peptides.df[protein.i,"rs"],
      protein.peptides.df[protein.i,"gs"],
      protein.peptides.df[protein.i,"us"]))


  if (show.header) {
  cat("  \\matrix[ampersand replacement=\\&,matrix anchor=us node.east,anchor=base] at (0,0) {\n")
      cat("    \\node{}; \\& \\node {rs};  \\& \\node {gs}; \\& \\node (us node){us};  \\\\\n")
  for (protein.i in seq_len(n.groupmember)) {
     if (show.pos) cat(sprintf("    \\node[draw,gray]{%s};",protein.i))
     cat(sprintf("    \\& \\node{%s};\\& \\node{%s};\\& \\node{%s}; \\\\\n",
                 protein.peptides.df[protein.i,"rs"],
                 protein.peptides.df[protein.i,"gs"],
                 protein.peptides.df[protein.i,"us"]))
  }
  cat("  };\n")
  }
  cat("\\end{tikzpicture}\n")

}

.print.proteinrow <- function(protein.idx,peptides.sel) {
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
