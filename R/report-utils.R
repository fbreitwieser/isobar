
create.meta.reports <- function(report.type="protein",properties.file="meta-properties.R",args=NULL) {
  source(system.file("report","meta-functions.R",package="isobar"))
  if (!exists("properties.env")) {
    properties.env <- load.properties(properties.file,
                                      system.file("report","meta-properties.R",package="isobar"),
                                      args=args)
  }

  protein.group <- .get.or.load("protein.group",properties.env,"protein group object","ProteinGroup")
  if (!is.null(properties.env$protein.info.f)) 
    proteinInfo(protein.group) <- 
      .create.or.load("protein.info",envir=properties.env,
                      f=properties.env$protein.info.f,
                      x=protein.group,
                      error=warning,default.value=proteinInfo(protein.group))
  
  if (is.null(proteinInfo(protein.group)) || length(proteinInfo(protein.group))==0)
    stop("No protein information available.")

  #if (properties.env$calculate.dnsaf)
  #  dnsaf <- .create.or.load("dnsaf",envir=properties.env,f=calculate.dNSAF,protein.group=protein.group)

  message("Merging tables ...")
  ac.vars <- switch(report.type,
                    protein = "ac",
                    peptide = c("peptide","modif"),
                    stop("report type ",report.type," unknown"))

  merge.cols <- c("class1","class2")
  merged.table <- get.merged.table(properties.env$samples,
                                   cols=c(ac.vars,merge.cols,"lratio","variance"),
                                   merge.by=c(ac.vars,merge.cols),format="wide")
  merged.table$comp <- paste(merged.table[,merge.cols[2]],sep="/",merged.table[,merge.cols[1]])

  tbl.wide <- reshape(merged.table,idvar=ac.vars,timevar="comp",direction="wide",drop=merge.cols)
  sel.zscore.both.f <- function(zscores) (zscores[,1] > 2.5 & zscores[,2] > 2.5) | (zscores[,1] < -2.5 & zscores[,2] < -2.5)
  sel.zscore.one.f <- function(zscores) (abs(zscores[,1]) > 2.5 & is.na(zscores[,2])) | (abs(zscores[,2]) > 2.5 & is.na(zscores[,1]))

  sel.zscore.one  <- sel.zscore.one.f(tbl.wide[,.grep_columns(tbl.wide,"zscore.*.RapaTNF/15mTNF")])
  sel.zscore.both <- sel.zscore.both.f(tbl.wide[,.grep_columns(tbl.wide,"zscore.*.RapaTNF/15mTNF")])

  tbl.wide$selected <- sprintf("not selected [n=%s]",sum(!sel.zscore.one & !sel.zscore.both,na.rm=T))
  tbl.wide$selected[sel.zscore.one] <- sprintf("selected once [n=%s]",sum(sel.zscore.one,na.rm=T))
  tbl.wide$selected[sel.zscore.both] <- sprintf("selected twice [n=%s]",sum(sel.zscore.both,na.rm=T))
  tbl.wide$size = 2
  tbl.wide$size[sel.zscore.one] <- 2.5
  tbl.wide$size[sel.zscore.both] <- 3

  pdf("/media/work/analysis/tk/scatter_ratio.pdf")
  ggplot(tbl.wide,aes(x=`lratio.15mTNF/CTRL`,y=`lratio.RapaTNF/CTRL`)) + 
    geom_point(aes(color=selected,size=size),alpha=0.8) +
#    scale_size_discrete() +
    xlab("log10 ratio (15m TNF/CTRL) [mean of tech reps]") +
    ylab("log10 ratio (TNF+Rapamycin/CTRL) [mean of tech reps]") +
    opts(title="RapaTNF/15mTNF significantly regulated proteins")
  dev.off()

  pdf("/media/work/analysis/tk/distr_zscores.pdf")
  print(ggplot(tbl.wide) + 
    geom_histogram(fill="#0000A0",binwidth=0.3,data=tbl.wide,aes(x=`zscore.P4041-1.RapaTNF/15mTNF`),alpha=0.5,color="black") + 
    geom_histogram(color="black",fill="#A00000",binwidth=0.3,aes(x=`zscore.P4041-2.RapaTNF/15mTNF`),alpha=0.5) +
    geom_vline(xintercept=c(-2.5,2.5),linetype="dashed",color="red") +
    geom_text(aes(x,y,label=label),data=data.frame(x=c(-2.5,2.5),y=-5,label=c("-2.5 MAD","2.5 MAD")),size=4,color="red",hjust=0,vjust=0) +
    #geom_line(aes(x,y),data=data.frame(x=range(pretty(c(tbl.wide[,"zscore.P4041-1.RapaTNF/15mTNF"],tbl.wide[,"zscore.P4041-2.RapaTNF/15mTNF"]),na.rm=TRUE)),y=50)) +
    xlab("log10 (TNF+Rapamycin/TNF only) z-score") + ylab("Frequency") +
    opts(title="Distribution of robust z-scores of RapaTNF/15mTNF",
         subtitle="(Median absolute deviations (MADs) away from the median ratio)",
         plot.title=theme_text(size = 14)))
  dev.off()

  table(tbl.wide$selected)
  #not selected  selected once selected twice 
  #       2912             67              5 

  t(sapply(.grep_columns(tbl.wide,"lratio.P",value=TRUE,logical=FALSE),function(i) isobar:::.sum.bool(!is.na(tbl.wide[,i]))))
  tbl.pg <- merge(tbl.wide,pg.df)
  sel <- apply(!is.na(tbl.pg[,.grep_columns(tbl.wide,"lratio.P")]),1,any)
  sel.1 <- !is.na(tbl.pg[,"lratio.P4041-1.RapaTNF/15mTNF"])
  sel.2 <- !is.na(tbl.pg[,"lratio.P4041-2.RapaTNF/15mTNF"])
  sel.selected <- tbl.pg$selected != "not selected"
  data.frame(
             protein=c(total=length(unique(tbl.pg[,"proteins"])),
                       quant.any=length(unique(tbl.pg[sel.1|sel.2,"proteins"])),
                       quant.1=length(unique(tbl.pg[sel.1,"proteins"])),
                       quant.2=length(unique(tbl.pg[sel.2,"proteins"])),
                       quant.both=length(unique(tbl.pg[sel.1&sel.2,"proteins"])),
                       selected=length(unique(tbl.pg[sel.selected,"proteins"]))),
             peptide=c(total=length(unique(tbl.pg[,"peptide"])),
                       quant.any=length(unique(tbl.pg[sel.1|sel.2,"peptide"])),
                       quant.1=length(unique(tbl.pg[sel.1,"peptide"])),
                       quant.2=length(unique(tbl.pg[sel.2,"peptide"])),
                       quant.both=length(unique(tbl.pg[sel.1&sel.2,"peptide"])),
                       selected=length(unique(tbl.pg[sel.selected,"peptide"]))),
             phospopep=c(total=nrow(tbl.wide),
                         quant.any=sum(sel.1|sel.2),
                         quant.1=sum(sel.1),
                         quant.2=sum(sel.2),
                         quant.both=sum(sel.1&sel.2),
                         selected=sum(sel.selected)))







  merge.cols <- c("class1","class2")
  merged.table <- get.merged.table(properties.env$samples,
                                   cols=c(ac.vars,merge.cols,"lratio","variance"),
                                   merge.by=c(ac.vars,merge.cols))
  merged.table$comp <- paste(merge.cols[2],sep="/",merge.cols[1])


  #merged.table <- subset(merged.table,r1==merged.table$r1[1])
  tbl.wide <- reshape(merged.table,idvar=ac.vars,timevar="comp",direction="wide",drop=merge.cols)
 # rownames(tbl.wide) <- tbl.wide$ac
  #all.names <- do.call(rbind,lapply(tbl.wide[,"ac"],get.names,protein.group=protein.group))
  #tbl.wide$dNSAF <- dnsaf[as.character(tbl.wide$ac)]

  if (report.type=="peptide") {
    pg.df <- isobar:::.proteinGroupAsConciseDataFrame(protein.group)  
    rownames(pg.df) <- do.call(paste,pg.df[,ac.vars,drop=FALSE])
  }

  p <- ggplot(sdf, aes(y,x))
  p <- p + geom <- tile(aes(fill=height), colour="white")
  p <- p + scale <- fill <- gradientn(colours = c("dark red", "white", "dark green" ), breaks=breaks, labels=format(breaks))
  p <- p + opts(title='Day 1') 



  ratio.matrix <- as.matrix(tbl.wide[,grep("lratio",colnames(tbl.wide))])
  variance.matrix <- as.matrix(tbl.wide[,grep("var",colnames(tbl.wide))])
  rownames(ratio.matrix)  <- do.call(paste,tbl.wide[,ac.vars,drop=FALSE])
  rownames(variance.matrix)  <- do.call(paste,tbl.wide[,ac.vars,drop=FALSE])
  sel <- !apply(is.na(ratio.matrix),1,any)
  ratio.matrix <- ratio.matrix[sel,]
  variance.matrix <- variance.matrix[sel,]

  m.median <- apply(ratio.matrix,2,median)
  normalized.ratio.matrix <- ratio.matrix-matrix(m.median,nrow=nrow(ratio.matrix),ncol=ncol(ratio.matrix),byrow=T)

  plot.heatmaps(ratio.matrix,properties.env$name)
  plot.pairs(properties.env$name)

}

create.reports <- function(properties.file="properties.R",args=NULL,
                           report.type="protein",compile=FALSE,zip=FALSE) {
  ow <- options("warn")
  options(warn=1)
  if (!exists("properties.env")) {
    properties.env <- load.properties(properties.file,
                                      system.file("report","properties.R",package="isobar"),
                                      args=args)
  }
  
  if (!exists("report.env")) {
    report.env <- .GlobalEnv
    initialize.env(report.env,report.type,properties.env)
  }

  zip.files <- c(properties.file)

  ## generate XLS report
  if(properties.env$write.xls.report) {
    message("Writing isobar-analysis.xls")
    write.xls.report(report.type,properties.env,report.env)
    zip.files <- c(zip.files,"isobar-analysis.xls")
  }
  
 
  ## generate Latex/Sweave report
  if(properties.env$write.qc.report) {
    message("Weaving isobar-qc report")
    Sweave(system.file("report","isobar-qc.Rnw",package="isobar"))
    if (properties.env$use.name.for.report) {
	qc.name <- sprintf("%s.qc",properties.env$name)
    	file.rename("isobar-qc.tex",sprintf("%s.tex",qc.name))
    } else {
        qc.name <- "isobar-qc"
    }

    zip.files <- c(zip.files,sprintf("%s.tex",qc.name))
    if (compile) 
      zip.files <- .compile.tex(qc.name,zip.files)
  }

  if(properties.env$write.report) {
    message("Weaving isobar-analysis report")
    name <- switch(report.type,
                   protein="isobar-analysis",
                   peptide="isobar-peptide-analysis",
                   stop(report.type," report type not known",
                        " - choose protein or peptide"))
    Sweave(system.file("report",paste(name,".Rnw",sep=""),package="isobar"))

    if (properties.env$use.name.for.report) {
    	tex.name <- sprintf("%s.tex",name)
	name <- sprintf("%s.quant",properties.env$name)
    	file.rename(tex.name,sprintf("%s.tex",name))
    } else {
    }


    zip.files <- c(zip.files,sprintf("%s.tex",name))
    if (compile)
      zip.files <- .compile.tex(name,zip.files)
  }

  if (zip) {
    zip.f <- sprintf("%s.zip",properties.env$name)
    zip(zip.f,zip.files)
    message("Created zip archive ",zip.f)
  }

  options(ow) 
  message("\nSUCCESSFULLY CREATED REPORTS\n")
}

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
      #protein.id.df <- as(get('ibspectra',report.env),"data.frame.concise")
      protein.id.df <- .IBSpectraAsConciseDataFrame(get('ibspectra',report.env))

      ## make columns w/ multiple groups gray
      sel.1group  <- protein.id.df$n.groups == 1
      #sel.1ac  <- protein.id.df$n.acs == 1
      #sel.1variant  <- protein.id.df$n.variants == 1
      #protein.id.df[!sel.1ac & sel.1group,1] <- paste("#color silver#",protein.id.df[!sel.1ac & sel.1group,1],sep="")
      #protein.id.df[!sel.1group | !protein.id.df$use.for.quant,1] <- paste("#color=gray#",protein.id.df[!sel.1group,1],sep="")

    } else {
      #protein.id.df <- as(get('ibspectra',report.env),"data.frame.concise")
      protein.id.df <- .IBSpectraAsConciseDataFrame(get('ibspectra',report.env))
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
                                  xlsx=system.file("pl","tab2xlsx.pl",package="isobar"),
                                  xls=system.file("pl","tab2xls.pl",package="isobar"),
                                  stop("spreadsheet.format property must be either 'xlsx' or 'xls'."))

    perl.cl <- paste(tab2spreadsheet.cmd," ",
                     ifelse(properties.env$use.name.for.report,sprintf("%s.quant",properties.env$name),"isobar-analysis"),
                     ".",properties.env$spreadsheet.format,
                     " ':autofilter,freeze_col=3,name=Quantifications:",protein.quant.f,"'",
                     ifelse(identical(report.type,"peptide") && !is.null(modificationSites),
                            paste(" ':autofilter,freeze_col=3,name=Modification Sites:",modifsites.f,"'",sep=""),""),
                     " ':autofilter,freeze_col=3,name=Identifications:",protein.id.f,"'",
                     " ':name=Analysis Properties:",analysis.properties.f,"'",
                     " ':name=Log:",log.f,"'",sep="")
    
    ## generate Excel report (using Spreadsheet::WriteExcel)
    message(perl.cl)
    system(perl.cl)
}

  .compile.tex <- function(name,zip.files) {
    .call.cmd <- function(cmd,stdout.to=NULL) 
      if (is.null(stdout.to)) {
        if (system(cmd) != 0) stop("\nError executing [",cmd,"]")
      } else {
        if (system(paste(cmd,">",stdout.to)) != 0) 
          stop("\nError executing [",cmd,"]: \n\n ...\n",
               paste(tail(readLines(stdout.to),n=10),collapse="\n"))
      }
    dir <- tempdir()
    cat("compiling ",name,".tex ...  1",sep="")
    .call.cmd(sprintf("R CMD pdflatex -halt-on-error -output-directory=%s %s.tex",dir,name),
              paste(dir,"/",basename(name),".stdout",sep=""))
    cat(" 2")
    .call.cmd(sprintf("R CMD pdflatex -halt-on-error -output-directory=%s %s.tex",dir,name),
              paste(dir,"/",basename(name),".stdout",sep=""))
    cat(" done!\n\n")
    .call.cmd(sprintf("mv %s/%s.pdf .",dir,name))
    c(zip.files,sprintf("%s.pdf",name))
  }
 
#- load properties

load.properties <- function(properties.file="properties.R",
                            global.properties.file=system.file("report","properties.R",package="isobar"),
                            args=NULL) {

  properties.env <- new.env()
  tmp.properties.env <- new.env()

  message("Loading global properties file ",global.properties.file," ...")
  tryCatch(sys.source(global.properties.file,properties.env),
           error=function(e) stop("\nCould not read properties file:\n", e) )

  message("Looking for local properties file ",properties.file," ...")
  if (!is.null(properties.file) && file.exists(properties.file)) {
    message("  Loading local properties file ...")
    ## load properties file
    tryCatch(sys.source(properties.file,tmp.properties.env),
             error=function(e) stop("\nCould not read properties file:\n", e) )
    .env.copy(properties.env,tmp.properties.env)
    properties.env$properties.file <- properties.file
    properties.env$properties.file.content <- readLines(properties.file)
  } else {
    message("  No local properties file.")
  }
  
  ## command argument parsing
  tmp.properties.env <- new.env()
  message("parsing command line arguments ...")
  for (arg in args) {
    if (grepl("^--",arg)) {
      arg.n.val <- strsplit(substring(arg,3),"=")[[1]]
      if (length(arg.n.val) == 1)
        tmp.properties.env[[arg.n.val]] <- TRUE
      else if (length(arg.n.val) == 2) 
        tmp.properties.env[[arg.n.val[1]]] <- switch(arg.n.val[2],'TRUE'=TRUE,'FALSE'=FALSE,arg.n.val[2])
      else
        stop("Could not parse command line argument ",arg)
    } else {
        stop("Could not parse command line argument ",arg)
    }
    .env.copy(properties.env,tmp.properties.env)
  }

  return(properties.env)
}

#- initialize environment
initialize.env <- function(env,report.type="protein",properties.env) {
  ## get property and exists property convenience functions
  get.property <- function(name) .get.property(name,properties.env)

  if (file.exists(get.property('cachedir'))) {
    if (!file.info(get.property('cachedir'))$isdir)
      stop("Cannot write into cachedir [",get.property('cachedir'),"]",
           " - it is a file and no directory")
  } else {
    ret <- dir.create(get.property('cachedir'))
    if (!ret) stop("Error creating cachedir [",get.property('cachedir'),"]")
  }

  env$ibspectra <- .create.or.load.ibspectra(properties.env)
  env$noise.model <- .create.or.load.noise.model(env,properties.env)
  env$ratiodistr <- .create.or.load.ratiodistr(env,properties.env,level=report.type)
  if (identical(report.type,"peptide") )
    env$ptm.info  <- .create.or.load.ptm.info(env,properties.env)
  env$quant.tbl <- .create.or.load.quant.table(env,properties.env,level=report.type)
  if (!"ac" %in% colnames(env$quant.tbl) && "protein" %in% colnames(env$quant.tbl))
    env$quant.tbl$ac <- env$quant.tbl$protein

  ## required for TeX
  if (identical(report.type,"protein"))
    env$my.protein.infos <- .create.or.load.my.protein.infos(env,properties.env)
  env$xls.quant.tbl <- .create.or.load.xls.quant.tbl(report.type,env,properties.env)
  #if (report.type == "protein") {
  #  env$my.protein.infos <- .create.or.load.my.protein.infos(env,properties.env)
  #  env$xls.quant.tbl <- .create.or.load.xls.quant.tbl(env,properties.env)
  #} else if (report.type == "peptide") {
    ## compute peptide ratios
  #  env$xls.peptide.tbl <- .create.or.load.xls.peptide.tbl(env,properties.env)
  #} else {
  #  stop("report type [",report.type,"] not known - choose protein or peptide")
  #}
}

#- property loading helper functions
.DOES.NOT.EXIST = "NOT AVAILABLE"

.get.or.load <- function(name,envir,msg.f=name,class=NULL,null.ok=FALSE,do.load=FALSE) {
  if (.exists.property(name,envir,null.ok=null.ok)) {
    o <- .get.property(name,envir)
    if (is.null(o) && null.ok) return(NULL)
    if (!is.null(class) && is(o,class)) {
      return(o)
    } else if (is(o,"character")) {
      file.name <- o
    } else {
      stop("property ",name," is neither 'character' nor of a specified class")
    }
  } else {
    file.name <- sprintf("%s/%s.rda",.get.property('cachedir',envir),name)
  }
  if (file.exists(file.name) & (!envir$regen || do.load)) {
    message(sprintf("loading %s from %s ...",msg.f, file.name))
    x <- .load.property(name,file.name)
    return(x)
  } else {
    stop("Cannot get or load property ",name)
  }
}

.create.or.load <- function(name,envir,f,msg.f=name,do.load=FALSE,class=NULL,error=stop,default.value=NULL,...) {
  x <- tryCatch(.get.or.load(name,envir,msg.f,class,do.load=do.load),error=function(e) .DOES.NOT.EXIST)
  if (identical(x,.DOES.NOT.EXIST)) {
    message(paste("creating",msg.f,"..."))
    tryCatch({
      x <- f(...)
      assign(name,x)
      file.name <- sprintf("%s/%s.rda",.get.property('cachedir',envir),name)
      save(list=c(name),file=file.name)
    },error=error)

    if (!exists(name,inherits=FALSE)) {
      x <- default.value
      assign(name,x)
    }
  }
  if (!is.null(class) && !is(x,class))
    stop("property [",name,"] should be of class [",class,"] but is of class [",class(x),"]")

  return(x)
}
 
.get.property <- function(x,envir) get(x,envir=envir,inherits=FALSE)

.exists.property <- function(x,envir,null.ok=TRUE) 
  exists(x,envir=envir,inherits=FALSE) &&
    (null.ok || !is.null(.get.property(x,envir)))

.check.property <- function(x,envir,inherits=FALSE,msg="",print=TRUE,valid=NULL,check.file=FALSE,...) {
  def <- grep(x,envir$properties.file.content,value=TRUE)
  if (length(def) > 0)
    def <- paste("\n  Corresponding line in ",envir$properties.file,":\n\t",
                 paste(def,collapse="\n\t"),"\n\n",sep="")

  does.not.exist <- !exists(x,envir=envir,inherits=inherits,...)
  is.null.p <- does.not.exist || is.null(.get.property(x,envir=envir))
  length.0 <- does.not.exist || length(.get.property(x,envir=envir))==0
  if (does.not.exist || is.null.p || length.0) 
    stop ("  property '",x,"' not defined or did not evaluate in the properties file, but it is required",
          ifelse(is.null(valid),"",paste(" to be one of \n\t",paste(valid,collapse="\n\t"),
          sep="")),def)
  
  if (check.file && !file.exists(.get.property(x,envir=envir)))
    stop("  property '",x,"' should be assigned to a file, but it is not:",def)

  if (!is.null(valid)) {
    p <- .get.property(x,envir=envir)
    isnt.valid <- !p %in% valid
    if (isnt.valid) stop(" property '",x,"' must be one of \n\t",paste(valid,collapse="\n\t"),",\n\n",
                         "  but it is ",p )
  }
  if (print) {
    message(paste("  property '",x,"' defined: ",
                  paste(.get.property(x,envir=envir),collapse=","),sep=""))
  }
}

.load.property <- function(name,file) {
  tmp.env <- new.env()
  load(file,envir=tmp.env)
  if (exists(name,envir=tmp.env))
    return(get(name,envir=tmp.env))
  else if (length(ls(envir=tmp.env)) == 1)
    return(get(ls(envir=tmp.env)[1],envir=tmp.env))
  else
    stop("Could not load property ",name," from file ",file)
}

.create.or.load.ibspectra <- function(properties.env) {
  get.property <- function(name) .get.property(name,properties.env)
  check.property <- function(name,...) .check.property(name,properties.env,...)
  ## check that required properties are defined
  check.property('ibspectra')
  check.property('type',valid=IBSpectraTypes())

  readIBSpectra.args <- get.property('readIBSpectra.args')
  if (is.null(readIBSpectra.args)) readIBSpectra.args <- list()
  for (name in names(readIBSpectra.args)) {
      arg <- readIBSpectra.args[[name]]
      if (!is.null(arg)) {
          if (!is.null(names(arg)))
            arg <- paste(names(arg),arg,collapse="=")
          if (is.function(arg)) { message("    ",name,": function")
          } else {
          message("    ",name,": ",
                  paste(arg,collapse=ifelse(length(arg)>2,"\n\t",", ")))
          }
      }
  }
  readIBSpectra.args$type=get.property('type')
  readIBSpectra.args$id.file=get.property('ibspectra')
  readIBSpectra.args$fragment.precision=get.property('fragment.precision')
  readIBSpectra.args$fragment.outlier.prob=get.property('fragment.outlier.prob')
  readIBSpectra.args$proteinGroupTemplate=.get.or.load('protein.group.template',properties.env,"ProteinGroup",null.ok=TRUE)

  if (all(file.exists(get.property('ibspectra')))) {
    if (grepl(".csv",get.property('ibspectra'))) {
        message("ibspectra ends on .csv:\n",
                sprintf('ibspectra <- readIBSpectra("%s",%s) ...',
                        get.property('type'),
                        ifelse(length(get.property('ibspectra'))==1,
                               paste0('"',get.property('ibspectra'),'"'),
                               paste0('c("',get.property('ibspectra'),'")',collapse='","')
                               )))
        ibspectra <- do.call(readIBSpectra,readIBSpectra.args)
    } else if (grepl(".rda",get.property("ibspectra"))) {
        message("ibspectra ends on .rda: ","loading")
        ibspectra <- .load.property("ibspectra",get.property("ibspectra"))
    } else {
        stop("weird naming of ibspectra: ",get.property("ibspectra"),". use XXX.csv or XXX.rda")
    }

  } else {
    # create ibspectra file
    message(sprintf('\n  file %s does not exist, creating ibspectra',get.property('ibspectra')))
    check.property('peaklist',print=FALSE,check.file=TRUE)
    check.property('identifications',print=FALSE,check.file=TRUE)
    message(sprintf('    id.files = %s%s\n    peaklist.files = %s%s\n',
          ifelse(length(get.property('identifications'))>1,"\n      ",""),
          paste(get.property('identifications'),collapse="\n      "),
          ifelse(length(get.property('peaklist'))>1,"\n      ",""),
          paste(get.property('peaklist'),collapse="\n      ")))

    readIBSpectra.args$id.file=get.property('identifications')
    readIBSpectra.args$peaklist.file=get.property('peaklist')

    ibspectra <- do.call(readIBSpectra,readIBSpectra.args)
    if (grepl(".csv",get.property('ibspectra'))) {
        write.table(as.data.frame(ibspectra),sep="\t",row.names=F,file=get.property('ibspectra'))
    } else if (grepl(".rda",get.property("ibspectra"))) {
        save(ibspectra,file=get.property("ibspectra"),compress=TRUE)
    }  else {
        stop("weird naming of ibspectra: ",get.property("ibspectra"),". use XXX.csv or XXX.rda")
    }
  }

  if (!is.null(get.property("isotope.impurities"))) {
    isotopeImpurities(ibspectra) <- get.property("isotope.impurities")
  }

 
  if (properties.env$correct.isotope.impurities)
    ibspectra <- correctIsotopeImpurities(ibspectra)
  if (properties.env$normalize) {
    ibspectra <-
      normalize(ibspectra,
                use.protein=properties.env$normalize.use.protein,
                exclude.protein=properties.env$normalize.exclude.protein,
                f=properties.env$normalize.function,
                channels=properties.env$normalize.channels,
                na.rm=properties.env$normalize.na.rm)
  }

  class.labels <- as.character(c(1,rep(0,length(sampleNames(ibspectra))-1)))
  if (.exists.property('class.labels',properties.env,null.ok=FALSE))
    class.labels <- get.property('class.labels')
  if (!any(table(class.labels)>1) && properties.env$summarize) {
    stop("When summarize=TRUE, the must be more then one channel per class")
  }
  if (!is.character(class.labels)) {
    stop("Please provide class.labels of class character!")
  }
  classLabels(ibspectra) <- class.labels

  if (!is.null(properties.env$protein.info.f))
    proteinInfo(proteinGroup(ibspectra)) <- 
      .create.or.load("protein.info",envir=properties.env,
                      f=properties.env$protein.info.f,
                      x=proteinGroup(ibspectra),
                      do.load=TRUE, msg.f="protein.info",
                      error=warning,default.value=proteinInfo(proteinGroup(ibspectra)))
 
  return(ibspectra)
}

.create.or.load.noise.model <- function(env,properties.env) {
  noise.model.channels <- .get.property("noise.model.channels",properties.env)

  noise.model <- .get.property("noise.model",properties.env)
  if (!is(noise.model,"NoiseModel")) {
    noise.model.f <- .get.property("noise.model",properties.env)
    if (!file.exists(noise.model.f)) {
      message("estimating noise model as non one-to-one ...")
      noise.model <- new("ExponentialNoiseModel",env,one.to.one=F,
                         reporterTagNames=noise.model.channels,
                         min.spectra=properties.env$noise.model.minspectra)
      save(noise.model,file=noise.model.f,compress=TRUE)
    } else {
      message(sprintf("loading noise model from %s ...",noise.model.f))
      noise.model <- .load.property("noise.model",noise.model.f)
    }
  }
  return(noise.model)
}

.create.or.load.ptm.info <- function(env,properties.env) {
  return(.create.or.load("ptm.info",envir=properties.env,
                         f=getPtmInfoFromNextprot,
                         protein.group=proteinGroup(env$ibspectra)))
}

.create.or.load.ratiodistr <- function(env,properties.env,level) {
  
  return(.create.or.load("ratiodistr",envir=properties.env,class="Distribution",
                         msg.f="protein ratio distribution",f=function(){
    if (properties.env$summarize) {
      method <- "intraclass"
   } else {
      message(" WARNING: ratiodistr will be computed based on global ratios!")
      method <- "global"
    }

    if (identical(level,"peptide"))
      all.ratios <- peptideRatios(env$ibspectra,noise.model=env$noise.model,do.warn=FALSE,
                                  peptide=peptides(proteinGroup(env$ibspectra)),
                                  cl=classLabels(env$ibspectra),combn.method=method,symmetry=TRUE)
    else
      all.ratios <- proteinRatios(env$ibspectra,noise.model=env$noise.model,do.warn=FALSE,
                                      proteins=reporterProteins(proteinGroup(env$ibspectra)),peptide=NULL,
                                      cl=classLabels(env$ibspectra),combn.method=method,symmetry=TRUE)

    if (all(is.nan(all.ratios$lratio)))
      stop("Cannot compute protein ratio distribution - no ratios available.\n",
           "Probably due to missing reporter intensities.")
    fitCauchy(all.ratios[,'lratio'])
  }))
}

.set <- function(x,name,list) {
  if (name %in% names(list)) {
    stop(name," already assigned in list!")
  }
  x
}

.create.or.load.quant.table <- function(env,properties.env,level) {
  protein.group <- proteinGroup(env$ibspectra)
  protein.info <- proteinInfo(protein.group)
  isoforms <- protein.group@isoformToGeneProduct
  
  .create.or.load("quant.tbl",envir=properties.env,class="data.frame",
                  msg.f=paste("table of ratios of",level),f=function() {
    if (!is.null(properties.env$ratios.opts$summarize)) {
        message("WARNING: ratio.opts$summarize will be overwritten,",
                " define it outside of ratio.opts!")
        warning("ratio.opts$summarize will be overwritten,",
                " define it outside of ratio.opts!")
    }
    ratios.opts <- properties.env$ratios.opts
    set.ratioopts <- function(x,name=names(x)) {
      if (is.null(intersect(name,names(ratios.opts))))
          stop("property '",intersect(name,names(ratios.opts)),
               "' already assigned in ratios.opts list - check your properties file!")

      if (is.list(x))
        ratios.opts <<- c(ratios.opts,x)
      else
        ratios.opts[[name]] <<- x
    }

    set.ratioopts(name="ibspectra",env$ibspectra)
    set.ratioopts(name="noise.model",env$noise.model)
    set.ratioopts(name="ratiodistr",env$ratiodistr)
    
    if(identical(level,"peptide")){
      set.ratioopts(list(
                         peptide=unique(apply(fData(env$ibspectra)[,c("peptide","modif")],2,cbind)),
                         proteins=NULL))
    } else if (identical(level,"protein")) {
      set.ratioopts(list(peptide=NULL,
                         proteins=reporterProteins(proteinGroup(env$ibspectra)),
                         quant.w.grouppeptides=properties.env$quant.w.grouppeptides))
    } else {
      stop("don't known level ",level)
    }
    if (is.null(properties.env$combn) & !is.null(properties.env$vs.class))
      properties.env$combn <- combn.matrix(reporterTagNames(env$ibspectra),
					   "versus.class",
					   properties.env$class.labels,
					   vs=properties.env$vs.class)

    set.ratioopts(list(combn.method=properties.env$combn.method,
                       cl=classLabels(env$ibspectra),
                       summarize=properties.env$summarize,
                       combn=properties.env$combn,
                       use.na=properties.env$use.na,
                       do.warn=FALSE))

    quant.tbl <- do.call("proteinRatios",ratios.opts)
    quant.tbl[,"sd"] <- sqrt(quant.tbl[,"variance"])
    
#    quant.tbl$sign.string <- "not significant"
#    quant.tbl$sign.string[quant.tbl$is.significant] <- "is significant"
    
#    if (length(properties.env$preselected) > 0) {
#      preselected <- unique(ip[sub("-.*","",names(ip)) %in% properties.env$preselected])
#      quant.tbl$is.preselected <- quant.tbl$protein %in% preselected

#      quant.tbl$sign.string <- "not significant"
#      quant.tbl$sign.string[is.sign & !quant.tbl$is.preselected] <- "is significant [ours]"
#      quant.tbl$sign.string[is.sign & quant.tbl$is.preselected] <- "is significant [both]"
#      quant.tbl$sign.string[!is.sign & quant.tbl$is.preselected] <- "is significant [theirs]"
#    }

    if (identical(level,"protein")) {
      quant.tbl[,"gene_names"] <- sapply(quant.tbl[,"ac"], function(x) {
        if (length(protein.info) == 0) return("")
        allreporter <- indistinguishableProteins(protein.group,protein.g=x)
        acs <- unique(isoforms[allreporter,"proteinac.wo.splicevariant"])
        paste(sort(unique(protein.info[protein.info$accession %in% acs,"gene_name"])),
              collapse=", ")
      })
      sort.genenames <- quant.tbl[,"gene_names"]
      sort.genenames[sort.genenames==""] <- quant.tbl[sort.genenames=="","ac"]
      quant.tbl <- quant.tbl[order(sort.genenames,quant.tbl[,"r1"],quant.tbl[,"r2"]),]
      quant.tbl[,"group"] <- as.numeric(factor(quant.tbl[,"ac"],levels=unique(quant.tbl[,"ac"])))
    }
    return(quant.tbl)
  })
}

.create.or.load.xls.quant.tbl <- function(report.type,env,properties.env) {
  .create.or.load("xls.quant.tbl",envir=properties.env,
                  msg.f="protein table for Excel export",f=function() {
    message("XLS report format: ",properties.env$xls.report.format)
                    
    protein.group <- proteinGroup(env$ibspectra)
    indist.proteins <- indistinguishableProteins(protein.group)

    if (!is.null(properties.env$compare.to.quant))
      compare.to.quant <- properties.env$compare.to.quant
    else
      compare.to.quant <- NULL

    if (isTRUE(properties.env$xls.report.format=="wide")) {
      xls.quant.tbl.tmp  <- ratiosReshapeWide(env$quant.tbl,vs.class=properties.env$vs.class,sep="###",cmbn=properties.env$combn)
      if (!is.null(compare.to.quant))
        compare.to.quant <- lapply(compare.to.quant,ratiosReshapeWide,vs.class=properties.env$vs.class,sep="###",cmbn=properties.env$combn)
    } else {
      xls.quant.tbl.tmp <- env$quant.tbl
    }

    round.digits <- 4;
    if (identical(report.type,"protein")) {
      xls.quant.tbl.tmp$i  <- seq_len(nrow(xls.quant.tbl.tmp))
      xls.quant.tbl <- data.frame(i=xls.quant.tbl.tmp$i,
                                  group=xls.quant.tbl.tmp[,"group"],
                 AC=.protein.acc(xls.quant.tbl.tmp[,"ac"],ip=indist.proteins),
                 ID=proteinInfo(protein.group,xls.quant.tbl.tmp[,"ac"],do.warn=FALSE),
                 n=sapply(xls.quant.tbl.tmp[,"ac"],function(p) {length(names(indist.proteins)[indist.proteins == p])}),
                 Description=proteinInfo(protein.group,xls.quant.tbl.tmp[,"ac"],select="protein_name",do.warn=FALSE),
                 Gene=proteinInfo(protein.group,xls.quant.tbl.tmp[,"ac"],select="gene_name",do.warn=FALSE),
                 "@comment=Number of specific peptides@Peptide Count"= peptide.count(protein.group,xls.quant.tbl.tmp$ac,specificity=c(GROUPSPECIFIC,REPORTERSPECIFIC),do.warn=FALSE),
                 "@comment=Number of specific spectra@Spectral Count"= spectra.count(protein.group,xls.quant.tbl.tmp$ac,specificity=c(GROUPSPECIFIC,REPORTERSPECIFIC),do.warn=FALSE),
                 "@comment=Coverage ratio of the sequence compared to the full protein sequence@Sequence Coverage"=round(sequence.coverage(protein.group,xls.quant.tbl.tmp$ac,do.warn=FALSE),round.digits),
                 check.names=FALSE)

    } else {
      ## PEPTIDE REPORT
      pnp  <- subset(as.data.frame(peptideNProtein(protein.group),stringsAsFactors=FALSE),
                     protein.g %in% reporterProteins(protein.group))

      my.ptm <- "PTM"
      if ("PHOS" %in% properties.env$ptm) my.ptm="Phosphorylation"
      if ("METH" %in% properties.env$ptm) my.ptm="Methylation"

      xls.quant.tbl.tmp$ac  <- NULL
      t <- table(pnp$peptide)
      pnp <- pnp[pnp$peptide %in% names(t)[t==1],]
      colnames(pnp)  <- c("peptide","ac")

      xls.quant.tbl.tmp <- merge(pnp,xls.quant.tbl.tmp,by="peptide")
      xls.quant.tbl.tmp$i  <- seq_len(nrow(xls.quant.tbl.tmp))
      xls.quant.tbl.tmp$Spectra <- apply(xls.quant.tbl.tmp,1,function(x) nrow(subset(fData(env$ibspectra),peptide==x['peptide'] & modif==x['modif'])))
      pg.df <- .proteinGroupAsConciseDataFrame(protein.group,modif.pos=properties.env$ptm,ptm.info=env$ptm.info)
      
      xls.quant.tbl <- merge(pg.df,xls.quant.tbl.tmp[,c("peptide","modif","i","Spectra")],by=c("peptide","modif"),all.y=TRUE)
      xls.quant.tbl$peptide <- .convertPeptideModif(xls.quant.tbl[,"peptide"],xls.quant.tbl[,"modif"])
      colnames(xls.quant.tbl)[colnames(xls.quant.tbl)=="peptide"] <- "Sequence"
      colnames(xls.quant.tbl)[colnames(xls.quant.tbl)=="proteins"] <- "ACs"
      colnames(xls.quant.tbl)[colnames(xls.quant.tbl)=="modif.pos"] <- sprintf("@comment=Absolute modification position in protein. Modifications in the same protein are separated by '&', in different proteins by ';'. Stars denote positions which are annotated as phosphorylated in NextProt.@%s Position",my.ptm)
      xls.quant.tbl$modif <- NULL
      xls.quant.tbl$pos <- NULL
      xls.quant.tbl$n.groups <- NULL
      xls.quant.tbl <- xls.quant.tbl[order(xls.quant.tbl$i),]
    }
    if (!is.null(compare.to.quant))
      for (ii in seq_along(compare.to.quant)) 
        xls.quant.tbl.tmp=merge(xls.quant.tbl.tmp,compare.to.quant[[ii]],by="ac",
                                all.x=TRUE,suffixes=c("",paste("###",names(compare.to.quant)[ii])))
    xls.quant.tbl.tmp <- xls.quant.tbl.tmp[order(xls.quant.tbl.tmp$i),]
#        xls.quant.tbl.tmp=merge(xls.quant.tbl.tmp,compare.to.quant[[ii]],by="ac",all.x=TRUE,suffixes=c("",".proteome"))


    get.cols <- function(df,cc,cc.new=NULL,f=NULL,...) {
      data.cc <- df[,grep(cc,colnames(df)),drop=FALSE]
      if (!is.null(f)) data.cc <- f(data.cc,...)
      if (!is.null(cc.new)) colnames(data.cc) <- gsub(cc,cc.new,colnames(data.cc))
      data.cc
    }
    combine.n.append.xls.tbl <- function(cc1,cc2,cc.new,f) {
      A <- get.cols(xls.quant.tbl.tmp,cc1)
      B <- get.cols(xls.quant.tbl.tmp,cc2)
      data.cc <- sapply(1:ncol(A), function(i) f(A[,i],B[,i]) )
      data.cc <- round(data.cc,round.digits)
      colnames(data.cc) <- gsub(cc1,cc.new,colnames(A))
      data.cc
      #xls.quant.tbl <<- cbind(xls.quant.tbl,data.cc)
    }
    append.xls.tbl <- function(...)
      get.cols(xls.quant.tbl.tmp,...)
      #xls.quant.tbl <<- cbind(xls.quant.tbl,get.cols(xls.quant.tbl.tmp,...))

    round.n.append.xls.tbl <- function(...,digits=round.digits)
      round(get.cols(xls.quant.tbl.tmp,...),digits=digits)
      #xls.quant.tbl <<- cbind(xls.quant.tbl,round(get.cols(xls.quant.tbl.tmp,...),digits=digits))


    if (properties.env$sum.intensities) {
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
      
      xls.quant.tbl <- cbind(xls.quant.tbl,protein.intensities(ibspectra,protein.tbl$protein))
    } else {
      if (isTRUE(properties.env$xls.report.format=="long")) {
       xls.quant.tbl <-cbind(xls.quant.tbl,
                             "Channels"=paste(xls.quant.tbl.tmp$r2,"/",xls.quant.tbl.tmp$r1))
      }

      if ("zscore" %in% properties.env$xls.report.columns) {
        ## TODO: zscore is calculated across all classes - 
        ##       it is probably more appropriate to calculate it individual for each class
        xls.quant.tbl.tmp$zscore <- calc.zscore(xls.quant.tbl.tmp$lratio)
      }
    
      # TODO: Add z score?
      for (cc in properties.env$xls.report.columns) {
        res <- switch(cc,
              log10.ratio =    round.n.append.xls.tbl("lratio","log10.ratio"),
              log2.ratio =     round.n.append.xls.tbl("lratio","log2.ratio",f=function(x) x/log10(2)),
              log10.variance = append.xls.tbl("variance","log10.var"),
              log2.variance =  append.xls.tbl("variance","log2.var",f=function(x) (sqrt(x)/log10(2)^2)),
              is.significant = append.xls.tbl("is.significant"),
              n.na1 =          append.xls.tbl("n.na1"),
              n.na2 =          append.xls.tbl("n.na2"),
              p.value.ratio =  append.xls.tbl("p.value.rat"),
              p.value.ratio.adjusted =  append.xls.tbl("p.value.rat.adjusted"),
              p.value.sample = append.xls.tbl("p.value.sample"),
              z.score =        round.n.append.xls.tbl("zscore"),
              ratio =          round.n.append.xls.tbl("lratio","ratio",f=function(x) 10^x),
              CI95.lower =     combine.n.append.xls.tbl("lratio","variance","CI95.lower",f=function(x,y) 10^qnorm(0.025,x,sqrt(y))),
              CI95.upper =     combine.n.append.xls.tbl("lratio","variance","CI95.upper",f=function(x,y) 10^qnorm(0.975,x,sqrt(y))),
              ratio.minus.sd = combine.n.append.xls.tbl("lratio","variance","ratio.minus.sd",f=function(x,y) 10^(x-sqrt(y))),
              ratio.plus.sd = combine.n.append.xls.tbl("lratio","variance","ratio.plus.sd",f=function(x,y) 10^(x+sqrt(y))),
              warning("ignoring unknown column ",cc," in Excel report"))
        
        if (is(res,"data.frame") || is(res,"matrix"))
          xls.quant.tbl <- cbind(xls.quant.tbl,res)
        else
          warning("ignore res",cc)

      }
    }
    if (properties.env$summarize) {
      append.xls.tbl("n.pos")
      append.xls.tbl("n.neg")
    }

    if (length(properties.env$preselected) > 0) {
      ## xls.quant.tbl <- cbind(xls.quant.tbl,"is.preselected"=xls.quant.tbl$is.preselected)
    }
    xls.quant.tbl$i <- NULL

    if (identical(report.type,"protein")) {
      return(xls.quant.tbl[order(xls.quant.tbl[,"group"]),])
    } else {
      return(xls.quant.tbl[order(xls.quant.tbl$ID,xls.quant.tbl$Sequence),])
    }
  })
}

.protein.acc <- function(prots,protein.info=NULL,ip=NULL) {
  if (is.null(ip)) {
    proteins <- list(prots)
  } else {
    proteins <- lapply(prots,function(p) {names(ip)[ip == p]})
  }

  sapply(proteins,function(prots) {
         ## consider ACs with -[0-9]*$ as splice variants (ACs w/ more than one dash are not considered)
         pos.splice <- grepl("^[^-]*-[0-9]*$",prots)
         df <- data.frame(protein=prots,accession=prots,splice=0,stringsAsFactors=FALSE)

         if (any(pos.splice))
           df[pos.splice,c("accession","splice")] <- 
             do.call(rbind,strsplit(prots[pos.splice],"-"))

         res <- 
           ddply(df,"accession",function(y) {
                 if(sum(y$splice>0) <= 1)
                   return(data.frame(protein=unique(y$protein)))
                 else 
                   return(data.frame(protein=sprintf("%s-[%s]",unique(y$accession),
                                                     paste(sort(y[y$splice>0,'splice']),collapse=","))))
                                 })
         return(paste(res$protein,collapse=", "))
  })
}

.create.or.load.xls.peptide.tbl <- function(env,properties.env) {
  .create.or.load("xls.peptide.tbl",envir=properties.env,
                  msg.f="peptide table for Excel export",f=function() {
                    
    protein.group <- proteinGroup(env$ibspectra)
    indist.proteins <- indistinguishableProteins(protein.group)
    df.pg <- as(protein.group,"data.frame.concise")    

    xls.peptide.tbl <- merge(df.pg,quant.tbl,by="peptide")
    xls.peptide.tbl$Channels <- paste(xls.peptide.tbl$r2,"/",xls.peptide.tbl$r1)
    message(paste(colnames(xls.peptide.tbl),collapse=":"))
    return(xls.peptide.tbl)
#    return(xls.peptide.tbl[,order(xls.peptide.tbl[,"proteins"],
#                                  xls.peptide.tbl[,"peptide"],
#                                  xls.peptide.tbl[,"modif"],
#                                  xls.peptide.tbl[,"Channels"]
#                              )])
  })
}

.create.or.load.my.protein.infos <- function(env,properties.env) {
  .create.or.load("my.protein.infos",envir=properties.env,
                  f=function() {

    protein.group <- proteinGroup(env$ibspectra)
    protein.group.table <- proteinGroupTable(protein.group)
    protein.groupnames <-unique(env$quant.tbl[,"ac"])
    ## if (is.null(protein.info)) { stop("protein info is null!")}                
    my.protein.infos <- lapply(protein.groupnames, function(x) {
      allgroupmember <- indistinguishableProteins(protein.group, protein.g =
                                                  protein.group.table$protein.g[protein.group.table$reporter.protein%in%x])
     
      reporter.protein.info <- my.protein.info(protein.group,x)
      collapsed.gene_name <- human.protein.names(reporter.protein.info)
      peptides <- peptides(protein.group,protein=x,do.warn=FALSE)
      peptides.gs <- peptides(protein.group,protein=x,
                              specificity=GROUPSPECIFIC,do.warn=FALSE)
      peptides.rs <- peptides(protein.group,protein=x,
                              specificity=REPORTERSPECIFIC,do.warn=FALSE)
      n.spectra <- length(names(spectrumToPeptide(protein.group))[spectrumToPeptide(protein.group)%in%peptides])
      
      tbl.protein.name <- sort(collapsed.gene_name$protein_name)[1];
      if (length(unique(collapsed.gene_name$protein_name)) > 1)
        tbl.protein.name <- paste(tbl.protein.name,", ...",sep="")
      
      list(n.reporter = nrow(reporter.protein.info),
           n.groupmember = length(allgroupmember),
           reporter.protein.info = reporter.protein.info,
           n.peptides=length(peptides),
           n.spectra=n.spectra,
           collapsed.gene_name = collapsed.gene_name,
           table.name = ifelse(
             length(collapsed.gene_name$ac_link)>3,
             paste(paste(collapsed.gene_name$ac_link[1:3],collapse=", "),
                   ", \\dots",sep=""),
             paste(collapsed.gene_name$ac_link,collapse=", ")),
           section.name = sanitize(paste(collapsed.gene_name$name_nolink,
             collapse=", ")),
           table.protein.name = tbl.protein.name,
           gene.name = paste(sort(unique(reporter.protein.info$gene_name)),collapse=", ")
           )
    })
    names(my.protein.infos) <- protein.groupnames
    return(my.protein.infos)
  })
}
  

## copys objects from env into parentenv.
## Those objects MUST exist in parentenv before.
.env.copy <- function (parentenv,env) {
  for (object in ls(envir=env)) {
    if (exists(object,envir=parentenv)) {
      ## assign object to parent env
      assign(object, value=get(object,envir=env),
             envir=parentenv)
      ## remove(object,envir=env)
    } else {
      stop("Illegal property ",object,"!\n",
           "Call script with --help for usage details.\n\n",
           "Available properties:\n",
           "\t",paste(ls(envir=parentenv),collapse="\n\t"),"\n")
                      
      ## more verbose list:
      ##cat(paste(unlist(lapply(objects(),function(o) {sprintf("%s [%s]",o,class(get(o)))})),collapse="\n\t"))
    }
  }
}


##########################
## latex helper functions

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
    human.protein.name <- human.protein.names(my.protein.info)
  } 
  cat("\\end{tabular}\n")

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
   
     .print.proteinrow(
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
