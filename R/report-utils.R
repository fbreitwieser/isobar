
testPdflatex <- function(pdflatex.cmd) {
  if (system(pdflatex.cmd) != 0) 
    stop("pdflatex does not seem to be installed. ",
         "Install LaTeX using the TeXLive or MiKTex distribution to generate ",
         "PDF reports with isobar. Set 'write.qc.report=FALSE' and 'write.report=FALSE', otherwise.")
}
    

create.reports <- function(properties.file="properties.R",
                           global.properties.file=system.file("report","properties.R",package="isobar"),
                           args=NULL,...,recreate.properties.env=TRUE,recreate.report.env=TRUE) {
  ow <- options("warn")
  if (!exists("properties.env") || recreate.properties.env) {
    properties.env <- load.properties(properties.file,
                                      global.properties.file,
                                      args=args,...)
    assign("properties.env",properties.env,envir=.GlobalEnv)
  }
  options(warn=properties.env[["warning.level"]])

  if (!exists("report.env") || recreate.report.env) {
    report.env <- .GlobalEnv
    initialize.env(report.env,properties.env)
  }

  zip.files <- c(properties.file)

  ## generate XLS report
  if(property('write.xls.report',properties.env)) {
    message("Writing isobar-analysis.xls")
    write.xls.report(properties.env,report.env)
    zip.files <- c(zip.files,"isobar-analysis.xls")
  }
 
  ## generate Latex/Sweave report
  if(property('write.qc.report',properties.env)) {
    message("Weaving isobar-qc report")
    Sweave(system.file("report","isobar-qc.Rnw",package="isobar"))
    if (property('use.name.for.report',properties.env)) {
      qc.name <- .sanitize.sh(sprintf("%s.qc",property('name',properties.env)))
    	file.rename("isobar-qc.tex",sprintf("%s.tex",qc.name))
    } else {
        qc.name <- "isobar-qc"
    }

    zip.files <- c(zip.files,sprintf("%s.tex",qc.name))
    if (properties.env[["compile"]]) 
      zip.files <- .compile.tex(qc.name,zip.files)
  }

  if(property('write.report',properties.env) && properties.env[["report.level"]] != 'peptide') {
    message("Weaving isobar-analysis report")
    name <- switch(properties.env[["report.level"]],
                   protein="isobar-analysis",
                   peptide="isobar-peptide-analysis",
                   stop(properties.env[["report.level"]]," report.level not known",
                        " - choose protein or peptide"))
    Sweave(system.file("report",paste(name,".Rnw",sep=""),package="isobar"))

    if (property('use.name.for.report',properties.env)) {
    	tex.name <- sprintf("%s.tex",name)
      name <- .sanitize.sh(sprintf("%s.quant",property('name',properties.env)))
    	file.rename(tex.name,sprintf("%s.tex",name))
    } else {
    }


    zip.files <- c(zip.files,sprintf("%s.tex",name))
    if (properties.env[["compile"]])
      zip.files <- .compile.tex(name,zip.files)
  }

  if (properties.env[["zip"]]) {
    zip.f <- sprintf("%s.zip",property('name',properties.env))
    zip(zip.f,zip.files)
    message("Created zip archive ",zip.f)
  }

  options(ow) 
  message("\nSUCCESSFULLY CREATED REPORTS\n")
}

  .compile.tex <- function(name,zip.files) {
    r.cmd = shQuote(file.path(R.home("bin"),"R"))
    testPdflatex(paste(r.cmd,"CMD","pdflatex -version"))
    
    dir <- tempdir()
    cat("compiling ",name,".tex ...  1",sep="")
    .call.cmd(sprintf("%s CMD pdflatex -halt-on-error -output-directory=%s %s.tex",r.cmd,dir,name),
              paste(dir,"/",basename(name),".stdout",sep=""))
    cat(" 2")
    .call.cmd(sprintf("%s CMD pdflatex -halt-on-error -output-directory=%s %s.tex",r.cmd,dir,name),
              paste(dir,"/",basename(name),".stdout",sep=""))
    cat(" done!\n\n")
    file.copy(file.path(dir,paste0(name,".pdf")),file.path(getwd(),paste0(name,".pdf")),overwrite=TRUE)
    file.remove(file.path(dir,paste0(name,".pdf")))
    c(zip.files,sprintf("%s.pdf",name))
  }
 
#- load properties

load.properties <- function(properties.file="properties.R",
                            global.properties.file=system.file("report","properties.R",package="isobar"),
                            args=NULL,...) {

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
    properties.env[["properties.file"]] <- properties.file
    properties.env[["properties.file.content"]] <- readLines(properties.file)
  } else {
    message("  No local properties file.")
  }
  
  ## command argument parsing
  tmp.properties.env <- new.env()
  if (length(args) > 0)
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

  ## ... parsing
  dotargs <- list(...)
  if (length(dotargs) > 0) {
    message("parsing function arguments")
    .env.copy(properties.env,list2env(dotargs))
  }

  return(properties.env)
}

#- initialize environment
initialize.env <- function(env,properties.env) {
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

  ib.name <- file.path(.get.property('cachedir',properties.env),"ibspectra.rda")
  if (file.exists(ib.name)) {
    load(ib.name)
    env[["ibspectra"]] <- get("ibspectra")
  } else {
    env[["ibspectra"]] <- .create.or.load.ibspectra(properties.env)
    save(list='ibspectra',envir=env,file=ib.name)
  }
  if ("site.probs" %in% colnames(fData(env[["ibspectra"]])) 
      && ! "pep.siteprobs" %in% colnames(fData(env[["ibspectra"]])))
    fData(env[["ibspectra"]])$pep.siteprobs <- 
      .convertPhosphoRSPepProb(fData(env[["ibspectra"]])[,'peptide'],fData(env[["ibspectra"]])[,'site.probs'],
                               round.to.frac=20)


  env[["noise.model"]] <- .create.or.load.noise.model(env,properties.env)
  env[["ratiodistr"]] <- .create.or.load.ratiodistr(env,properties.env)
  if (identical(properties.env[["report.level"]],"peptide") ) {
    env[["ptm.info"]]  <- .create.or.load.ptm.info(properties.env,proteinGroup(env[["ibspectra"]]))
    if (all(c('use.for.quant','pep.siteprobs') %in% colnames(fData(env[['ibspectra']]))) && !all(fData(env[['ibspectra']])[['use.for.quant']]))
      env[["quant.tbl.notlocalized"]] <- 
        .create.or.load.quant.table(env,properties.env,name="quant.tbl.notlocalized",type="other-sites")
  }
  env[["quant.tbl"]] <- .create.or.load.quant.table(env,properties.env)
  if (!"ac" %in% colnames(env[["quant.tbl"]]) && "protein" %in% colnames(env[["quant.tbl"]]))
    env[["quant.tbl"]][,'ac'] <- env[["quant.tbl"]]$protein

  ## required for TeX
  if (property('write.report',properties.env) && identical(properties.env[["report.level"]],"protein"))
    env[["my.protein.infos"]] <- .create.or.load.my.protein.infos(env,properties.env)

  if (property('write.xls.report',properties.env)) {
    xls.table.name <- paste0("xls.quant.tbl.",properties.env[["xls.report.format"]])
    env[["xls.quant.tbl"]] <- 
      .create.or.load.xls.quant.tbl(env,properties.env,xls.table.name,"quant.tbl")
    if (exists("quant.tbl.notlocalized",envir=env)) {
      xls.ntable.name <- paste0("xls.quant.tbl.notlocalized.",properties.env[["xls.report.format"]])
      env[["xls.quant.tbl.notlocalized"]] <- 
        .create.or.load.xls.quant.tbl(env,properties.env,xls.ntable.name,"quant.tbl.notlocalized")
    }
  }
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
  if (file.exists(file.name) && (!envir[["regen"]] || do.load)) {
    message(sprintf("  loading %s from %s ... ",msg.f, file.name),appendLF=FALSE)
    x <- .load.property(name,file.name)
    message("done")
    return(x)
  } else {
    stop("Cannot get or load property [",name,"] - would expect it in [",file.name,"], but file does not exist.")
  }
}

.create.or.load <- function(name,envir,f,msg.f=name,regenerate=FALSE,do.load=FALSE,class=NULL,error=stop,default.value=NULL,...) {
  x <- if ( !regenerate ) tryCatch(.get.or.load(name,envir,msg.f,class,do.load=do.load),error=function(e) .DOES.NOT.EXIST) else .DOES.NOT.EXIST
  if (identical(x,.DOES.NOT.EXIST)) {
    message(paste("  creating",msg.f,"... "),appendLF=FALSE)
    tryCatch({
      x <- f(...)
      assign(name,x)
      file.name <- sprintf("%s/%s.rda",.get.property('cachedir',envir),name)
      message("done ",appendLF=FALSE)
      save(list=c(name),file=file.name)
      message("& saved to ",file.name)

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

property <- function(x, envir, null.ok=TRUE,class=NULL) {
  if (!.exists.property(x,envir)) 
    stop("property ",x," does not exist!")

  o <- get(x,envir=envir,inherits=FALSE)

  if (!null.ok && is.null(o))
    stop("property ",x," may not be NULL, but it is!")

  if (!is.null(class) && !is(o,class))
    stop("property ",x," does not have class ",class,"!")

  return(o)
}

.exists.property <- function(x,envir,null.ok=TRUE) 
  exists(x,envir=envir,inherits=FALSE) &&
    (null.ok || !is.null(.get.property(x,envir)))

.check.property <- function(x,envir,inherits=FALSE,msg="",print=TRUE,valid=NULL,check.file=FALSE,...) {
  def <- grep(x,envir[["properties.file.content"]],value=TRUE)
  if (length(def) > 0)
    def <- paste("\n  Corresponding line in ",envir[["properties.file"]],":\n\t",
                 paste(def,collapse="\n\t"),"\n\n",sep="")

  p <- .get.property(x,envir=envir)
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
    isnt.valid <- !p %in% valid
    if (isnt.valid) stop(" property '",x,"' must be one of \n\t",paste(valid,collapse="\n\t"),",\n\n",
                         "  but it is ",p )
  }
  if (print) {
    message("  property '",x,"' defined",ifelse(is.character(p),paste0(": ",paste(p,collapse="; ")),"."))
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
  readIBSpectra.args[["type"]]=get.property('type')
  readIBSpectra.args[["id.file"]]=get.property('ibspectra')
  readIBSpectra.args[["fragment.precision"]]=get.property('fragment.precision')
  readIBSpectra.args[["fragment.outlier.prob"]]=get.property('fragment.outlier.prob')
  readIBSpectra.args[["proteinGroupTemplate"]]=.get.or.load('protein.group.template',properties.env,"ProteinGroup",null.ok=TRUE)

  if (is.data.frame(get.property('ibspectra')) || all(file.exists(get.property('ibspectra')))) {
    if (is.data.frame(get.property('ibspectra'))) {
      readIBSpectra.args[["id.file"]] <- get.property('ibspectra')
      ibspectra <- do.call(readIBSpectra,readIBSpectra.args)
    } else if (grepl(".csv",get.property('ibspectra'))) {
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

    readIBSpectra.args[["id.file"]]=get.property('identifications')
    readIBSpectra.args[["peaklist.file"]]=get.property('peaklist')

    ibspectra <- do.call(readIBSpectra,readIBSpectra.args)
    if (grepl(".csv",get.property('ibspectra'))) {
        write.table(as.data.frame(ibspectra),sep="\t",row.names=F,file=get.property('ibspectra'),quote=FALSE)
    } else if (grepl(".rda",get.property("ibspectra"))) {
        save(ibspectra,file=get.property("ibspectra"),compress=TRUE)
    }  else {
        stop("weird naming of ibspectra: ",get.property("ibspectra"),". use XXX.csv or XXX.rda")
    }
  }

  if (!is.null(get.property("isotope.impurities"))) {
    isotopeImpurities(ibspectra) <- get.property("isotope.impurities")
  }

  if (property('correct.isotope.impurities',properties.env))
    ibspectra <- correctIsotopeImpurities(ibspectra)
  if (property('normalize',properties.env)) {
    ibspectra <-
      normalize(ibspectra,
                use.protein=property('normalize.use.protein',properties.env),
                exclude.protein=property('normalize.exclude.protein',properties.env),
                peptide.specificity=property('peptide.specificity',properties.env),
                f=property('normalize.function',properties.env),
                channels=property('normalize.channels',properties.env),
                na.rm=property('normalize.na.rm',properties.env),
                normalize.factors=property('normalize.factors',properties.env))
  }

  class.labels <- LETTERS[1:length(reporterTagNames(ibspectra))]
  if (.exists.property('class.labels',properties.env,null.ok=FALSE))
    class.labels <- get.property('class.labels')
  if (!any(table(class.labels)>1) && property('summarize',properties.env)) {
    stop("When summarize=TRUE, the must be more then one channel per class")
  }
  if (!is.character(class.labels)) {
    stop("Please provide class.labels of class character!")
  }
  classLabels(ibspectra) <- class.labels

  if (!is.null(property('protein.info.f',properties.env)))
    tryCatch({
    proteinInfo(proteinGroup(ibspectra)) <- 
      .create.or.load("protein.info",envir=properties.env,
                      f=property('protein.info.f',properties.env),
                      x=proteinGroup(ibspectra),
                      do.load=TRUE, msg.f="protein.info",
                      error=warning,default.value=proteinInfo(proteinGroup(ibspectra)))
    },error=function(e) stop("Error creating proteinInfo using function defined in [protein.info.f]: ",e,"Set protein.info.f to NULL if you want to create a report anyway."))

  if ("site.probs" %in% colnames(fData(ibspectra)) 
      && ! "pep.siteprobs" %in% colnames(fData(ibspectra)))
     fData(ibspectra)$pep.siteprobs <- 
        .convertPhosphoRSPepProb(fData(ibspectra)[,'peptide'],fData(ibspectra)[,'site.probs'],
                                 round.to.frac=20)
 
  return(ibspectra)
}

.create.or.load.noise.model <- function(env,properties.env) {
  noise.model.channels <- .get.property("noise.model.channels",properties.env)
  if (is.null(noise.model.channels)) 
    noise.model.channels <- reporterTagNames(env[["ibspectra"]])[!is.na(classLabels(env[["ibspectra"]]))]
  is.one.to.one <- isTRUE(.get.property("noise.model.is.technicalreplicates",properties.env))

  noise.model <- .get.property("noise.model",properties.env)
  if (!is(noise.model,"NoiseModel")) {
    noise.model.f <- .get.property("noise.model",properties.env)
    if (!file.exists(noise.model.f)) {
      message("estimating noise model as ",ifelse(is.one.to.one,"","non"),
              " one-to-one from channels ",paste0(noise.model.channels,collapse=", ")," ...")
      noise.model <- new("ExponentialNoiseModel",env[["ibspectra"]],one.to.one=is.one.to.one,
                         reporterTagNames=noise.model.channels,
                         min.spectra=property('noise.model.minspectra',properties.env))
      save(noise.model,file=noise.model.f,compress=TRUE)
    } else {
      message(sprintf("  loading noise model from %s ...",noise.model.f))
      noise.model <- .load.property("noise.model",noise.model.f)
    }
  }
  return(noise.model)
}

.create.or.load.ptm.info <- function(properties.env,protein.group) {
  if (is.null(property('ptm.info.f',properties.env))) 
    properties.env[["ptm.info.f"]] <- getPtmInfoFromNextprot

  return(.create.or.load("ptm.info",envir=properties.env,
                         f=property('ptm.info.f',properties.env),
                         protein.group=protein.group))
}

.create.or.load.ratiodistr <- function(env,properties.env) {
  
  return(.create.or.load("ratiodistr",envir=properties.env,class="Distribution",
                         msg.f="biological variability ratio distribution",f=function(){

    cl <- classLabels(env[["ibspectra"]])
    if (!is.null(property('ratiodistr.class.labels',properties.env)))
      cl <- property('ratiodistr.class.labels',properties.env)

    ratios.for.distr.fitting <- .create.or.load.ratiodistr.ratios(env,properties.env,cl)

    if (all(is.na(ratios.for.distr.fitting[,'lratio']))) 
      stop("All ratios for distr fitting are NA - are the correct class labels used?")

    if (!is.function(property('ratiodistr.fitting.f',properties.env)))
      stop("ratiodistr.fitting.f must be set to a function [e.g. fitCauchy or fitTd]")

    ratiodistr <- property('ratiodistr.fitting.f',properties.env)(ratios.for.distr.fitting[,'lratio'])
    ratiodistr <- .round.distr(ratiodistr,digits=5)
    attr(ratiodistr,"combn.method") <- attr(ratios.for.distr.fitting,"combn.method")
    attr(ratiodistr,"cl") <- attr(ratios.for.distr.fitting,"cl")
    attr(ratiodistr,"tagNames") <- reporterTagNames(env[["ibspectra"]])
    ratiodistr
  }))
}

.round.distr <- function(distr,digits) {
  for (s in slotNames(param(distr))) 
    if (s != "name" && is.numeric(slot(param(distr),s))) 
      slot(distr@param,s) <- round(slot(param(distr),s),digits)
  distr
}

.create.or.load.ratiodistr.ratios <- function (env,properties.env,cl) {
  .create.or.load("ratios.for.distr.fitting",envir=properties.env,class="data.frame",
                  msg.f="ratios for biological variability distribution fitting",f=function() {

    message("  Using class labels ",paste(cl,collapse=", "))
    if (property('summarize',properties.env) || any(table(cl)>=2)) {
      method <- "intraclass"
    } else {
      message(" WARNING: ratiodistr will be computed based on global ratios")
      method <- "global"
    }
    if (any(table(properties.env[["class.labels"]])>2))
      do.summarize <- properties.env[["summarize"]]
    else
      do.summarize <- FALSE

    if (identical(properties.env[["report.level"]],"peptide"))
      ratios.for.distr.fitting <- peptideRatios(env[["ibspectra"]],noise.model=env[["noise.model"]],do.warn=FALSE,
                                  cl=cl,combn.method=method,symmetry=isTRUE(properties.env[["ratiodistr.symmetry"]]),summarize=do.summarize)
    else
      ratios.for.distr.fitting <- proteinRatios(env[["ibspectra"]],noise.model=env[["noise.model"]],do.warn=FALSE,
                                      cl=cl,combn.method=method,symmetry=isTRUE(properties.env[["ratiodistr.symmetry"]]),summarize=do.summarize)

    if (all(is.nan(ratios.for.distr.fitting[["lratio"]])))
      stop("Cannot compute protein ratio distribution - no ratios available.\n",
           "Probably due to missing reporter intensities.")

    attr(ratios.for.distr.fitting,"combn.method") <- method
    attr(ratios.for.distr.fitting,"cl") <- cl

    ratios.for.distr.fitting
  })
}


.set <- function(x,name,list) {
  if (name %in% names(list)) {
    stop(name," already assigned in list!")
  }
  x
}

.create.or.load.quant.table <- function(env,properties.env,name="quant.tbl",type='confident-sites') {
  protein.group <- proteinGroup(env[["ibspectra"]])
  protein.info <- proteinInfo(protein.group)
  isoforms <- protein.group@isoformToGeneProduct
  
  .create.or.load(name,envir=properties.env,class="data.frame",
                  msg.f=paste("table of ratios of",properties.env[["report.level"]],"as",name),f=function() {
    if (!is.null(property('ratios.opts',properties.env)$summarize)) {
        message("WARNING: ratio.opts$summarize will be overwritten,",
                " define it outside of ratio.opts!")
        warning("ratio.opts$summarize will be overwritten,",
                " define it outside of ratio.opts!")
    }
    ratios.opts <- property('ratios.opts',properties.env)
    set.ratioopts <- function(x,name=names(x)) {
      if (is.null(intersect(name,names(ratios.opts))))
          stop("property '",intersect(name,names(ratios.opts)),
               "' already assigned in ratios.opts list - check your properties file!")

      if (is.list(x))
        ratios.opts <<- c(ratios.opts,x)
      else
        ratios.opts[[name]] <<- x
    }

    set.ratioopts(name="ibspectra",env[["ibspectra"]])
    set.ratioopts(name="noise.model",env[["noise.model"]])
    set.ratioopts(name="ratiodistr",env[["ratiodistr"]])
    
    if(identical(properties.env[["report.level"]],"peptide")) {
      if (!"use.for.quant" %in% colnames(fData(env[["ibspectra"]])))
        fData(env[["ibspectra"]])[["use.for.quant"]] <- TRUE
      if (type=='confident-sites')
        pep.n.modif <- unique(apply(fData(env[["ibspectra"]])[fData(env[["ibspectra"]])[["use.for.quant"]],
                                    c("peptide","modif")],2,cbind))
      else {
        pep.n.modif <- unique(apply(fData(env[["ibspectra"]])[!fData(env[["ibspectra"]])[["use.for.quant"]],
                                    c("peptide","modif","pep.siteprobs")],2,cbind))
        set.ratioopts(list(use.for.quant.only=FALSE))
      }


      set.ratioopts(list(peptide=pep.n.modif,
                         proteins=NULL))

    } else if (identical(properties.env[["report.level"]],"protein")) {
      set.ratioopts(list(peptide=NULL,
                         proteins=reporterProteins(proteinGroup(env[["ibspectra"]])),
                         quant.w.grouppeptides=property('quant.w.grouppeptides',properties.env)))
    } else {
      stop("don't known level ",properties.env[["report.level"]])
    }

    if (is.null(property('cmbn',properties.env)) & 
	!is.null(property('vs.class',properties.env)))
      properties.env[["cmbn"]] <- combn.matrix(reporterTagNames(env[["ibspectra"]]),
					   "versus.class",
					   property('class.labels',properties.env),
					   vs=property('vs.class',properties.env))

    set.ratioopts(list(combn.method=property('combn.method',properties.env),
                       cl=classLabels(env[["ibspectra"]]),
                       summarize=property('summarize',properties.env),
                       cmbn=property('cmbn',properties.env),
                       use.na=property('use.na',properties.env),
                       zscore.threshold=properties.env[["zscore.threshold"]],
                       do.warn=FALSE))

    if (!is.null(property('correct.peptide.ratios.with',properties.env))) {
      protein.quant.tbl <- .get.or.load("correct.peptide.ratios.with",properties.env)
      correct.protein.group <- if ( !is.null(property('correct.peptide.ratios.with_protein.group',properties.env)) ) {
        .get.or.load("correct.peptide.ratios.with_protein.group",properties.env)
      } else {
        warning( '"correct.peptide.ratios.with_protein.group" property not found, ',
                 'using the existing protein group' )
        .get.or.load("protein.group.template",properties.env)
      }
      ratios.opts[["before.summarize.f"]] <- function(...)
        correct.peptide.ratios(..., protein.quant.tbl=protein.quant.tbl,
                               protein.group.combined = correct.protein.group,
                               correlation = property('peptide.protein.correlation',properties.env))
    }

    quant.tbl <- do.call("proteinRatios",ratios.opts)

    quant.tbl[,"sd"] <- sqrt(quant.tbl[,"variance"])
    
#    quant.tbl$sign.string <- "not significant"
#    quant.tbl$sign.string[quant.tbl$is.significant] <- "is significant"
    
#    if (length(property('preselected',properties.env)) > 0) {
#      preselected <- unique(ip[sub("-.*","",names(ip)) %in% property('preselected',properties.env)])
#      quant.tbl$is.preselected <- quant.tbl$protein %in% preselected

#      quant.tbl$sign.string <- "not significant"
#      quant.tbl$sign.string[is.sign & !quant.tbl$is.preselected] <- "is significant [ours]"
#      quant.tbl$sign.string[is.sign & quant.tbl$is.preselected] <- "is significant [both]"
#      quant.tbl$sign.string[!is.sign & quant.tbl$is.preselected] <- "is significant [theirs]"
#    }

    if (identical(properties.env[["report.level"]],"protein")) {
      quant.tbl[,"gene_names"] <- sapply(quant.tbl[,"ac"], function(x) {
        if (length(protein.info) == 0) return("")
        allreporter <- indistinguishableProteins(protein.group,protein.g=x)
        acs <- unique(isoforms[allreporter,"proteinac.wo.splicevariant"])
        paste(sort(unique(protein.info[protein.info[["accession"]] %in% acs,"gene_name"])),
              collapse=", ")
      })
      sort.genenames <- quant.tbl[,"gene_names"]
      sort.genenames[sort.genenames==""] <- quant.tbl[sort.genenames=="","ac"]
      quant.tbl <- quant.tbl[order(sort.genenames,quant.tbl[,"r1"],quant.tbl[,"r2"]),]
      quant.tbl[,"group"] <- as.numeric(factor(quant.tbl[,"ac"],levels=unique(quant.tbl[,"ac"])))
    }
  
    if (all(is.na(quant.tbl[["lratio"]])))
      stop("All ratios are NA")

    return(quant.tbl)
  })
}

.create.or.load.my.protein.infos <- function(env,properties.env) {
  .create.or.load("my.protein.infos",envir=properties.env,
                  f=function() {

    protein.groupnames <-unique(env[["quant.tbl"]][,"ac"])
    ## if (is.null(protein.info)) { stop("protein info is null!")}                
    my.protein.infos <- llply(protein.groupnames, .do.create.protein.info, protein.group=proteinGroup(env[["ibspectra"]]), 
                              .parallel=isTRUE(getOption('isobar.parallel')))
    #my.protein.infos <- lapply(protein.groupnames, .do.create.protein.info, protein.group=protein.group)
    names(my.protein.infos) <- protein.groupnames
    return(my.protein.infos)
  })
}
  
.do.create.protein.info <- function(x,protein.group) {
    protein.group.table <- proteinGroupTable(protein.group)
      allgroupmember <- indistinguishableProteins(protein.group, protein.g =
                                                  protein.group.table[protein.group.table[,'reporter.protein'] %in% x,"protein.g"])
     
      reporter.protein.info <- my.protein.info(protein.group,x)
      collapsed.gene_name <- human.protein.names(reporter.protein.info)
      peptides <- peptides(protein.group,protein=x,do.warn=FALSE)
      peptides.gs <- peptides(protein.group,protein=x,
                              specificity=GROUPSPECIFIC,do.warn=FALSE)
      peptides.rs <- peptides(protein.group,protein=x,
                              specificity=REPORTERSPECIFIC,do.warn=FALSE)
      n.spectra <- length(names(spectrumToPeptide(protein.group))[spectrumToPeptide(protein.group)%in%peptides])
      
      tbl.protein.name <- sort(collapsed.gene_name[,'protein_name'])[1];
      if (length(unique(collapsed.gene_name[,'protein_name'])) > 1)
        tbl.protein.name <- paste(tbl.protein.name,", ...",sep="")
      
      list(n.reporter = nrow(reporter.protein.info),
           n.groupmember = length(allgroupmember),
           reporter.protein.info = reporter.protein.info,
           n.peptides=length(peptides),
           n.spectra=n.spectra,
           collapsed.gene_name = collapsed.gene_name,
           table.name = ifelse(
             length(collapsed.gene_name[,'ac_link'])>3,
             paste(paste(collapsed.gene_name[1:3,'ac_link'],collapse=", "),
                   ", \\dots",sep=""),
             paste(collapsed.gene_name[,'ac_link'],collapse=", ")),
           section.name = sanitize(paste(unique(collapsed.gene_name[,'name_nolink']),
             collapse=", ")),
           table.protein.name = tbl.protein.name,
           gene.name = paste(sort(unique(reporter.protein.info[,'gene_name'])),collapse=", ")
           )
    }

## copys objects from env into parentenv.
## Those objects MUST exist in parentenv before.
.env.copy <- function (parentenv,env) {
  printval <- function(x) if(is.character(x)) capture.output(dput(x)) else toString(x)

  for (object in ls(envir=env)) {
    if (exists(object,envir=parentenv)) {
      ## assign object to parent env
      val <- get(object,envir=env)
      if (is.list(val)) {
        val.s <- paste0("list(",paste(names(val),sapply(val,printval),sep=": ",collapse=";"),")")
      } else {
        val.s <- tryCatch(printval(val),error=function(e) capture.output(print(val)))
      }
      message(sprintf("%30s = %s",object,paste(val.s,collapse=paste0("\n",paste(rep(" ",35),collapse="")))))
      assign(object, value=get(object,envir=env),envir=parentenv)
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



