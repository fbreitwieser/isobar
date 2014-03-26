### =========================================================================
### NoiseModel objects.
### -------------------------------------------------------------------------
###
### Class definition.

# A virtual class representing the interface to noise models.
# Use ExponentialNoiseModel or InverseNoiseModel
setClass("NoiseModel",
    representation(
        na.region="numeric",
        low.intensity="numeric",
        f="function",
        parameter="numeric",
        "VIRTUAL"),
    contains = "VersionedBiobase",
    prototype = prototype(
        new("VersionedBiobase",versions=c(NoiseModel="1.0.0"))
    )
)

#b*e^(cx)
setClass("ExponentialNoANoiseModel",
    contains = "NoiseModel",
    prototype = prototype(
        new("VersionedBiobase",versions=c(classVersion("NoiseModel"),ExponentialNoANoiseModel="1.0.0")),
        f = function(data,parameter) { parameter[1] * exp (- data * parameter[2]) },
        parameter = c(10,1),
        na.region = 0,
        low.intensity = 0
    )
)

#a + b*e^(cx)
setClass("ExponentialNoiseModel",
    contains = "NoiseModel",
    prototype = prototype(
        new("VersionedBiobase",versions=c(classVersion("NoiseModel"),ExponentialNoiseModel="1.0.0")),
        f = function(data,parameter) { parameter[1] + parameter[2] * exp (- data * parameter[3]) },
        parameter = c(0,10,1),
        na.region = 0,
        low.intensity = 0
    )
)

#a + bx^c
setClass("InverseNoiseModel",
    contains = "NoiseModel",
    prototype = prototype(
        new("VersionedBiobase",versions=c(classVersion("NoiseModel"),InverseNoiseModel="1.0.0")),
        f = function(data,parameter) { parameter[1] + parameter[2]*data^-parameter[3] },
        parameter = c(0,1,1),
        na.region = 0,
        low.intensity = 0
    )
)

# bx^c
setClass("InverseNoANoiseModel",
    contains = "NoiseModel",
    prototype = prototype(
        new("VersionedBiobase",versions=c(classVersion("NoiseModel"),InverseNoiseModel="1.0.0")),
        f = function(data,parameter) {  parameter[1]*data^-parameter[2] },
        parameter = c(1,1),
        na.region = 0,
        low.intensity = 0
    )
)

setClass("GeneralNoiseModel",
    contains = "NoiseModel",
    prototype = prototype(
        new("VersionedBiobase",versions=c(classVersion("NoiseModel"),
                                          GeneralNoiseModel="1.0.0")),
        f = function(data,parameter) {   },
        parameter = c(),
        na.region = 0,
        low.intensity = 0
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.NoiseModel.slots <- function(object) {
  NULL
}

.valid.NoiseModel <- function(object) {
	.valid.NoiseModel.slots(object)
}

setValidity("NoiseModel",.valid.NoiseModel)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

setMethod("initialize","NoiseModel",
    function(.Object,ibspectra=NULL,reporterTagNames=NULL,one.to.one=TRUE,
             min.spectra=10,max.n.proteins=50,plot=F,pool=FALSE,
             nm.col=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999"),nm.cmp=NULL,...) {
      if (!is.null(ibspectra))
      {
        if (is.null(reporterTagNames)) reporterTagNames <- reporterTagNames(ibspectra)
        if (!is.character(reporterTagNames))
          stop("reporterTagNames has to be of class character.")
        
        if (is.matrix(reporterTagNames)) {
          # a combination matrix is supplied
          if (nrow(reporterTagNames) != 2)
            stop("reporterTagNames has to be either a vector or a combination matrix (with two rows).")
          reporterTagNames.combn <- reporterTagNames
        } else {
          reporterTagNames.combn <- combn(reporterTagNames,2)
        }

        ri <- reporterIntensities(ibspectra,na.rm=FALSE)
        
        if (plot & !pool)
          par(mfrow=rep(ceiling(sqrt(ncol(reporterTagNames.combn))),2))

        # matrix to store fitted paramters of reporter combinations
        # will be averaged in the end
        all.parameters <-
          matrix(nrow=ncol(reporterTagNames.combn),ncol=length(parameter(.Object)))
        
        pool.ch1 <- c()
        pool.ch2 <- c()

        if (one.to.one) {
          
          # for each reporter combination, fit a noise function on the ratios
          for (ch_i in seq_len(ncol(reporterTagNames.combn))) {
            ch1 <- ri[,reporterTagNames.combn[1,ch_i]]
            ch2 <- ri[,reporterTagNames.combn[2,ch_i]]
            pool.ch1 <- c(pool.ch1,ch1)
            pool.ch2 <- c(pool.ch2,ch2)

            all.parameters[ch_i,] <-
              .fitNoiseFunction(ch1,ch2,noiseFunction(.Object),parameter(.Object),...)

             if (plot) {
                x.data <- log10(sqrt(ch1*ch2))
                plot(x.data,ch1/ch2,log="y",
                       main=sprintf("%s vs %s",reporterTagNames.combn[1,ch_i],
                       reporterTagNames.combn[2,ch_i]),ylim=c(0.2,5))
                .lines.nf(noiseFunction(.Object),parameter=all.parameters[ch_i,],xlim=range(x.data,na.rm=TRUE))
              }

          }
        } else {

          # scale the ratios of proteins down to 1, and learn noise model on
          # the differences in their spectra
          
          pg <- proteinGroup(ibspectra)
          pnp <- peptideNProtein(pg)
          pepnprot <- as.data.frame(pnp[pnp[,'peptide'] %in% 
              peptides(pg,protein=reporterProteins(pg),
                       specificity=REPORTERSPECIFIC),],
              stringsAsFactors=FALSE)
          
          pepnspec <-
            as.data.frame(.as.matrix(spectrumToPeptide(pg),
                                     c("spectrum","peptide")),stringsAsFactors=FALSE)
          pgdf <- merge(pepnprot,pepnspec)

          # only take non-NA spectra
          pgdf <- pgdf[pgdf$spectrum %in%
                       rownames(ri)[apply(ri,1,function(d) !any(is.na(d)))],]
          pgdf$protein.g <- as.character(pgdf$protein.g)
          
          t <- table(pgdf$protein.g)
          sel.proteins <-
            names(t)[t>=min.spectra][order(t[t>=min.spectra],decreasing=TRUE)]
          sel.proteins <- sel.proteins[seq_len(min(length(sel.proteins),max.n.proteins))]
          sel = pgdf$protein.g %in% sel.proteins
          cat(sprintf("%s proteins with more than %s spectra, taking top %s.\n",
                  length(unique(pgdf$protein.g)),
                  min.spectra,max.n.proteins))
          pgdf <- cbind(pgdf[sel,],ri[pgdf[sel,"spectrum"],])
          

          for (ch_i in seq_len(ncol(reporterTagNames.combn))) {
            ch1 <- c()
            ch2 <- c()
            pool.ch1 <- c(pool.ch1,ch1)
            pool.ch2 <- c(pool.ch2,ch2)
            # compute protein ratios
            protein.ratios <- estimateRatio(ibspectra=ibspectra, 
                channel1=reporterTagNames.combn[1,ch_i],
                channel2=reporterTagNames.combn[2,ch_i],
                protein=unique(pgdf$protein.g),combine=F)
            
            # bring protein ratio to ratio 1
            for (protein in rownames(protein.ratios)) {
              if (abs(protein.ratios[protein,"lratio"]) < 1) {
              ch1 <- c(ch1,pgdf[pgdf$protein.g==protein,reporterTagNames.combn[1,ch_i]])
              ch2 <- c(ch2,pgdf[pgdf$protein.g==protein,reporterTagNames.combn[2,ch_i]] * 
                  10^protein.ratios[protein,"lratio"])
              }
            }
            
            if (!pool) {
              all.parameters[ch_i,] <- .fitNoiseFunction(ch1,ch2,
                  noiseFunction(.Object),parameter(.Object),...)
            
             if (plot) {
                x.data <- log10(sqrt(ch1*ch2))
                plot(x.data,ch1/ch2,log="y",
                       main=sprintf("%s vs %s",reporterTagNames.combn[1,ch_i],
                       reporterTagNames.combn[2,ch_i]),ylim=c(0.2,5))
                .lines.nf(noiseFunction(.Object),parameter=all.parameters[ch_i,],xlim=range(x.data))
              }
              
              print(all.parameters[ch_i,])
            }
            pool.ch1 <- c(pool.ch1,ch1)
            pool.ch2 <- c(pool.ch2,ch2)
          }
        }
        if (!pool) {
          parameter(.Object) <- .averageParameters(.Object,all.parameters) 
        } else {
          parameter(.Object) <- .fitNoiseFunction(pool.ch1,pool.ch2,noiseFunction(.Object),parameter(.Object),...)
        }
        if (plot) {
          x.data <- log10(sqrt(pool.ch1*pool.ch2))
          plot(x.data,pool.ch1/pool.ch2,log="y",
               main="All channels",ylim=c(0.2,5))
          .lines.nf(noiseFunction(.Object),parameter=parameter(.Object),xlim=range(x.data,na.rm=TRUE),col=nm.col[1])
          if (!is.null(nm.cmp)) {
            if (is(nm.cmp,"NoiseModel")) nm.cmp <- c(nm.cmp)
            for (nm_i in seq_along(nm.cmp)) {
              .lines.nm(nm.cmp[[nm_i]],xlim=range(x.data,na.rm=TRUE),col=nm.col[nm_i+1])
            }
          }

        }

        sel.allna <- apply(is.na(ri),1,all)
        sel.anyna <- apply(is.na(ri),1,any)
        ri.anyna <- as.numeric(ri[!sel.allna & sel.anyna,])
        na.val <- log10(as.numeric(ri.anyna))
        na.val <- na.val[!is.na(na.val)]

        print(parameter(.Object))
        lowIntensity(.Object) <- quantile(unlist(ri[,reporterTagNames]),prob=0.05,na.rm=T)
        naRegion(.Object) <- na.val
      }
      callNextMethod(.Object)
    }
)

.lines.nm <- function(nm,xlim,col="red",parameter=NULL,lwd=2,...) {
  xx <- seq(from=xlim[1],to=xlim[2],length.out=100)
  lines(xx,10^(1.96*sqrt(variance(nm,xx))),col=col,lwd=lwd)
  lines(xx,10^(-1.96*sqrt(variance(nm,xx))),col=col,lwd=lwd)
}


.lines.nf <- function(f,xlim,col="red",parameter=NULL,lwd=2,...) {
  xx <- seq(from=xlim[1],to=xlim[2],length.out=100)
  lines(xx,10^(1.96*sqrt(0.5)*f(xx,parameter=parameter)),col=col,lwd=lwd,...)
  lines(xx,10^(-1.96*sqrt(0.5)*f(xx,parameter=parameter)),col=col,lwd=lwd,...)
}

setGeneric("NoiseModel",function(ibspectra,...) standardGeneric("NoiseModel"))

setMethod("NoiseModel",signature(ibspectra="IBSpectra"),
  function(ibspectra,type="exponential",...) {
    if (type=="exponential") new("ExponentialNoiseModel",ibspectra,...)
    else if (type=="exponential.noa") new("ExponentialNoANoiseModel",ibspectra,...)
    else if (type=="inverse") new("InverseNoiseModel",ibspectra,...)
  }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .fitNoiseFunction and .averageParameters.
### (called by initialize)

setGeneric(".fitNoiseFunction",function(data1,data2,f,initial.parameters,...)
  standardGeneric(".fitNoiseFunction"))
setGeneric(".averageParameters",function(x,parameters)
  standardGeneric(".averageParameters"))

setMethod(".fitNoiseFunction", 
    signature=c(data1="numeric",data2="numeric",f="function",
      initial.parameters="numeric"),
    function(data1,data2,f,initial.parameters,n.bins=30,min.n.bin=50,
        na.rm=is.null(set.na.to),set.na.to=NULL,...) {

      fit.f <- function(param,x,y) 
           -sum(dnorm(x,mean=0,sd=f(data=y,parameter=param),log=T))
      
      # remove points which are NA in both channels
      sel <- is.na(data1) & is.na(data2)
      data1 <- data1[!sel]; data2 <- data2[!sel]
      
      if (na.rm) {
        sel <- is.na(data1) | is.na(data2)
        data1 <- data1[!sel]; data2 <- data2[!sel]
      } else {
        if (is.null(set.na.to)) 
          set.na.to <- quantile(10^abs(log10(data1/data2)),probs=0.99,na.rm=T)
      }
      
      ch1 <- log10(data1)
      ch2 <- log10(data2)
      
      lratio <- ch1 - ch2
      avg <- 0.5*(ch1 + ch2) 
      
      if (!is.null(set.na.to)) {
        sel <- is.na(data1);
        lratio[sel] <- log10(set.na.to)
        avg[sel] <- ch2[sel]
        sel <- is.na(data2); 
        lratio[sel] <- -log10(set.na.to)
        avg[sel] <- ch1[sel]
      }
      
      result <- nlminb(initial.parameters,fit.f,
          lower=rep(1e-10,length(initial.parameters)),x=lratio,y=avg,...)
      
      result$par
    }
)


setMethod(".averageParameters",
    signature(x="ExponentialNoiseModel",parameters="matrix"),
    function(x,parameters) {
      parameter <- numeric()
      parameter[1] <- mean(parameters[,1])
      parameter[2] <- mean(parameters[,2])
      parameter[3] <-
        -log(mean(parameters[,2]*exp(-parameters[,3]))) + log(parameter[2])
      parameter
    }
)

setMethod(".averageParameters",
    signature(x="ExponentialNoANoiseModel",parameters="matrix"),
    function(x,parameters) {
      parameter <- numeric()
      parameter[1] <- mean(parameters[,1])
      parameter[2] <-
        -log(mean(parameters[,1]*exp(-parameters[,2]))) + log(parameter[1])
      parameter
    }
)


setMethod(".averageParameters",
    signature(x="InverseNoiseModel",parameters="matrix"),
    function(x,parameters) { # Coincidentally the same average formula as for the exponential model!
      parameter <- numeric()
      parameter[1] <- mean(parameters[,1])
      parameter[2] <- mean(parameters[,2])
      parameter[3] <-
        -log(mean(parameters[,2]*exp(-parameters[,3]))) + log(parameter[2])
      parameter
    }
)

setMethod(".averageParameters",
    signature(x="InverseNoANoiseModel",parameters="matrix"),
    function(x,parameters) { # Coincidentally the same average formula as for the exponential model!
      parameter <- numeric()
      parameter[1] <- mean(parameters[,1])
      parameter[2] <- -log(mean(parameters[,1]*exp(-parameters[,2]))) + log(parameter[1])
      parameter
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("variance",function(x,channel1,channel2) standardGeneric("variance"))
setGeneric("stddev",function(x,...) standardGeneric("stddev"))
setGeneric("noiseFunction",function(x) standardGeneric("noiseFunction"))
setGeneric("parameter",function(x) standardGeneric("parameter"))
setGeneric("parameter<-",function(x,value) standardGeneric("parameter<-"))
setGeneric("lowIntensity",function(x) standardGeneric("lowIntensity"))
setGeneric("lowIntensity<-",function(x,value) standardGeneric("lowIntensity<-"))
setGeneric("naRegion",function(x) standardGeneric("naRegion"))
setGeneric("naRegion<-",function(x,value) standardGeneric("naRegion<-"))

## could support switching to fitting of the square-root standard deviation 
##  instead of the standard deviations themselves (as done in limma):
#setMethod("noiseFunction","NoiseModel",function(x) function(...) x@f(...)**2)
## as far as it was tested, it had little impact on the fit of the noise model

setMethod("noiseFunction","NoiseModel",function(x) x@f)

setMethod("variance",signature(x="NoiseModel",channel1="numeric",channel2="missing"),
    function(x,channel1) 0.5*(noiseFunction(x)(data=channel1,parameter=parameter(x))^2)
)
setMethod("variance",signature(x="NoiseModel",channel1="numeric",channel2="numeric"),
    function(x,channel1,channel2) {
      var1 <- 0.5*(noiseFunction(x)(data=channel1,parameter=parameter(x))^2)
      var2 <- 0.5*(noiseFunction(x)(data=channel2,parameter=parameter(x))^2)
      var1 + var2
    }
)
setMethod("stddev",signature(x="NoiseModel"),
    function(x,...) sqrt(variance(x,...))
)


setMethod("parameter","NoiseModel",function(x) x@parameter)
setReplaceMethod("parameter","NoiseModel",function(x,value) {
      x@parameter <- value
      x }
)

setMethod("lowIntensity","NoiseModel",function(x) x@low.intensity)
setReplaceMethod("lowIntensity","NoiseModel",function(x,value) {
      x@low.intensity <- value
      x }
)

setMethod("naRegion","NoiseModel",function(x) x@na.region)
setReplaceMethod("naRegion","NoiseModel",function(x,value) {
      x@na.region <- value
      x }
)
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show",signature(object="NoiseModel"),
    function(object){
      cat(class(object),"\n")
      cat("  definition: ")
      print(noiseFunction(object))
      cat("  parameter:     [",paste(sprintf("%.2f",object@parameter),collapse="; "),"]\n")
      q <- quantile(object@na.region)
      cat("  na region:     [",paste(sprintf("%s",names(q)),sprintf("%.2f",q),collapse="; "),"]\n")
      cat("  low.intensity: ",sprintf("%.2f",object@low.intensity),"\n")
    }
)


plot.NoiseModel <- function(x,y=NULL,...) {
  args = list(...)
  if (is.null(args$xlim)) args$xlim = c(2,10)
  ss <- seq(from=args$xlim[0],to=args$xlim[1],length.out=100)
  nm.sd <- 1.96*sqrt(variance(x,ss))
  plot(ss,10^nm.sd,ylim=c(10^-max(nm.sd),10^max(nm.sd)),col="red",lwd=2,type="l",log="y",xlab="log10 intensity",ylab="ratio",...)
  lines(ss,10^(-nm.sd),col="red",lwd=2)
}
