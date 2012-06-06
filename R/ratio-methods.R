###  =========================================================================
### Ratio estimation and ratio distribution functions.
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Distribution fit functions
###
# Distributions used are Norm and Cauchy from package 'distr'.
# See class?Norm and class?Cauchy

fitWeightedNorm <- function(x,weights) {
  new("Norm",
      mean=weightedMean(data=x,weights=weights),
      sd=sqrt(weightedVariance(data=x,weights=weights))
  )
}

fitNorm <- function(x,portion=0.75) {
  x <- x[!is.na(x)]
  if (is.null(portion) || (portion <= 0)){
    bp <- boxplot.stats(x,coef=1)
    limit <- 0.5*(abs(bp$stats[1])+abs(bp$stats[5]))
  }
  else
    limit <- quantile(abs(x),prob=portion)
  good <- (x >= -limit) & (x <= limit)
  new("Norm",
      mean=mean(x[good]),
      sd=sd(x[good])
  )
}

fitGumbel <- function(x) {
  gumbel.fit <- function(theta,x){
    -sum(d(Gumbel(loc=theta[1],scale=theta[2]))(x,log=TRUE),na.rm=T)
  }
  good <- !is.na(x) & !is.nan(x)
  theta.start <- c(median(x[good]),IQR(x[good])/2)
  res <- nlminb(theta.start,gumbel.fit,x=x[good])
#                lower=c(-10,1e-20),upper=c(10,10)) 
  new("Gumbel",loc=res$par[1],scale=res$par[2])
}


fitCauchy <- function(x) {
  cauchy.fit <- function(theta,x){
    -sum(dcauchy(x,location=theta[1],scale=theta[2],log=TRUE),na.rm=T)
  }
  good <- !is.na(x) & !is.nan(x)
  theta.start <- c(median(x[good]),IQR(x[good])/2)
  res <- nlminb(theta.start,cauchy.fit,x=x[good],
                lower=c(-10,1e-20),upper=c(10,10)) 
  new("Cauchy",location=res$par[1],scale=res$par[2])
}

fitTd <- function(x) {
  t.fit <- function(theta,x){
    -sum(dt(x,df=theta[1],log=TRUE),na.rm=T)
  }
  good <- !is.na(x) & !is.nan(x)
  theta.start <- c(1) # TODO: find good starting value
  res <- nlminb(theta.start,t.fit,x=x[good])
  new("Td",df=res$par[1])
}



fitNormalCauchyMixture <- function(x) {
  gc.fit <- function(theta,x){
    -sum(log(theta[4]*dcauchy(x,location=theta[1],scale=theta[2])+(1-theta[4])*dnorm(x,mean=theta[1],sd=theta[3])))
  }
  good <- !is.na(x) & !is.nan(x)
  theta.start <- c(median(x[good]),IQR(x[good])/2,mad(x[good]),0.5)
  res <- nlminb(theta.start,gc.fit,x=x[good],
                lower=c(-10,1e-20,1e-20,0),upper=c(10,10,10,1)) 
  d1 <- Cauchy(location=res$par[1],scale=res$par[2])
  d2 <- Norm(mean=res$par[1],sd=res$par[3])
  dd <- UnivarMixingDistribution(d1,d2,mixCoeff=c(res$par[4],1-res$par[4]))
  list(mixture=dd,cauchy=d1,normal=d2)
}

fitGaussianMixture <- function(x,n=500) {
  x <- x[!is.na(x)]
  ## EM Algorithm
  ## Initialization: Guessing parameters
  theta <- list(
                weights = c(0.5,0.5),
                means = c(mean(x),mean(x)),
                sds = c(sd(x)-0.5*sd(x), sd(x)+0.5*sd(x)))

  em.gauss.mixd = function(x0,theta,same.mean=TRUE) {
    N <- length(x0)

    ## Expectation-Step
    ##prob for each datapoint to come from each distr
    p_iks <- do.call(cbind,lapply(seq_along(theta$weights),
                                  function(m) dnorm(x0,theta$means[m],theta$sds[m])*theta$weights[m]))
  
    ## membership weights
    Expectation <- p_iks / rowSums(p_iks)
    sum.weights <- colSums(Expectation)
  
    ## Maximization-Step
    theta$weights <- 1/N * colSums(Expectation) # new weights
    if (!same.mean)
      theta$means <- 1/sum.weights * colSums(Expectation*x0) # weighted mean
    theta$sds <- sqrt(1/sum.weights *
                      colSums(Expectation*(x - matrix(theta$means,byrow=T,nrow=N,ncol=2))^2))
    
    return(theta)
  }

  for(i in seq_len(n)){ theta=em.gauss.mixd(x,theta,same.mean=T) }

  d1 <- Norm(theta$means[1],theta$sds[1])
  d2 <- Norm(theta$means[2],theta$sds[2])
  dd <- UnivarMixingDistribution(d1,d2,mixCoeff=theta$weights)
  return(dd)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### estimateRatio generic and functions.
###
setGeneric("estimateRatioNumeric",
           function(channel1, channel2, noise.model, ...)
             standardGeneric("estimateRatioNumeric"))

setGeneric("estimateRatio",
           function(ibspectra, noise.model=NULL,
                    channel1, channel2, protein, peptide,...)
           standardGeneric("estimateRatio") )

## method definitions
setMethod("estimateRatioNumeric",
    signature(channel1="numeric",channel2="numeric",noise.model="NULL"),
    function(channel1,channel2,noise.model=NULL,...) 
      estimateRatioNumeric(channel1,channel2, ...)
      )
 
setMethod("estimateRatioNumeric",
    signature(channel1="numeric",channel2="numeric",noise.model="missing"),
    function(channel1,channel2,summarize.f=median, ...) {
      if (length(channel1) != length(channel2))
        stop("length of channel 1 does not equal length of channel 2")
      if (length(channel1)==0)
        return(c(lratio=NA))
      sel <- !is.na(channel1) & !is.na(channel2) & channel1 > 0 & channel2 > 0
      
      return(c(lratio=summarize.f(log10(channel1[sel])-log10(channel2[sel]))))
    }
)

setMethod("estimateRatioNumeric",signature(channel1="numeric",channel2="numeric",
                                           noise.model="NoiseModel"),
    function(channel1,channel2,noise.model,ratiodistr=NULL,
             variance.function="maxi",
             sign.level=0.05,sign.level.rat=sign.level,sign.level.sample=sign.level,
             remove.outliers=TRUE,outliers.coef=1,outliers.trim=0,n.sample=NULL, 
             method="isobar",fc.threshold=1.3,channel1.raw=NULL,channel2.raw=NULL,
             use.na=FALSE,preweights=NULL) {
      
      if (length(channel1) != length(channel2))
        stop("length of channel 1 does not equal length of channel 2")
      if (!is.null(channel1.raw) && length(channel1) != length(channel1.raw) ||
          !is.null(channel2.raw) && length(channel2) != length(channel2.raw)) 
        stop("length of orig channels [",length(channel1.raw),"] is not the same as channels [",length(channel1),"]")
      
      if (length(channel1)==0) {
        if (method == "isobar")
          return(c(lratio=NA,variance=NA,n.spectra=NA,n.na1=NA,n.na2=NA,
                   p.value.rat=NA,p.value.sample=NA,is.significant=NA))
        else if (method=="libra" | method == "pep")
          return(c(ch1=NA,ch2=NA))
        else if (method=="multiq")
          return(c(lratio=NA,variance=NA,n.spectra=0,isum=NA))
        else if (method=="test"|method=="compare.all")
          return(c(lratio=NA, variance=NA,
                   n.spectra=NA,
                   unweighted.ratio=NA,
                   is.sign.isobar    = NA,
                   is.sign.isobar.ev = NA,
                   is.sign.rat       = NA,
                   is.sign.sample    = NA,
                   is.sign.ttest     = NA,
                   ## is.significant.wttest    = res.isobar$is.significant,
                   is.sign.fc        = NA,
                   is.sign.fc.nw     = NA
                   ))
        else warning("no spectra, unknown method")
      }
      
      sel <- !is.na(channel1) & !is.na(channel2) & channel1 > 0 & channel2 > 0
      sel.notna <- !is.na(channel1) & !is.na(channel2) & channel1 > 0 & channel2 > 0
      ## Implementation of multiq and libra methods
      if (method=="multiq") {
        lratio <- log10(sum(channel1[sel])/sum(channel2[sel]))
        return(c(lratio, variance=NA,n.spectra=length(lratio),
                 isum=sum(channel1[sel]+channel2[sel])))
      }
      if (method=="libra" | method=="pep") {
        return(c(ch1=sum(channel1[sel]),ch2=sum(channel2[sel])))
      }

      i1 <- log10(channel1)
      i2 <- log10(channel2)

      ## Compute all the spectrum ratios and their individual
      ## variance based on the noise model
      log.ratio = i2[sel.notna] - i1[sel.notna]
      if (is.null(channel1.raw))
        var.i <- variance(noise.model,i1,i2)
      else 
        var.i <- variance(noise.model,
                          log10(channel1.raw),log10(channel2.raw))
      
      ## linear regression estimation
      if (method=="lm" || method=="compare.all") {
        fm <- lm(channel2~channel1+0,
                 weights=1/var.i)
        summary.fm <- summary.lm(fm)$coefficients
        ci <- confint(fm,level=1-sign.level.rat)
        res.lm <- c(lratio=log10(summary.fm[1]),
                    ratio=summary.fm[1],
                    stderr=summary.fm[2],
                    ratio.ci.l=ci[1],
                    ratio.ci.u=ci[2])

        if (!is.null(ratiodistr)) {
          rat.neg <- log10(res.lm['ratio']) < distr::q(ratiodistr)(0.5)
          ci[ci < 0] <- 0
          # TODO: fix confidence interval which goes to minus
          #message(paste(ci,collapse=":"))
          lrat.p <- log10(ifelse(rat.neg,ci[2],ci[1]))
          res.lm['p.value'] <- distr::p(ratiodistr)(lrat.p,lower.tail=rat.neg)
          res.lm['is.significant'] <- res.lm['p.value'] < sign.level.sample
        }
        
        if (method == "lm") return (res.lm)
      }      
       
      # First, compute ratios on spectra with intensities for both reporter ions
      lratio.n.var <-
        .calc.weighted.ratio(sel.notna,i1,i2,var.i,
                             remove.outliers,outliers.coef,outliers.trim,
                             variance.function,preweights)
      
      if (method == "ttest" || method == "compare.all") {
        if (length(log.ratio) < 2) p.value <- 1  else p.value <- t.test(log.ratio)$p.value
        res.ttest <- c(
                       lratio=lratio.n.var['lratio'],variance=NA,n.spectra=length(log.ratio),
            p.value=p.value,is.significant=p.value<sign.level)
        if (method != "compare.all") return(res.ttest)
      }
#if (method == "wttest" || method == "compare.all") {
#        p.value <- (length(log.ratio) < 2)? 1 : weighted.t.test(log.ratio,w=weights,weighted.ratio)$p.value
#        res.wttest <- c(
#            lratio=weighted.ratio,variance=NA,n.spectra=length(log.ratio),
#            p.value=p.value,is.significant=p.value<sign.level)
#        if (method != "compare.all") return(res.wttest)
#      }
      if (method == "fc" || method == "compare.all") {
        res.fc <- c(lratio=lratio.n.var['lratio'],variance=NA,
                    n.spectra=length(log.ratio),
                    is.significant=abs(lratio.n.var['lratio'])>log10(fc.threshold))
        ratio.nw <- mean(log.ratio,na.rm=TRUE)
        res.fc.nw <- c(lratio=ratio.nw,variance=NA,n.spectra=length(log.ratio),
            is.significant=abs(ratio.nw)>log10(fc.threshold))
        if (method != "compare.all") return(res.fc)
      }
 
      if (method == "isobar" || method == "compare.all") {
        # Check for channels where one is NA
        sel.ch1na <- is.na(channel1) & !is.na(channel2)
        sel.ch2na <- is.na(channel2) & !is.na(channel1)

        if (use.na) {
        # Require channel intensity to be outside of the 'NA region'
        sel.ch1na <- sel.ch1na & i2 > quantile(naRegion(noise.model),prob=0.99)
        sel.ch2na <- sel.ch2na & i1 > quantile(naRegion(noise.model),prob=0.99)

        # Set the value of that channel to a rather high one
        val <- quantile(naRegion(noise.model),prob=0.5)

        # Require the ratio to be higher than the previously observed one
        # TODO: handle case when we have no ratio
        if (!is.na(lratio.n.var['lratio'])) {
          if (lratio.n.var['lratio'] < 0) {
            
            sel.ch1na <- sel.ch1na & i2 - val < lratio.n.var['lratio']
            sel.ch2na <- sel.ch2na & val - i1 < lratio.n.var['lratio']
          } else {
            sel.ch1na <- sel.ch1na & i2 - val > lratio.n.var['lratio']
            sel.ch2na <- sel.ch2na & val - i1 > lratio.n.var['lratio']
          }
        }

#        i1[sel.ch1na] <- max(val,i2[sel.ch1na]-10)
#        i2[sel.ch2na] <- max(val,i1[sel.ch2na]-10)
#        i1.raw[sel.ch1na] <- max(val,i2[sel.ch1na]-10)
#        i2.raw[sel.ch2na] <- max(val,i1[sel.ch2na]-10)
        
        i1[sel.ch1na] <- val
        i2[sel.ch2na] <- val
        i1.raw[sel.ch1na] <- i2.raw[sel.ch1na]
        i2.raw[sel.ch2na] <- i1.raw[sel.ch2na]

        sel <- sel.notna | sel.ch1na | sel.ch2na
#        message("ch1na: ",sum(sel.ch1na),"; ch2na: ",sum(sel.ch2na))

        # calculate final ratio
        lratio.n.var <- 
          .calc.weighted.ratio(sel,i1,i2,var.i,
                               remove.outliers,outliers.coef,outliers.trim,
                               variance.function,preweights)
        } # use.na

        weighted.ratio <- as.numeric(lratio.n.var['lratio'])
        calc.variance <- as.numeric(lratio.n.var['calc.variance'])

        res.isobar <- 
          c(lratio=weighted.ratio, variance=calc.variance,
            n.spectra=length(log.ratio),
            n.na1=sum(sel.ch1na),n.na2=sum(sel.ch2na),
            p.value.rat=pnorm(weighted.ratio,mean=0,sd=sqrt(calc.variance),
                              lower.tail=weighted.ratio<0),
            p.value.sample=NA,is.significant=NA)
        
        if (!is.null(ratiodistr)) 
          res.isobar['p.value.sample'] <- 
            p(ratiodistr)(weighted.ratio,lower.tail=weighted.ratio<distr::q(ratiodistr)(0.5))

        if (method=="compare.all") 
          res.isobar['is.significant.ev'] <- 
            (res.isobar['p.value.sample'] <= sign.level.sample) &&
            pnorm(weighted.ratio,mean=distr::q(ratiodistr)(0.5),
                  sd=as.numeric(sqrt(lratio.n.var['estimator.variance'])),
                  lower.tail=weighted.ratio<distr::q(ratiodistr)(0.5)) <= sign.level.rat

        res.isobar['is.significant'] <- 
          (res.isobar['p.value.sample'] <= sign.level.sample) && 
          (res.isobar['p.value.rat'] <= sign.level.rat)

        if (method != "compare.all") return(res.isobar)
      }
      if (method != "compare.all") stop(paste("method",method,"not available"))
        return(c(lratio=weighted.ratio,
                 ratio=10^weighted.ratio,
                 ratio.lm=as.numeric(res.lm['ratio']),
                 variance=calc.variance,
                 var.lm=as.numeric(res.lm['stderr']**2),
                 n.spectra=length(log.ratio),
                 unweighted.ratio=mean(log.ratio,na.rm=TRUE),
                 is.sign.isobar    = as.numeric(res.isobar['is.significant']),
                 is.sign.isobar.ev = as.numeric(res.isobar['is.significant.ev']),
                 is.sign.lm        = as.numeric(res.lm['is.significant']),
                 is.sign.rat       = as.numeric(res.isobar['p.value.rat']<sign.level),
                 is.sign.sample    = as.numeric(res.isobar['p.value.sample']<sign.level),
                 is.sign.ttest     = as.numeric(res.ttest['is.significant']),
                 ## is.significant.wttest    = res.isobar$is.significant,
                 is.sign.fc        = as.numeric(res.fc['is.significant']),
                 is.sign.fc.nw     = as.numeric(res.fc.nw['is.significant'])
                 ))
    }
)
.get.ri <- function(ri,ch) {
  if (ch == "ALL")
    rowSums(ri,na.rm=TRUE)
  else if (ch == "AVG")
    apply(ri,1,mean,na.rm=TRUE)
  else
    ri[,ch]
}

                           
.calc.weighted.ratio <- function(sel,i1,i2,variance,
                                 remove.outliers,outliers.coef,outliers.trim,
                                 variance.function,preweights=NULL) {
  
  log.ratio <- i2[sel] - i1[sel]
  variance <- variance[sel]
  preweights <- preweights[sel]
  
  if (remove.outliers) {
    # TODO: Weighted outlier removal
    if (outliers.trim == 0) {
      # use box-plot method
      bp <- boxplot.stats(log.ratio,coef=outliers.coef)
      sel.or <- (log.ratio >= bp$stats[1]) & (log.ratio <=bp$stats[5])
    } else {
      # use trim method
      sel.or <- log.ratio > quantile(log.ratio,outliers.trim) & 
                log.ratio < quantile(log.ratio,1-outliers.trim)
    }
    log.ratio <- log.ratio[sel.or]
    variance <- variance[sel.or]
    preweights <- preweights[sel.or]
  }

  weights <- 1/variance
  if (!is.null(preweights))
    weights <- weights*preweights
  sum.weights <- sum(weights)
  weighted.ratio <- weightedMean(log.ratio,weights)

  # variance calculation
  estimator.variance <- 1/sum.weights
  if (length(log.ratio) == 1)
    sample.variance <- estimator.variance^0.75
  else
    if (length(log.ratio) == 2)
      sample.variance <- max(c(weightedVariance(log.ratio,weights,weighted.ratio),
                               estimator.variance^0.75))
    else
      sample.variance <- weightedVariance(log.ratio,weights,weighted.ratio)

  calc.variance <- switch(variance.function,
      maxi = max(estimator.variance,sample.variance,na.rm=T),
      ev = estimator.variance,
      wsv = sample.variance
  )

  return(c(lratio=weighted.ratio,
           estimator.variance=estimator.variance,
           sample.variance=sample.variance,
           calc.variance=calc.variance))
}

estimateRatioForProtein <- function(protein,ibspectra,noise.model,channel1,channel2,
        combine=TRUE,method="isobar",specificity=REPORTERSPECIFIC,quant.w.grouppeptides=NULL,...) {
      if (combine) {
        if (method == "multiq" || method == "libra" || method=="pep") {
          ## first compute peptide ratios, summarize then
          peptide.ratios <-
            estimateRatioForPeptide(peptides(proteinGroup(ibspectra),protein=protein),
                                    ibspectra,noise.model=noise.model,
                                    channel1=channel1,channel2=channel2,
                                    combine=FALSE,method=method)

          if (method == "libra" | method=="pep") {
            # normalize to sum
            sum.i <- peptide.ratios$ch1+peptide.ratios$ch2
            p.ch1 <- peptide.ratios$ch1/sum.i
            p.ch2 <- peptide.ratios$ch2/sum.i
            if (method=="libra") {
              ##remove outliers +2sd
              sd1 <- sd(p.ch1)
              sd2 <- sd(p.ch2)
              
              p.ch1 <- p.ch1[p.ch1 > mean(p.ch1)-2*sd1 & p.ch1 < mean(p.ch1)+2*sd1]
              p.ch2 <- p.ch2[p.ch2 > mean(p.ch2)-2*sd2 & p.ch2 < mean(p.ch2)+2*sd2]
            }
            return(lratio=log10(p.ch2/p.ch1),variance=var(log10(p.ch2/p.ch1)),
                   ch1=p.ch1,ch2=p.ch2)
             
          }
          if (method == "multiq") {
            ## TODO: filtration of peptides with low confidence score
            ## TODO: Dynamic range comp.

            ## weighted average of peptide ratios
            return(c(
              lratio=weightedMean(peptide.ratios$lratio,weights=peptide.ratios$isum),
              variance=var(peptide.ratios$lratio)))
          }

        } else if (method %in% c("isobar","lm","ttest","fc","compare.all")) {
          .call.estimateRatio(protein,"protein",ibspectra,
                             noise.model,channel1,channel2,method=method,
                             specificity=specificity,...)
        } else {
          stop("method ",method," not known")
        }
      } else {
        res <- do.call(rbind,lapply(protein,function(individual.protein) {
             if (individual.protein %in% quant.w.grouppeptides) 
               specificity <- c(GROUPSPECIFIC,specificity)

             .call.estimateRatio(individual.protein,"protein",ibspectra,
                                noise.model,channel1,channel2,method=method,...,
                                specificity=specificity)
            }
        ))
        rownames(res) <- protein
        res
        #res[apply(res,2,!function(r) all(is.na(r))),]
      }
    }

estimateRatioForPeptide <- function(peptide,ibspectra,noise.model,channel1,channel2,combine=TRUE,...) {
      if (combine) {
        r <- .call.estimateRatio(peptide,"peptide",ibspectra,noise.model,
                                 channel1,channel2,...)
      } else {
        if (is.matrix(peptide)) {
          r <- t(apply(peptide,1,function(individual.peptide) 
                  .call.estimateRatio(matrix(individual.peptide,ncol=2),"peptide",ibspectra,noise.model,
                                      channel1,channel2,...)))
        
        } else {
          r <- t(sapply(peptide,function(individual.peptide) 
                        .call.estimateRatio(individual.peptide,"peptide",ibspectra,noise.model,
                                            channel1,channel2,...)))
        }
      }
      attr(r,"input") <- peptide
      return(r)
}


### Handling NULL protein or peptide argument

setMethod("estimateRatio",
          signature(ibspectra="IBSpectra",noise.model="ANY",
                    channel1="missing",channel2="missing",
                    protein="character",peptide="missing"),
          function(ibspectra,noise.model=NULL,protein,val="lratio",summarize=FALSE,...) {
            channels <- reporterTagNames(ibspectra)
            res <- matrix(NA,nrow=length(channels),ncol=length(channels),dimnames=list(channels,channels))
            for(channel1 in channels)
              for (channel2 in channels) {
                rat <- estimateRatio(ibspectra,noise.model=noise.model,
                                     channel1=channel1,channel2=channel2,
                                     protein=protein,...)
                res[channel1,channel2] <- rat[val]
              }
            return(res)
          }
)


setMethod("estimateRatio",
    signature(ibspectra="IBSpectra",noise.model="ANY",
        channel1="character",channel2="character",
        protein="character",peptide="NULL"),
    function(ibspectra,noise.model,channel1,channel2,protein,peptide=NULL,...) {
        estimateRatio(ibspectra,noise.model,channel1,channel2,protein=protein,...)
    }
)

setMethod("estimateRatio",
    signature(ibspectra="IBSpectra",noise.model="ANY",
              channel1="character",channel2="character",
              protein="NULL",peptide="character"),
    function(ibspectra,noise.model,channel1,channel2,protein=NULL,peptide,...) {
      estimateRatio(ibspectra,noise.model,channel1,channel2,peptide=peptide,...)
    }
)

setMethod("estimateRatio",
    signature(ibspectra="IBSpectra",noise.model="ANY",
              channel1="character",channel2="character",
              protein="NULL",peptide="matrix"),
    function(ibspectra,noise.model,channel1,channel2,protein=NULL,peptide,...) {
      estimateRatio(ibspectra,noise.model,channel1,channel2,peptide=peptide,...)
    }
)

### Calling estimateRatioForPeptide and calcRatioForProtein, resp.
setMethod("estimateRatio",
    signature(ibspectra="IBSpectra",noise.model="ANY",
        channel1="character",channel2="character",
        protein="character",peptide="missing"),
    function(ibspectra,noise.model,channel1,channel2,protein,...) {
        estimateRatioForProtein(protein,ibspectra,noise.model,channel1,channel2,...)
    }
)

setMethod("estimateRatio",
    signature(ibspectra="IBSpectra",noise.model="ANY",
              channel1="character",channel2="character",
              protein="missing",peptide="character"),
    function(ibspectra,noise.model,channel1,channel2,peptide,...) {
      estimateRatioForPeptide(peptide,ibspectra,noise.model,channel1,channel2,...)
    }
)

setMethod("estimateRatio",
    signature(ibspectra="IBSpectra",noise.model="ANY",
              channel1="character",channel2="character",
              protein="missing",peptide="matrix"),
    function(ibspectra,noise.model,channel1,channel2,peptide,...) {
      estimateRatioForPeptide(peptide,ibspectra,noise.model,channel1,channel2,...)
    }
)

### Helper function to estimateRatioNumeric
.call.estimateRatio <- function(x,level,ibspectra,noise.model,
                                channel1,channel2,
                                specificity=REPORTERSPECIFIC,modif=NULL,
                                n.sample=NULL,groupspecific.if.same.ac=FALSE,
                                use.precursor.purity=FALSE,do.warn=TRUE,...) {
  allowed.channels <- c(reporterTagNames(ibspectra),'AVG','ALL')
  if (is.null(channel1) || is.null(channel2))
    stop("channel1 and channel2 must not be NULL, but one of [",paste(allowed.channels,collapse=", "),"] !")
  #if (length(channel1) == 0 || length(channel1) > 1 || length(channel2) == 0 || length(channel2) > 1)
  #  stop("channel1 and channel2 must be of length one! Lengths: [",length(channel1),",",length(channel2),"]")
  if (!all(channel1 %in% allowed.channels) || !all(channel2 %in% allowed.channels))
    stop("channel1 and channel2 must be one of the reporter names: \n\t",paste(allowed.channels,collapse=", "),".")
  
  if (length(channel1) > 1 || length(channel2) > 1) {
    res  <- c()
    for (c1 in channel1) {
      for (c2 in channel2) {
        res <- rbind(res,data.frame(channel1=c1,channel2=c2,
                                    t(as.data.frame(.call.estimateRatio(x,level,ibspectra,noise.model,channel1=c1,channel2=c2,
                                                                        specificity,modif,n.sample,groupspecific.if.same.ac,
                                                                        use.precursor.purity,do.warn=do.warn,...))),
                                    stringsAsFactors=FALSE))
      }
    }
    rownames(res) <- NULL
    return(res)
  }
  if (level=="protein")
    sel <- spectrumSel(ibspectra,protein=x,specificity=specificity,do.warn=do.warn,
                       modif=modif,groupspecific.if.same.ac=groupspecific.if.same.ac)
  if (level=="peptide") 
    sel <- spectrumSel(ibspectra,peptide=x,modif=modif,do.warn=do.warn)
  
  ri <- reporterIntensities(ibspectra)[sel,,drop=FALSE]
  ri.raw <- reporterData(ibspectra,element="ions_not_normalized")[sel,,drop=FALSE]
  if (use.precursor.purity && .SPECTRUM.COLS['PRECURSOR.PURITY'] %in% colnames(fData(ibspectra)))
    precursor.purity <- fData(ibspectra)[sel,"precursor.purity"]
  else
    precursor.purity <- NULL

  ## TODO: implement n.sample
  i1 <- .get.ri(ri,channel1)
  i2 <- .get.ri(ri,channel2)

  if (is.null(ri.raw)) {
    estimateRatioNumeric(channel1=i1,channel2=i2,noise.model=noise.model,...,preweights=precursor.purity)
  } else {
    estimateRatioNumeric(channel1=i1,channel2=i2,noise.model=noise.model,...,preweights=precursor.purity,
                  channel1.raw=.get.ri(ri.raw,channel1),
                  channel2.raw=.get.ri(ri.raw,channel2))
  }
}


getMultUnifPValues <- function(product,pvals=NULL,n=NULL){
  
  if (is.null(n))
    return(NULL)
  
  if (!is.null(pvals))
    if (n == length(pvals))
      product <- prod(pvals)
    else
      return(NULL)
  
  if (n == 1)
    return(product)
  
  s <- 1
  for (i in 1:(n-1))
    s <- s + (-log(product))^i/factorial(i)
  product*s
  
} # getMultUnifPValues


getMultUnifDensity <- function(x,n=NULL){
  
  if (is.null(n))
    return(NULL)
  if (n == 1)
    return(1)
  
  i <- n-1
  (-log(x))^i/factorial(i)
  
} # getMultUnifDensity


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### proteinRatios generic and function.
###

## combn.matrix generates pairwise combinations 
##  of the values in 'x':
##  - all against all (using combn, methd="global")
##  - all combinations across classes (method="interclass")
##     used when computing ratios class A against class B
##  - all combinations within classes (method="intraclass")
##     used when computing a intra-class ratio distribution 
##     for estimation of biological variability
##  - method "versus.class": all combinations against class vs
##  - method "verus.channel": all combinations against channel vs
combn.matrix <- function(x,method="global",cl=NULL,vs=NULL) {
  if (method != "global" & is.null(cl))
    stop("No class labels cl provided.")
  if (method != "global" & !is.character(cl)) {
    warning("Class labels should be of type character - using as.character to convert.")
    cl <- as.character(cl)
  }
  if (method != "global" & length(cl) != length(x))
    stop(sprintf("cl argument does not have the same length as x (cl: %s, x: %s).",
                 length(cl),length(x)))

  if (!is.null(vs) && !is.character(vs))
    vs <- as.character(vs)

  # create a combn matrix with all combinations of channels to consider 
  if (method == "versus.class" || method == "versus.channel") {
    if (is.null(vs)) stop("vs argument may not be null when method is versus")
    if (method == "versus.channel") {
      if (!vs %in% x) stop("vs argument must be one of [",paste(x,collapse=", "),"]")
      pos <- which(x==vs)
      cmbn <- rbind(vs,x[-pos])
      if (!is.null(cl)) {
        vs.class <- cl[pos]
        cmbn <- rbind(cmbn,vs.class,cl[-pos])
      }
    }
    if (method == "versus.class") {
      if (!all(vs %in% cl)) stop("vs argument must be one of [",paste(cl,collapse=", "),"]")
      if (is.null(cl)) stop("class labels must be given with method versus.class")
      pos <- which(cl==vs)
      cmbn <- rbind(x[pos],rep(x[-pos],each=length(pos)))
      cmbn <- rbind(cmbn,vs,rep(cl[-pos],each=length(pos)))

    }
  } else if (method == "global") {
    cmbn <- combn(x,2) # take all combinations
    if (!is.null(cl))
      cmbn <- rbind(cmbn,combn(cl,2))
  } else {
    x <- x[!is.na(cl)]
    cl <- cl[!is.na(cl)]
    t <- table(cl)[unique(cl)]

    if (method == "intraclass") {
      if (length(t[t==1]) > 0)
        warning("Some class labels are not repeated - ",
                "those are ignored in intraclass ratios.")
          
      cmbn <- do.call(cbind,lapply(names(t)[t>1],
                                    function(xx) rbind(combn(x[which(cl==xx)],2),class1=xx,class2=xx)))
      
    } else if (method == "interclass") {
      if (length(t) == 1) {
        warning("Cannot compute interclass ratios when there is only one class - taking ratios vs ALL")
        cmbn <- matrix(c(x,rep("ALL",length(x))),nrow=2,byrow=TRUE)
      } else {
        cmbn <- matrix(nrow=4,ncol=0)
        for (name in names(t)) {
          pos=which(cl==name);posn=which(cl!=name);
          for (i in pos) 
            for (j in posn) {
              cc <- c(x[i],x[j],cl[i],cl[j])
              if (ncol(cmbn) > 0 &
                  any(apply(cmbn,2,function(xx) identical(rev(cc[1:2]),xx[1:2]))))
                next;
              cmbn <- cbind(cmbn,cc)
            }
        }
      }
    } else {
      stop(paste("method",method,"not implemented."))
    }
  }
  if (nrow(cmbn) == 2) rownames(cmbn) <- c("r1","r2")
  else if (nrow(cmbn) == 4) rownames(cmbn) <- c("r1","r2","class1","class2")
  return(cmbn)
}

## create a table with all protein ratios
combn.protein.tbl <- function(ibspectra,noise.model,ratiodistr,
                              proteins=NULL,cmbn,peptide=NULL,modif=NULL,
                              symmetry=FALSE,reverse=FALSE,variance.function="maxi",...) {
  
  ratios <- do.call(rbind,apply(cmbn,2,function(x) {
    if (reverse)
      if (length(x) == 4)
        x <- x[c(2,1,4,3)]
      else
        x <- rev(x)

    r <- estimateRatio(ibspectra=ibspectra,noise.model=noise.model,
                       channel1=x[1],channel2=x[2],protein=proteins,
                       peptide=peptide,modif=modif,
                       ratiodistr=ratiodistr,variance.function=variance.function,...)
    if (class(r)=="numeric") {
      r <- t(r)
      rownames(r) <- "prot1"
    }
    df <- as.data.frame(r,stringsAsFactors=FALSE)
    df$ac <- rownames(df)

    if (is.matrix(attr(r,"input")))
      df <- cbind(as.data.frame(attr(r,"input"),stringsAsFactors=FALSE),df)
    
    #rownames(df) <- paste(df$ac,paste(x[2],x[1],sep=":"),sep="_")
    rownames(df) <- NULL
    df$r1 <- x[1]
    df$r2 <- x[2]
    if (length(x) == 4) {
      df$class1 <- x[3]
      df$class2 <- x[4]
    }
    return(df)
  }))
  ratios <- ratios[order(ratios$ac,ratios$r1,ratios$r2),]
  
  if (symmetry) {
    ratios.inv <- ratios
    ratios.inv[,'lratio'] <- -ratios.inv[,'lratio']
    r1 <- ratios.inv[,'r1']
    ratios.inv[,'r1'] <- ratios.inv[,'r2']
    ratios.inv[,'r1'] <- r1
    ratios <- rbind(ratios,ratios.inv)
  }
  return(ratios)  
}


peptideRatios <- function(ibspectra,...,protein=NULL,peptide=peptides(proteinGroup(ibspectra))) {
  proteinRatios(ibspectra,...,proteins=protein,peptide=peptide)
}

ratiosReshapeWide <- function(quant.tbl,grouped.cols=TRUE,vs.class=NULL,sep=".") {
  attrs <- attributes(quant.tbl)

  classes.unique <- "class1" %in% colnames(quant.tbl) &&
                    !any(is.na(quant.tbl$class1)) && !any(is.null(quant.tbl$class1)) && 
                    all(table(unique(quant.tbl[,c("r1","class1")])$class1)==1) &&
                    all(table(unique(quant.tbl[,c("r2","class2")])$class2)==1)
  if (!is.null(vs.class)) {
    if (length(vs.class)==1)
      quant.tbl$comp <- paste(quant.tbl$class2)
    else 
      quant.tbl$comp <- paste(quant.tbl$class2,quant.tbl$class1,sep="/")

    quant.tbl <- subset(quant.tbl,class1 %in% vs.class)
  } else {
    if (classes.unique) {
      quant.tbl$comp <- paste(quant.tbl$class2,quant.tbl$class1,sep="/")
    } else {
      quant.tbl$comp <- paste(quant.tbl$r2,quant.tbl$r1,sep="/")
    }
  }
  quant.tbl  <- quant.tbl[,-(c(which(colnames(quant.tbl) %in% c("r1","r2","class1","class2"))))]
  v.names <- c("lratio","variance","n.spectra","p.value.rat","p.value.sample","is.significant","sd","n.na1","n.na2")
  v.names <- v.names[v.names %in% colnames(quant.tbl)]
  if ("n.pos" %in% colnames(quant.tbl)) v.names <- c(v.names,"n.pos")
  if ("n.neg" %in% colnames(quant.tbl)) v.names <- c(v.names,"n.neg")
  timevar <- "comp"
  idvar <- colnames(quant.tbl)[!colnames(quant.tbl) %in% c(timevar,v.names)]

  res <- reshape(quant.tbl,v.names=v.names,idvar=idvar,timevar=timevar,direction="wide",sep=sep)
  if (grouped.cols) {
    col.order <- c(idvar,
                   unlist(lapply(v.names,function(n) grep(n,colnames(res),fixed=TRUE,value=TRUE))))
    res <- res[,col.order]
  }
  for (a in names(attrs)) 
    if (!a %in% c("row.names","names","class")) attr(res,a) <- attrs[[a]]
  res
}

proteinRatios <-
  function(ibspectra,noise.model,
           reporterTagNames=NULL,
           proteins=reporterProteins(proteinGroup(ibspectra)),peptide=NULL,
           cl=classLabels(ibspectra),
           method="global",symmetry=FALSE,
           summarize=FALSE,summarize.method="mult.pval",
           min.detect=NULL,strict.sample.pval=TRUE,strict.ratio.pval=TRUE,orient.div=0,
           sign.level=0.05,sign.level.rat=sign.level,sign.level.sample=sign.level,
           ratiodistr=NULL,variance.function="maxi",
           combine=FALSE,p.adjust=NULL,reverse=FALSE,
           combn=NULL,...) {

    if ((!is.null(proteins) && !is.null(peptide)) ||
        (is.null(proteins) && is.null(peptide)))
      stop("supply either protein or peptides!")      
    
    if (!is.null(p.adjust) && !p.adjust %in% p.adjust.methods)
      stop("p.adjust parameter must be one of '",paste(p.adjust.methods,collapse="','"),"'")
    
    #if (summarize && method=="interclass" && length(unique(cl[!is.na(cl)])) != 2)
    #  stop("Usage of inter-class ratios with summarize when having more than two classes is not supported ATM, sorry.",
    #       " class labels: ",paste(cl,collapse=", "))

    if (is.null(reporterTagNames)) reporterTagNames <- reporterTagNames(ibspectra)
    if (is.null(cl) && is.null(combn)) stop("please supply class labels as argument cl or a combn matrix")

    if (is.null(combn))
      combn <- combn.matrix(reporterTagNames,method,cl)

    if (ncol(combn) < 1) 
      stop("No possible combination for reporters [",paste(reporterTagNames,sep=","),"]",
           " w/ classes [",paste(cl,sep=","),"] and method ",method," possible.",
           " summarize=",ifelse(summarize,"TRUE","FALSE"))
    
    ratios <- combn.protein.tbl(ibspectra,noise.model,ratiodistr,proteins,combn,
                                peptide=peptide,symmetry=symmetry,reverse=reverse,
                                variance.function=variance.function,combine=combine,...)

    if (summarize) {
      if (method=="global")
        stop("summarization not meaningful with method='global'. ",
             "Use method='intraclass' or method='interclass' to use ratios in or between classes.")

      #n.combination <- length(unique(combn[1,]))
      n.combination <- table(combn["class1",],combn["class2",])
      if (max(n.combination) < 2)
        stop("Summarize=TRUE makes no sense with only one combination, set summarize to FALSE or class labels differently.")

      if (is.null(min.detect))
        min.detect <- n.combination
      
      df <- summarize.ratios(ratios,summarize.method,min.detect,n.combination,
                             strict.sample.pval,strict.ratio.pval,orient.div,
                             sign.level,sign.level.rat,sign.level.sample,
                             variance.function=variance.function,
                             ratiodistr=ratiodistr)

      attributes(df) = c(attributes(df),list(
          classLabels=cl,combn.method=method,symmetry=symmetry,
          summarize=TRUE,summarize.method=summarize.method,
          min.detect=min.detect,
          strict.sample.pval=strict.sample.pval,
          strict.ratio.pval=strict.ratio.pval,
          orient.div=orient.div,
          sign.level.rat=sign.level.rat,
          sign.level.sample=sign.level.sample,
          ratiodistr=ratiodistr,
          variance.function=variance.function,
          combine=combine,p.adjust=p.adjust,reverse=reverse))

      return(df)
    } else {

      if (!is.null(p.adjust)) {
        stop("p.adjust argument is not working")
        ##ratios$p.value.rat <- p.adjust(ratios$p.value.rat,p.adjust)
      }

      attributes(ratios) = c(attributes(ratios),list(
          classLabels=cl,combn.method=method,symmetry=symmetry,
          sign.level.rat=sign.level.rat,sign.level.sample=sign.level.sample,
          ratiodistr=ratiodistr,variance.function=variance.function,
          combine=combine,p.adjust=p.adjust,reverse=reverse))


      ratios
    }      
}

summarize.ratios <-
  function(ratios,summarize.method="mult.pval",min.detect=NULL,n.combination=NULL,strict.sample.pval=TRUE,
           strict.ratio.pval=TRUE,orient.div=0,
           sign.level=0.05,sign.level.rat=sign.level,sign.level.sample=sign.level,
           variance.function="maxi",ratiodistr=NULL,p.adjust=NULL) {

    if (!is.null(p.adjust) && !p.adjust %in% p.adjust.methods)
      stop("p.adjust parameter must be one of '",paste(p.adjust.methods,collapse="','"),"'")
    if (!summarize.method  %in% c("mult.pval","mean"))
      stop("Implemented summarize.methods: mult.pval, mean. Please choose one of them instead of ",summarize.method,".")
    if (!all(c("class1","class2") %in% colnames(ratios)))
      stop("ratios must specify classes w/ columns class1 and class2!")

    classes <- unique(ratios[,c("class1","class2")])
    
    if (is.null(n.combination)) {
      cc <- unique(ratios[c("r1","r2","class1","class2"),])
      n.combination <- table(cc["class1",],cc["class2",])
    }
    if (is.null(min.detect))
      min.detect <- n.combination

    mean.r <- ifelse(is.null(ratiodistr),0,distr::q(ratiodistr)(0.5))
    if (summarize.method == "mult.pval") {
#      usable <- ratios$p.value.rat <= orient.p.value & !is.na(ratios$lratio)
#      val.acs <- unique(ratios$ac[usable])
      
      result <- do.call(rbind,lapply(unique(ratios$ac),function(ac) {
        ac.sel.1 <- (ratios$ac==ac) & !is.na(ratios$lratio)
        do.call(rbind,lapply(seq_len(nrow(classes)),function(class_i) {
          class1 <- classes[class_i,1]
          class2 <- classes[class_i,2]
          n.combination.c <- ifelse(is.matrix(n.combination),n.combination[class2,class1],n.combination)
          min.detect.c <- ifelse(is.matrix(min.detect),min.detect[class2,class1],min.detect)
          ac.sel <-ac.sel.1 & ratios$class1 == class1 & ratios$class2 == class2
          if (!any(ac.sel)) {
            ## no data for ac
            return(list(ac=ac,lratio=NA,variance=NA,
                              n.spectra=0,n.pos=0,n.neg=0,
                              p.value.rat=1,p.value.sample=1,
                              is.significant=FALSE,r1=class1,r2=class2))
          }
          
          n.pos <- sum(ratios$lratio[ac.sel]>mean.r,na.rm=T)
          n.neg <- sum(ratios$lratio[ac.sel]<mean.r,na.rm=T)
          is.pos <- (n.pos > n.neg && n.neg <= orient.div)
          is.neg <- (n.neg > n.pos && n.pos <= orient.div)

          ## ratio summarization
          good.sel <- ac.sel
          if (!strict.ratio.pval) {
            ## take only positive/negative spectra for ratio summarization
            if (is.pos) 
              good.sel <- good.sel & ratios$lratio>mean.r
            if (is.neg)
              good.sel <- good.sel & ratios$lratio<mean.r
          }

          ac.ratios <- ratios$lratio[good.sel]
          ac.vars <- ratios$variance[good.sel]
          
          sample.var <- weightedVariance(ac.ratios,weights=1/ac.vars)
          estim.var <- 1/sum(1/ac.vars)
          
          variance <- switch(variance.function,
                             maxi = max(estim.var,sample.var,na.rm=T),
                             ev = estim.var,
                             wsv = sample.var
                             )  
          lratio <- weightedMean(ac.ratios,weights=1/ac.vars)
          
          if (is.null(ratiodistr))
            return(list(ac=ac,lratio=lratio,variance=variance,
                              n.spectra=min(ratios$n.spectra[ac.sel]),
                              n.pos=n.pos,n.neg=n.neg,p.value.rat=NA,p.value.sample=NA,
                              is.significant=NA,r1=class1,r2=class2))

          ## p-value computation
          p.value.rat <-
            pnorm(lratio,mean=mean.r,sd=sqrt(variance),lower.tail=lratio<mean.r)

          if (!is.null(p.adjust)) {
            p.value.rat <- p.adjust(p.value.rat,p.adjust)
          }
          
          product.p.vals <-
            prod(p(ratiodistr)(ac.ratios,lower.tail=n.neg>n.pos))
        
          if (strict.sample.pval) {
            p.value.sample <- getMultUnifPValues(product.p.vals*0.5^(n.combination.c-sum(ac.sel)),n=n.combination.c)
          } else {
            p.value.sample <- getMultUnifPValues(product.p.vals,
                                                 n=sum(ac.sel))
          }

          ## significance
          is.significant <- (p.value.sample <= sign.level.sample) &&
            (p.value.rat <= sign.level.rat) &&
            (sum(ac.sel) >= min.detect.c) &&
            (is.pos | is.neg)

          return(list(ac=ac,lratio=lratio,variance=variance,
                      n.spectra=min(ratios$n.spectra[ac.sel]),
                      n.pos=n.pos,n.neg=n.neg,
                      p.value.rat=p.value.rat,p.value.sample=p.value.sample,
                      is.significant=is.significant,r1=class1,r2=class2))
        }))
      }))
      
      return(do.call(rbind,apply(result,1,as.data.frame,stringsAsFactors=FALSE)))
      
    } else if (summarize.method=="mean") {
      do.call(rbind,lapply(unique(ratios$ac),function(ac) {
        ac.sel <- (ratios$ac==ac) & !is.na(ratios$lratio)
        if (!any(ac.sel)) {
            return(data.frame(ac=ac,lratio=NA,stringsAsFactors=FALSE))
        }
        return(data.frame(ac=ac,lratio=mean(ratios$lratio[ac.sel]),stringsAsFactors=FALSE))

      }))

    } else {
      stop (paste("summarize method",summarize.method,"not implemented."))
    }
  }
