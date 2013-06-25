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

fitTlsd <- function(x) {
  t.fit <- function(theta,x){
    #-sum(dtls(x,df=theta[3],location=theta[1],scale=theta[2],log=TRUE),na.rm=T)
    -sum(log(dtls(x,df=theta[3],location=theta[1],scale=theta[2])),na.rm=T)
  }
  dtls <- function(x,df,location,scale,log = FALSE)
            1/scale * dt((x - location)/scale, df, log = log)

  good <- !is.na(x) & !is.nan(x)
  theta.start <- c(median(x[good]),sd(x[good]),2) # TODO: find good starting value
  res <- nlminb(theta.start,t.fit,x=x[good],lower=c(-1,10^(-6),2),upper=c(1,10,100))
  new("Tlsd",df=res$par[3],location=res$par[1],scale=res$par[2])
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
             remove.outliers=TRUE,outliers.args=list(method="iqr",outliers.coef=1.5),
             method="isobar",fc.threshold=1.3,channel1.raw=NULL,channel2.raw=NULL,
             use.na=FALSE,preweights=NULL) {
      
      if (length(channel1) != length(channel2))
        stop("length of channel 1 does not equal length of channel 2")
      if (!is.null(channel1.raw) && length(channel1) != length(channel1.raw) ||
          !is.null(channel2.raw) && length(channel2) != length(channel2.raw)) 
        stop("length of orig channels [",length(channel1.raw),"] is not the same as channels [",length(channel1),"]")

      is.method <- function(m) identical(method,m)
      is.a.method <- function(m) m %in% method || is.method("compare.all")
      
      if (length(channel1)==0) {
        if (is.method("isobar"))
          return(c(lratio=NA,variance=NA,n.spectra=NA,n.na1=NA,n.na2=NA,
                   p.value.rat=NA,p.value.sample=NA,is.significant=NA))
        else if (is.method("libra") || is.method("pep"))
          return(c(ch1=NA,ch2=NA))
        else if (is.method("multiq"))
          return(c(lratio=NA,variance=NA,n.spectra=0,isum=NA))
        else if (is.method("test") || is.method("compare.all") || length(method) > 1)
          return(NULL)
        else warning("no spectra, unknown method")
      }
      
      sel <- !is.na(channel1) & !is.na(channel2) & channel1 > 0 & channel2 > 0
      sel.notna <- !is.na(channel1) & !is.na(channel2) & channel1 > 0 & channel2 > 0
      ## Implementation of multiq and libra methods
      if (is.method("multiq")) {
        lratio <- log10(sum(channel1[sel])/sum(channel2[sel]))
        return(c(lratio, variance=NA,n.spectra=length(lratio),
                 isum=sum(channel1[sel]+channel2[sel])))
      }
      if (is.method("libra") || is.method("pep")) {
        return(c(ch1=sum(channel1[sel]),ch2=sum(channel2[sel])))
      }

      i1 <- log10(channel1)
      i2 <- log10(channel2)

      ## Compute all the spectrum ratios and their individual
      ## variance based on the noise model
      log.ratios = i2 - i1
      if (is.null(channel1.raw))
        var.i <- variance(noise.model,i1,i2)
      else 
        var.i <- variance(noise.model,
                          log10(channel1.raw),log10(channel2.raw))
 
      if (remove.outliers) {
        if (identical(outliers.args$method,"wtd.iqr")) 
          outliers.args$weights <- 1/var.i
        outliers.args$log.ratio <- log.ratios
        sel.outliers <- do.call(.sel.outliers,outliers.args)
        sel <- sel & !sel.outliers
      }
      
      # Select non-NA outlier removed spectra
      channel1 <- channel1[sel]
      channel2 <- channel2[sel]
      i1 <- i1[sel]
      i2 <- i2[sel]
      log.ratios <- log.ratios[sel]
      var.i <- var.i[sel]
    
      ## linear regression estimation
      if (is.a.method("lm")) {
        res.lm <- .calc.lm(channel1,channel2,sign.level.sample,sign.level.rat,ratiodistr) 
        if (is.method("lm")) return (res.lm)
      }

      ## weighted linear regression estimation
      if (is.a.method("weighted lm")) {
        res.wlm <- .calc.weighted.lm(channel1,channel2,var.i,sign.level.sample,sign.level.rat,ratiodistr) 
        if (is.method("weighted lm")) return (res.wlm)
      }      
       
      # First, compute ratios on spectra with intensities for both reporter ions
      lratio.n.var <-
        .calc.weighted.ratio(log.ratios,var.i,variance.function,preweights[sel])
      
      if (is.a.method("ttest")) {
        if (length(log.ratios) < 2) p.value <- 1  else p.value <- t.test(log.ratios)$p.value
        res.ttest <- c(
                       lratio=lratio.n.var['lratio'],variance=NA,n.spectra=length(log.ratios),
            p.value=p.value,is.significant=p.value<sign.level)
        if (is.method("ttest")) return(res.ttest)
      }
#if (method == "wttest" || method == "compare.all") {
#        p.value <- (length(log.ratios) < 2)? 1 : weighted.t.test(log.ratios,w=weights,weighted.ratio)$p.value
#        res.wttest <- c(
#            lratio=weighted.ratio,variance=NA,n.spectra=length(log.ratios),
#            p.value=p.value,is.significant=p.value<sign.level)
#        if (method != "compare.all") return(res.wttest)
#      }
      if (is.a.method("fc")) {
        res.fc <- c(lratio=as.numeric(lratio.n.var['lratio']),variance=NA,
                    n.spectra=length(log.ratios),
                    is.significant=abs(as.numeric(lratio.n.var['lratio']))>log10(fc.threshold))
        ratio.nw <- mean(log.ratios,na.rm=TRUE)
        res.fc.nw <- c(lratio=ratio.nw,variance=NA,n.spectra=length(log.ratios),
            is.significant=abs(ratio.nw)>log10(fc.threshold))
        if (is.method("fc")) return(res.fc)
      }
 
      weighted.ratio <- as.numeric(lratio.n.var['lratio'])
      calc.variance <- as.numeric(lratio.n.var['calc.variance'])
      if (is.a.method("isobar")) {
        # Check for channels where one is NA
        sel.ch1na <- is.na(channel1) & !is.na(channel2)
        sel.ch2na <- is.na(channel2) & !is.na(channel1)

        if (use.na) {
          # TODO: stub
          .calc.w.na()
        } # use.na

        res.isobar <- 
          c(lratio=weighted.ratio, variance=calc.variance,
            n.spectra=length(log.ratios),
            n.na1=sum(sel.ch1na),n.na2=sum(sel.ch2na),
            p.value.rat=calculate.ratio.pvalue(weighted.ratio, calc.variance, ratiodistr),
            p.value.sample=calculate.sample.pvalue(weighted.ratio, ratiodistr),
            is.significant=NA)
        
        if (!is.method("isobar"))
          res.isobar['is.significant.ev'] <- 
            (res.isobar['p.value.sample'] <= sign.level.sample) &&
            calculate.ratio.pvalue(weighted.ratio,lratio.n.var['estimator.variance'],
                                   ratiodistr) <= sign.level.rat

        res.isobar['is.significant'] <- 
          (res.isobar['p.value.sample'] <= sign.level.sample) && 
          (res.isobar['p.value.rat'] <= sign.level.rat)

        if (is.method("isobar")) return(res.isobar)
      }
      if (length(method) == 1 && !is.method("compare.all")) stop(paste("method",method,"not available"))

      add.res <- function(x,name) {
        if (!exists(x)) return(NULL)
        o <- get(x)
        as.numeric(o[name])
      }

      res <- c(lratio=weighted.ratio,
               lratio.lm=add.res("res.lm",'lratio'),
               lratio.wlm=add.res('res.wlm','lratio'),
               unweighted.ratio  = mean(log.ratios,na.rm=TRUE),
               n.spectra=length(log.ratios),
               variance=calc.variance,
               var.ev=as.numeric(lratio.n.var['estimator.variance']),
               var.sv=as.numeric(lratio.n.var['sample.variance']),
               var.lm=add.res('res.lm','stderr')**2,
               var.lm.w=add.res('res.wlm','stderr')**2,
               is.sign.isobar    = add.res('res.isobar','is.significant'),
               is.sign.isobar.ev = add.res('res.isobar','is.significant.ev'),
               is.sign.lm        = add.res('res.lm','is.significant'),
               is.sign.wlm       = add.res('res.wlm','is.significant'),
               is.sign.rat       = add.res('res.isobar','p.value.rat')<sign.level,
               is.sign.sample    = add.res('res.isobar','p.value.sample')<sign.level,
               is.sign.ttest     = add.res('res.ttest','is.significant'),
               is.sign.fc        = add.res('res.fc','is.significant'),
               is.sign.fc.nw     = add.res('res.fc.nw','is.significant'))
                
    }
)

calculate.ratio.pvalue <- function(lratio, variance, ratiodistr = NULL) {
  center.val <-  ifelse(is.null(ratiodistr), 0 , distr::q(ratiodistr)(0.5))
  sapply(seq_along(lratio),function(r.i) 
    pnorm(lratio[r.i],mean=center.val,sd=sqrt(variance[r.i]),lower.tail=lratio[r.i]<center.val)
  )
}

calculate.sample.pvalue <- function(lratio,ratiodistr) {
  if (!is.null(ratiodistr))
    center.val <- distr::q(ratiodistr)(0.5)

  sapply(lratio,function(r) {
    if (is.null(ratiodistr) || is.na(lratio))
      return(NA)
    p(ratiodistr)(r,lower.tail=r<center.val)
  })
}


calculate.mult.sample.pvalue <- function(lratio,ratiodistr,strict.pval,lower.tail,
                                         n.possible.val, n.observed.val) {
  if (is.null(ratiodistr)) 
    return(NA)

  product.p.vals <- prod(distr::p(ratiodistr)(lratio,lower.tail=lower.tail))
  if (strict.pval)
    pval <- getMultUnifPValues(product.p.vals*0.5^(n.possible.val-n.observed.val),n=n.possible.val)
  else
    pval <- getMultUnifPValues(product.p.vals,n=n.observed.val)

  return(pval)
}

adjust.ratio.pvalue <- function(quant.tbl,p.adjust,sign.level.rat,globally=TRUE) {
  if (globally) {
    quant.tbl[,'p.value.rat.adjusted'] <- p.adjust(quant.tbl[,'p.value.rat'], p.adjust)
    quant.tbl[,'is.significant'] <- quant.tbl[,'is.significant'] & quant.tbl[,'p.value.rat.adjusted'] < sign.level.rat
  } else {
    comp.cols <- c("r1","r2","class1","class2")
    comp.cols <- comp.cols[comp.cols %in% colnames(quant.tbl),]
    quant.tbl <- ddply(quant.tbl,comp.cols,function(x) {
      x[,'p.value.rat.adjusted'] <- p.adjust(x[,'p.value.rat'], p.adjust)
      x[,'is.significant'] <- x[,'is.significant'] & x[,'p.value.rat.adjusted'] < sign.level.rat
    })
  }
  quant.tbl
}

correct.peptide.ratios <- function(ibspectra, peptide.quant.tbl, protein.quant.tbl,  
                                   adjust.variance=TRUE, correlation=0, recalculate.pvalue = TRUE) {

  attrs <- attributes(peptide.quant.tbl)
  protein.group <- proteinGroup(ibspectra)
  if (isTRUE(attr(peptide.quant.tbl,'summarize')) || isTRUE(attr(protein.quant.tbl,'summarize')))
    stop("Ratio correction should be done before summarization! Use proteinRatio with argument before.summarize.f!")

  # map from peptides to protein group identifier
  pnp <- peptideNProtein(protein.group)
  pi <- protein.group@peptideInfo
  #pnp <- pnp[pnp[,'protein.g'] %in% reporterProteins(protein.group) & ,]
  all.q.prots <- unique(protein.quant.tbl$ac)
  n.quant <- table(protein.quant.tbl[!is.na(protein.quant.tbl[,'lratio']),'ac'])
  q.protein.acs <- strsplit(all.q.prots,",")
  pep.to.ac <- sapply(unique(peptide.quant.tbl$peptide),function(pep) {
                              pep.prots <- pi[pi[,'peptide'] == pep,'protein']
                              prots.ok <- sapply(q.protein.acs,function(q.prots) any(pep.prots %in% q.prots))
                              if (sum(prots.ok) == 0) {
                                message("could not map peptide ",pep,"!")
                                ac <- "NA"
                              } else if (sum(prots.ok) > 1) {
                                n.quant.acs <- n.quant[all.q.prots[prots.ok]]
                                ac <- names(n.quant.acs)[which.max(n.quant.acs)]
                                message(pep," matches to multiple ACs: ", paste(all.q.prots[prots.ok],collapse=" & ")," [Using ",ac," with ",n.quant.acs[ac],"quantifications]")
                              } else {
                                ac <- all.q.prots[prots.ok]
                              }
                              return(ac)

  })
  peptide.quant.tbl[,'ac'] <- pep.to.ac[peptide.quant.tbl[,'peptide']]

  # merged peptide and protein quant table
  tbl <- merge(peptide.quant.tbl,protein.quant.tbl[,c("ac","r1","r2","lratio","variance","is.significant")],
               by = c("ac","r1","r2"), all.x = TRUE, suffixes=c(".modpep",".prot"))

  tbl[,'lratio'] <- .cn(tbl,'lratio.modpep') - .cn(tbl,'lratio.prot')

  attrs$adjust.variance <- adjust.variance
  if (adjust.variance) {
    attrs$adjust.variance.corralation <- correlation
    cov <- correlation * sqrt(.cn(tbl,'variance.modpep')) * sqrt(.cn(tbl,'variance.prot'))
    tbl[,'variance'] <- .cn(tbl,'variance.modpep') + .cn(tbl,'variance.prot') * 2 * cov
  }

  attrs$recalculate.pvalue <- recalculate.pvalue
  if (recalculate.pvalue) {
    ratiodistr <- attr(peptide.quant.tbl,'ratiodistr')
    tbl[,'p.value.rat'] <- calculate.ratio.pvalue(tbl[,'lratio'], tbl[,'variance'], ratiodistr)
    tbl[,'p.value.sample'] <- calculate.sample.pvalue(tbl[,'lratio'], ratiodistr)
    tbl[,'is.significant'] <- tbl[,'p.value.rat'] <= attr(peptide.quant.tbl,'sign.level.rat') &
                              tbl[,'p.value.sample'] <= attr(peptide.quant.tbl,'sign.level.sample') 
  }
  attrs$names <- colnames(tbl)
  attrs$row.names <- rownames(tbl)
  attributes(tbl) <- attrs

  return(tbl)
}

.calc.lm <- function(channel1,channel2,
                     sign.level.rat,sign.level.sample,ratiodistr) {
  res.lm <- NA
  if (any(!is.na(channel1) & !is.na(channel2))){
    fm <- lm(channel2~channel1+0)
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
  }
  res.lm
}

.calc.weighted.lm <- function(channel1,channel2,var.i,
                              sign.level.rat,sign.level.sample,ratiodistr) {
  res.lm <- NA
  if (any(!is.na(channel1) & !is.na(channel2))){
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
  }
  res.lm
}



.calc.w.na <- function() {
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

#  i1[sel.ch1na] <- max(val,i2[sel.ch1na]-10)
#  i2[sel.ch2na] <- max(val,i1[sel.ch2na]-10)
#  i1.raw[sel.ch1na] <- max(val,i2[sel.ch1na]-10)
#  i2.raw[sel.ch2na] <- max(val,i1[sel.ch2na]-10)
  
  i1[sel.ch1na] <- val
  i2[sel.ch2na] <- val
  i1.raw[sel.ch1na] <- i2.raw[sel.ch1na]
  i2.raw[sel.ch2na] <- i1.raw[sel.ch2na]

  sel <- sel.notna | sel.ch1na | sel.ch2na
#  message("ch1na: ",sum(sel.ch1na),"; ch2na: ",sum(sel.ch2na))

  # calculate final ratio
  lratio.n.var <- 
    .calc.weighted.ratio(sel,i1,i2,var.i,
                         remove.outliers,outliers.coef,outliers.trim,
                         variance.function,preweights)
}

.get.ri <- function(ri,ch) {
  if (is.null(ri))
    return(NULL)
  if (ch == "ALL")
    rowSums(ri,na.rm=TRUE)
  else if (ch == "AVG")
    apply(ri,1,mean,na.rm=TRUE)
  else
    ri[,ch]
}

.sel.outliers <- function(log.ratio,method="boxplot",
                          outliers.coef=1.5,outliers.trim=0.1,
                          weights=NULL) {
  # discussion: http://stats.stackexchange.com/questions/7155/rigorous-definition-of-an-outlier
  sel <- is.na(log.ratio)
  if (method == "boxplot") {
    # Tukey's (1977) method; see ?boxplot.stats and below at 'iqr'.
    # outliers.coef=1.5, typically
    bp <- boxplot.stats(log.ratio,coef=outliers.coef)
    sel | log.ratio<bp$stats[1] | log.ratio>bp$stats[5]
  } else if (method == "iqr") {
    # Tukey's (1977) method
    # applicable to skewed data (makes no distributional assumptions)
    # may not be appropriate for small sample sizes
    #   selects points which are more than outliers.coef times the interquartile range 
    #   above the third quartile or below the first quartile.  
    #  possible outliers: outlier.coef=1.5 (99.3% of the data within range in normal distribution)
    #  probable outliers: outlier.coef=3

    qs <- quantile(log.ratio,c(.25,.75),na.rm=TRUE)
    iqr <- qs[2] - qs[1]
    sel | log.ratio < qs[1]-outliers.coef*iqr | log.ratio > qs[2]+outliers.coef*iqr
  } else if (method == "wtd.iqr") {
    # weighted implementation of iqr method with weighted quantiles
    require(Hmisc)
    qs <- wtd.quantile(log.ratio,weights,c(.25,.75),na.rm=TRUE)
    iqr <- qs[2] - qs[1]
    sel | log.ratio < qs[1]-outliers.coef*iqr | log.ratio > qs[3]+outliers.coef*iqr
  } else if (method == "robust.zscore") {
    # modified z-score proposed by Iglewicz and Hoaglin (1993)
    #  68% have zscores between +/- 1, 95% have zscores between +/- 2
    #  99.7% have zscores between +/- 3
    # outliers.coef=3.5, typically
    sel | 0.6745*abs(log.ratio-mean(log.ratio,na.rm=TRUE))/mad(log.ratio,constant=1,na.rm=TRUE) > outliers.coef
  } else if (method == "zscore") {
    # 99.7% of the data lie within 3 standard deviations of the mean
    # outliers.coef=3, typically
    sel | abs(log.ratio-mean(log.ratio,na.rm=TRUE))/sd(log.ratio,na.rm=TRUE) > outliers.coef

  } else if (method == "trim") {
    sel | log.ratio < quantile(log.ratio,outliers.trim/2,na.rm=TRUE) |
          log.ratio > quantile(log.ratio,1-outliers.trim/2,na.rm=TRUE)
  } else {
    stop("method ",method," not known, use one of [boxplot,iqr,wtd.iqr,zscore,trim]")
  }
}


                           
.calc.weighted.ratio <- function(log.ratio,variance,
                                 variance.function,preweights=NULL) {
  
  variance <- variance
  preweights <- preweights
  
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
      #message("parallel processing: ",isTRUE(getOption('isobar.parallel')))
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
              lratio=weightedMean(peptide.ratios[,'lratio'],weights=peptide.ratios$isum),
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
        res <- ldply(protein,function(individual.protein) {
             if (individual.protein %in% quant.w.grouppeptides) 
               specificity <- c(GROUPSPECIFIC,specificity)

             .call.estimateRatio(individual.protein,"protein",ibspectra,
                                noise.model,channel1,channel2,method=method,...,
                                specificity=specificity)
            },.parallel=isTRUE(getOption('isobar.parallel'))
        )
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
        if (is.data.frame(peptide))
          peptide <- as.matrix(peptide)
        if (is.matrix(peptide)) {
          r <- ldply(seq_len(nrow(peptide)),function(p_i) 
                  .call.estimateRatio(peptide[p_i,,drop=FALSE],"peptide",ibspectra,noise.model,
                                      channel1,channel2,...),.parallel=isTRUE(getOption('isobar.parallel')))
        
        } else {
          r <- ldply(peptide,function(individual.peptide) 
                        .call.estimateRatio(individual.peptide,"peptide",ibspectra,noise.model,
                                            channel1,channel2,...),.parallel=isTRUE(getOption('isobar.parallel')))
        }
      }
      attr(r,"input") <- peptide
      attr(r,"combine") <- combine
      return(r)
}


### Handling NULL protein or peptide argument
setMethod("estimateRatio",
          signature(ibspectra="IBSpectra",noise.model="ANY",
                    channel1="missing",channel2="missing",
                    protein="character",peptide="missing"),
          function(ibspectra,noise.model=NULL,protein,
                   val="lratio",summarize=FALSE,combine=TRUE,...) {
            channels <- reporterTagNames(ibspectra)
            if (combine) {
              res <- matrix(NA,nrow=length(channels),ncol=length(channels),
                            dimnames=list(channels,channels))
              for(channel1 in channels)
                for (channel2 in channels) {
                  rat <- estimateRatio(ibspectra,noise.model=noise.model,
                                       channel1=channel1,channel2=channel2,
                                       protein=protein,...)
                  res[channel1,channel2] <- rat[val]
                }
              return(res)
            } else {
              cmbn <- t(combn(channels,2))
              res <- estimateRatio(ibspectra,noise.model=noise.model,
                                   channel1=cmbn[i,1],channel2=cmbn[i,2],
                                   protein=protein,combine=FALSE,...)

              apply(cmbn,1,function(i) 
                    cbind(r1=cmbn[i,1],r2=cmbn[i,2],ac=rownames(res),res))
            }
            }
)



setMethod("estimateRatio",
          signature(ibspectra="IBSpectra",noise.model="ANY",
                    channel1="missing",channel2="missing",
                    protein="missing",peptide="character"),
          function(ibspectra,noise.model=NULL,peptide,
                   val="lratio",summarize=FALSE,combine=TRUE,...) {
            channels <- reporterTagNames(ibspectra)
            if (combine) {
              res <- matrix(NA,nrow=length(channels),ncol=length(channels),
                            dimnames=list(channels,channels))
              for(channel1 in channels)
                for (channel2 in channels) {
                  rat <- estimateRatio(ibspectra,noise.model=noise.model,
                                       channel1=channel1,channel2=channel2,
                                       peptide=peptide,...)
                  res[channel1,channel2] <- rat[val]
                }
              return(res)
            } else {
              cmbn <- t(combn(channels,2))
              res <- estimateRatio(ibspectra,noise.model=noise.model,
                                   channel1=cmbn[i,1],channel2=cmbn[i,2],
                                   peptide=peptide,combine=FALSE,...)

              apply(cmbn,1,function(i) 
                    cbind(r1=cmbn[i,1],r2=cmbn[i,2],ac=rownames(res),res))
            }
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
              protein="missing",peptide="data.frame"),
    function(ibspectra,noise.model,channel1,channel2,peptide,...) {
      estimateRatioForPeptide(peptide,ibspectra,noise.model,channel1,channel2,...)
    }
)

setMethod("estimateRatio",
    signature(ibspectra="IBSpectra",noise.model="ANY",
              channel1="character",channel2="character",
              protein="NULL",peptide="data.frame"),
    function(ibspectra,noise.model,channel1,channel2,protein=NULL,peptide,...) {
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
                                                                        specificity=specificity,modif=modif,n.sample=n.sample,
                                                                        groupspecific.if.same.ac=groupspecific.if.same.ac,
                                                                        use.precursor.purity=use.precursor.purity,do.warn=do.warn,...))),
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

  i1 <- .get.ri(ri,channel1)
  i2 <- .get.ri(ri,channel2)
  i1.raw <- .get.ri(ri.raw,channel1)
  i2.raw <- .get.ri(ri.raw,channel2)

  # sample data for testing puposes (TP/FP estimation)
  if (!is.null(n.sample)) {
    if (n.sample <= length(i1)) {
      indices <- sample(seq_along(i1),n.sample)
      i1 <- i1[indices]
      i2 <- i2[indices]
      if (!is.null(ri.raw)) {
        i1.raw <- i1.raw[indices]
        i2.raw <- i2.raw[indices]
      }
    }
    else
      i1 <- i2 <- i1.raw <- i2.raw <- numeric()
   }

   estimateRatioNumeric(channel1=i1,channel2=i2,noise.model=noise.model,...,
                        preweights=precursor.purity,
                        channel1.raw=i1.raw,channel2.raw=i2.raw)
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

  # Filter NA channels
  if (!is.null(cl)) {
    x <- x[!is.na(cl)]
    cl <- cl[!is.na(cl)]
  }

  if (!is.null(vs) && !is.character(vs))
    vs <- as.character(vs)

  # Create a combn matrix with all combinations of channels to consider 
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
combn.protein.tbl <- function(cmbn, reverse=FALSE, ...) {
  
  ratios <- do.call(rbind,apply(cmbn,2,function(x) {
    if (reverse)
      if (length(x) == 4)
        x <- x[c(2,1,4,3)]
      else
        x <- rev(x)

    message("ratios ",x[2]," vs ",x[1])
    r <- estimateRatio(channel1=x[1],channel2=x[2],...)
    if (class(r)=="numeric") {
      r <- t(r)
      rownames(r) <- "prot1"
    }

    if (is.matrix(attr(r,"input")))
      df <- data.frame(attr(r,"input"),r,stringsAsFactors=FALSE)
    else
     df <- data.frame(r,stringsAsFactors=FALSE)

    if (!is.null(rownames(r))) 
      df[,'ac']<- rownames(df)
    rownames(df) <- NULL

    df$r1 <- x[1]; df$r2 <- x[2]
    if (length(x) == 4) {
      df$class1 <- x[3]; df$class2 <- x[4]
    }
    return(df)
  }))

  attr(ratios,"arguments") <- list(...)

  attr(ratios,"cmbn") <- cmbn
  attr(ratios,"reverse") <- reverse

  if (all(c("peptide","modif") %in% colnames(ratios))) 
    ratios <- ratios[order(ratios[,'peptide'],ratios[,'modif'],ratios$r1,ratios$r2),]
  else
    ratios <- ratios[order(ratios[,'ac'],ratios$r1,ratios$r2),]
  
  return(ratios)  
}


peptideRatios <- function(ibspectra,...,peptide=peptides(proteinGroup(ibspectra),columns=c("peptide","modif"))) {
  proteinRatios(ibspectra,...,proteins=NULL,peptide=peptide)
}

ratiosReshapeWide <- function(quant.tbl,grouped.cols=TRUE,vs.class=NULL,sep=".",cmbn=NULL,short.names=FALSE) {
  attrs <- attributes(quant.tbl)

  if (!is.null(cmbn)) {
    sel <- paste(quant.tbl$r1,quant.tbl$r2) %in% paste(cmbn[1,],cmbn[2,])
    quant.tbl <- quant.tbl[sel,]
  }
  classes.unique <- "class1" %in% colnames(quant.tbl) &&
                    !any(is.na(quant.tbl$class1)) && !any(is.null(quant.tbl$class1)) && 
                    all(table(unique(quant.tbl[,c("r1","class1")])$class1)==1) &&
                    all(table(unique(quant.tbl[,c("r2","class2")])$class2)==1)

  if (!is.null(vs.class)) {
    if (!any(quant.tbl[,"class1"]==vs.class)) stop("vs.class set to ",vs.class,", but it is not present in quant table")
    if (length(vs.class)==1 && short.names)
      quant.tbl[,'comp']<- paste(quant.tbl$class2)
    else 
      quant.tbl[,'comp']<- paste(quant.tbl$class2,quant.tbl$class1,sep="/")

    quant.tbl <- subset(quant.tbl,class1 %in% vs.class)
  } else {
    if (classes.unique) {
      quant.tbl[,'comp']<- paste(quant.tbl$class2,quant.tbl$class1,sep="/")
    } else {
      quant.tbl[,'comp']<- paste(quant.tbl$r2,quant.tbl$r1,sep="/")
    }
  }
  quant.tbl  <- quant.tbl[,-(c(which(colnames(quant.tbl) %in% c("r1","r2","class1","class2"))))]
  v.names <- c("lratio","lratio.modpep","lratio.prot","variance","variance.modpep","variance.prot","n.spectra",
               "p.value.rat","p.value.rat.adjusted","p.value.sample",
               "is.significant","is.significant.modpep","is.significant.prot","sd","n.na1","n.na2","n.pos","n.neg") # 
  v.names <- v.names[v.names %in% colnames(quant.tbl)]
  timevar <- "comp"
  idvar <- colnames(quant.tbl)[!colnames(quant.tbl) %in% c(timevar,v.names)]
  if ("ac" %in% colnames(quant.tbl))
    quant.tbl[is.na(quant.tbl[,"ac"]),"ac"] <- "NA"

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
  function(ibspectra, noise.model, reporterTagNames=NULL,
           proteins=reporterProteins(proteinGroup(ibspectra)),peptide=NULL,
           cl=classLabels(ibspectra), combn.method="global",symmetry=FALSE,
           summarize=FALSE,summarize.method="mult.pval",
           min.detect=NULL,strict.sample.pval=TRUE,strict.ratio.pval=TRUE,orient.div=0,
           sign.level=0.05,sign.level.rat=sign.level,sign.level.sample=sign.level,
           ratiodistr=NULL,zscore.threshold=NULL,variance.function="maxi",
           combine=FALSE,p.adjust=NULL,reverse=FALSE,
           cmbn=NULL,before.summarize.f=NULL,...) {

  if ((!is.null(proteins) && !is.null(peptide)) ||
      (is.null(proteins) && is.null(peptide)))
    stop("supply either protein or peptides!")      
  
  if (!is.null(p.adjust) && !p.adjust %in% p.adjust.methods)
    stop("p.adjust parameter must be one of '",paste(p.adjust.methods,collapse="','"),"'")
  
  if (is.null(reporterTagNames)) reporterTagNames <- reporterTagNames(ibspectra)
  if (is.null(cl) && is.null(cmbn)) stop("please supply class labels as argument cl or a cmbn matrix")

  if (is.null(cmbn))
    cmbn <- combn.matrix(reporterTagNames,combn.method,cl)

  if (ncol(cmbn) < 1) 
    stop("No possible combination for reporters [",paste(reporterTagNames,sep=","),"]",
         " w/ classes [",paste(cl,sep=","),"] and combn.method ",combn.method," possible.",
         " summarize=",ifelse(summarize,"TRUE","FALSE"))
  
  ratios <- combn.protein.tbl(cmbn,reverse=reverse,
                              ibspectra=ibspectra,noise.model=noise.model,
                              ratiodistr=ratiodistr,
                              protein=proteins,peptide=peptide,
                              sign.level=sign.level,sign.level.rat=sign.level.rat,
                              sign.level.sample=sign.level.sample,
                              variance.function=variance.function,combine=combine,...)

  attributes(ratios) = c(attributes(ratios),list(
                         classLabels=cl,combn.method=combn.method,symmetry=symmetry,summarize=FALSE,
                         sign.level.rat=sign.level.rat,sign.level.sample=sign.level.sample,
                         ratiodistr=ratiodistr,variance.function=variance.function,
                          combine=combine,p.adjust=p.adjust,reverse=reverse))

  if (!is.null(before.summarize.f)) {
    .check.isfunction(before.summarize.f)
    ratios <- before.summarize.f(ibspectra,ratios)
  }


  if (summarize) {
    if (combn.method=="global")
      stop("summarization not meaningful with combn.method='global'. ",
           "Use combn.method='intraclass' or combn.method='interclass' to use ratios in or between classes.")

    ## calculate the number of combinations for each class-combination
    n.combination <- table(cmbn["class1",],cmbn["class2",])
    if (nrow(n.combination)==1 & ncol(n.combination)==1) n.combination <- as.numeric(n.combination)
    if (max(n.combination) < 2)
      stop("Summarize=TRUE makes no sense with only one combination, set summarize to FALSE or class labels differently.")

    if (is.null(min.detect))
      min.detect <- n.combination

    if (!is.null(proteins)) {
      by.column <- "ac"
    } else {
      by.column <- "peptide"
      if (!is.null(dim(peptide)) && "modif" %in% colnames(peptide))
       by.column <- c("peptide","modif") 
    } 
    
    ratios <- summarize.ratios(ratios,by.column,
                           summarize.method,min.detect,n.combination,
                           strict.sample.pval,strict.ratio.pval,orient.div,
                           sign.level,sign.level.rat,sign.level.sample,
                           variance.function=variance.function,
                           ratiodistr=ratiodistr)
  }

  ratios[,'zscore'] <- .calc.zscore(ratios[,'lratio'])
  if (!is.null(zscore.threshold))
    ratios[,'is.significant'] <- ratios[,'p.value.rat'] < sign.level.rat & abs(ratios[,'zscore']) > zscore.threshold


  if (!is.null(p.adjust)) 
    ratios <- adjust.ratio.pvalue(ratios,p.adjust,sign.level.rat)

  if (symmetry) {
    ratios.inv <- ratios
    ratios.inv[,'lratio'] <- -ratios.inv[,'lratio']
    ratios.inv[,c('r1','r2')] <- ratios.inv[,c('r2','r1')]
    if ('class1' %in% colnames(ratios))
      ratios.inv[,c('class1','class2')] <- ratios.inv[,c('class2','class1')]
    ratios <- rbind(ratios,ratios.inv)
  }

  return(ratios)
} # end proteinRatios

summarize.ratios <-
  function(ratios,by.column="ac",summarize.method="mult.pval",min.detect=NULL,n.combination=NULL,
           strict.sample.pval=TRUE,strict.ratio.pval=TRUE,orient.div=0,
           sign.level=0.05,sign.level.rat=sign.level,sign.level.sample=sign.level,
           variance.function="maxi",ratiodistr=NULL) {

    if (!summarize.method  %in% c("mult.pval","mean"))
      stop("Implemented summarize.methods: mult.pval, mean. Please choose one of them instead of ",summarize.method,".")
    if (!all(c("class1","class2") %in% colnames(ratios)))
      stop("'ratios' argument must specify classes w/ columns class1 and class2!")

    classes <- unique(ratios[,c("class1","class2")])
    
    if (is.null(n.combination)) {
      cc <- unique(ratios[,c("r1","r2","class1","class2")])
      n.combination <- table(cc[,"class1"],cc[,"class2"])
      if (nrow(n.combination)==1 & ncol(n.combination)==1) n.combination <- as.numeric(n.combination)
    }
    if (is.null(min.detect))
      min.detect <- n.combination

    if (!all(by.column %in% colnames(ratios))) 
      stop("'by.column' ",paste(by.column,collapse=",")," defined, but no such columns in 'ratios'.")

    mean.r <- ifelse(is.null(ratiodistr),0,distr::q(ratiodistr)(0.5))
    if (summarize.method == "mult.pval") {
      
      result <- ddply(ratios,by.column,function(ratios.subset) {

        ldply(seq_len(nrow(classes)),function(class_i) {

          class1 <- classes[class_i,1]
          class2 <- classes[class_i,2]
          n.combination.c <- ifelse(is.matrix(n.combination),n.combination[class1,class2],n.combination)
          ac.sel <- !is.na(ratios.subset$lratio) & ratios.subset$class1 == class1 & ratios.subset$class2 == class2

          if (!any(ac.sel)) { ## no data for AC and classes
            return(data.frame(ratios.subset[1,by.column],lratio=NA,variance=NA,n.spectra=0,n.pos=0,n.neg=0,
                              p.value.rat=1,p.value.sample=1,is.significant=FALSE,r1=class1,r2=class2,
                              class1=class1,class2=class2,stringsAsFactors=FALSE))
          }
          
          n.pos <- sum(ratios.subset$lratio[ac.sel]>mean.r,na.rm=T)
          n.neg <- sum(ratios.subset$lratio[ac.sel]<mean.r,na.rm=T)
          is.pos <- (n.pos > n.neg && n.neg <= orient.div)
          is.neg <- (n.neg > n.pos && n.pos <= orient.div)

          ## ratio summarization
          good.sel <- ac.sel
          if (!strict.ratio.pval) {
            ## take only positive/negative spectra for ratio summarization
            if (is.pos) 
              good.sel <- good.sel & ratios.subset$lratio>mean.r
            if (is.neg)
              good.sel <- good.sel & ratios.subset$lratio<mean.r
          }

          ac.ratios <- ratios.subset$lratio[good.sel]
          ac.vars <- ratios.subset$variance[good.sel]
          
          sample.var <- weightedVariance(ac.ratios,weights=1/ac.vars)
          estim.var <- 1/sum(1/ac.vars)
          
          variance <- switch(variance.function,
                             maxi = max(estim.var,sample.var,na.rm=T),
                             ev = estim.var,
                             wsv = sample.var
                             )  
          lratio <- weightedMean(ac.ratios,weights=1/ac.vars)
          
          ## p-value computation
          p.value.rat <- pnorm(lratio,mean=mean.r,sd=sqrt(variance),lower.tail=lratio<mean.r)

          p.value.sample <- calculate.mult.sample.pvalue(
            ac.ratios, ratiodistr, strict.pval=strict.sample.pval, lower.tail=n.neg>n.pos,
            n.possible.val = n.combination.c, n.observed.val = sum(ac.sel))

          min.detect.c <- ifelse(is.matrix(min.detect),min.detect[class1,class2],min.detect)
          ## significance
          is.significant <- (p.value.sample <= sign.level.sample) &&
            (p.value.rat <= sign.level.rat) &&
            (sum(ac.sel) >= min.detect.c) &&
            (is.pos | is.neg)

          return(data.frame(ratios.subset[1,by.column],lratio=lratio,variance=variance,
                            n.spectra=min(ratios$n.spectra[ac.sel]),n.pos=n.pos,n.neg=n.neg,
                            p.value.rat=p.value.rat,p.value.sample=p.value.sample,
                            is.significant=is.significant,r1=class1,r2=class2,
                            class1=class1,class2=class2,stringsAsFactors=FALSE))
        })
      },.parallel=isTRUE(getOption('isobar.parallel')))
      
      if (is.null(result)) stop("Error summarizing.")

      attributes(result) = c(attributes(result),
                             list(summarize=TRUE,
                                  by.column=by.column,
                                  min.detect=min.detect,
                                  orient.div=orient.div))

      return(result)
      
    } else if (summarize.method=="mean") {
      stop("summarize method mean not anymore available.")
#      do.call(rbind,lapply(unique(ratios$ac),function(ac) {
#        ac.sel <- (ratios$ac==ac) & !is.na(ratios$lratio)
#        if (!any(ac.sel)) {
#            return(data.frame(ac=ac,lratio=NA,stringsAsFactors=FALSE))
#        }
#        return(data.frame(ac=ac,lratio=mean(ratios$lratio[ac.sel]),stringsAsFactors=FALSE))
#      }))

    } else {
      stop (paste("summarize method",summarize.method,"not implemented."))
    }
  }


.calc.zscore <- function(lratio) {
  s.median <- median(lratio,na.rm=TRUE)
  s.mad <- mad(lratio,na.rm=TRUE)
  (lratio-s.median)/s.mad
}
