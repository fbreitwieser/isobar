shrink.ratios2 <- function(pr,do.plot=0,nrand=10000,use.n=FALSE) {
  ## using the formulas developed in http://www.cs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf

  require(MASS)
  sel <- !is.na(pr$variance) 

  ## fit gamma distribution on the inverse of the sample variance (used as prior)
  x <- pr$variance[!is.na(pr$variance)&pr$n.spectra>2]
  (fit.gamma <- fitdistr(1/x,"gamma"))
  
  var.ks <- ks.test(1/x,pgamma,
          shape=fit.gamma$estimate['shape'],
          rate=fit.gamma$estimate['rate'])
  message(sprintf("fitting inverse gamma on variance (two-sided ks D=%.2f, p-value=%.3f)",
          var.ks$statistic,var.ks$p.value))

  ## fit normal distribution on sample mu (used as prior)
  prior.mu <- median(pr$lratio[sel])
  prior.var <- mad(pr$lratio[sel])**2

  mu.ks <- ks.test(pr$lratio[sel],"dnorm",
                   prior.mu,sqrt(prior.var))

  message(sprintf("fitting normal on mean (two-sided ks D=%.2f, p-value=%.3f)",
          mu.ks$statistic,mu.ks$p.value))


  ## calculate posterior variance
  prior.alpha <- fit.gamma$estimate[1]
  prior.beta <- fit.gamma$estimate[2]

  sample.n <- pr$n.spectra
  sample.var <- pr$variance

  # posterior alpha = prior.alpha + n/2
  # posterior beta = prior.beta + 1/2*sum(x_i-\bar{x})^2 + 
  if (use.n) {
    post.shape <- prior.alpha + sample.n/2
    post.scale <- prior.beta+1/2*sample.var*(pr$n.spectra-1) 
  } else {
    post.shape <- prior.alpha+sample.n/2 - 1/2
    post.scale <- prior.beta+1/2*sample.var*(pr$n.spectra-1) 
  }

  post.shape <- prior.alpha+sample.n/2 - 1/2
  post.scale <- prior.beta+1/2*sample.var*(pr$n.spectra-1) 

  post.var <- sapply(seq_len(nrand),function(i)
       1/rgamma(sum(sel),post.shape[sel],post.scale[sel]))

  if (do.plot > 0) {
    par(mfrow=c(do.plot,2))
    if (do.plot > 1) {
    plot(density(x),"variance prior")
    seqi <- seq(from=min(x),to=max(x),length.out=1000)
    lines(seqi,dgamma(1/seqi,
                      shape=fit.gamma$estimate['shape'],
                      scale=fit.gamma$estimate['rate']),
            col="red",type="l")
    }
    plot(pr$variance[sel],rowMeans(post.var),
         cex=sqrt(pr$n.spectra[sel])/10,
         log="xy",col="#000000A0",main="sample vs posterior variance",
         xlab="sample variance", ylab="posterior variance")
    abline(0,1,col="red")
  }
 
  ## calculate posterior mu
  post.mu <- sapply(seq_len(nrand),function(i) {
      mu.post.mu <- (prior.mu+pr$n.spectra[sel]*pr$lratio[sel])/
                     (pr$n.spectra[sel] + 1)
      mu.post.var <- post.var[,i]/(pr$n.spectra[sel]+1)
      rnorm(sum(sel),mu.post.mu,sqrt(mu.post.var))
  })


  if (do.plot) {
    if (do.plot > 1) {
    plot(density(pr$lratio[sel]),main="mean prior")
    seqi <- seq(from=-max(abs(pr$lratio[sel])),
                            to=max(abs(pr$lratio[sel])),length.out=10000)
    lines(seqi,dnorm(seqi,prior.mu,sqrt(prior.var)),
                col="red",type="l")
    }
    plot(pr$lratio[sel],rowMeans(post.mu),
         cex=sqrt(pr$n.spectra[sel])/10,
         log="xy",col="#000000A0",main="sample vs posterior mean",
         xlab="sample mean", ylab="posterior mean")
    abline(0,1,col="red")
  }
 
  pr$lratio[sel] <- rowMeans(post.mu)
  pr$variance[sel] <- rowMeans(post.var)

  pr
}


fit.normal <- function(x,robust=FALSE) {
  x <- x[!is.na(x)]

  if (robust) {
    mu <- median(x)
    variance <- mad(x)**2
  } else {
    mu <- mean(x)
    variance <- sd(x)**2
  }

  mu.ks <- ks.test(x,"dnorm",
                   mu,sqrt(variance))
  message(sprintf("fitting normal on mean: mean=%.2f, variance=%.2f (two-sided ks D=%.2f, p-value=%.3f)",
                  mu,variance,mu.ks$statistic,mu.ks$p.value))

  c(mean=mu,variance=variance)
}

fit.gamma <- function(x) {
  require(MASS)

  x <- x[!is.na(x)]
  (fit.gamma <- fitdistr(x,"gamma"))
  
  var.ks <- ks.test(x,pgamma,
          shape=fit.gamma$estimate['shape'],
          rate=fit.gamma$estimate['rate'])
  message(sprintf("fitting inverse gamma on variance: alpha=%.2f, beta=%.2f (two-sided ks D=%.2f, p-value=%.3f)",
                  fit.gamma$estimate[1],fit.gamma$estimate[2],var.ks$statistic,var.ks$p.value))


  c(alpha=as.numeric(fit.gamma$estimate[1]),
    beta=as.numeric(fit.gamma$estimate[2]))
}

fit.normal.gamma <- function(pr) {
  prior.mu <- fit.normal(pr$lratio[pr$n.spectra>3])
  prior.phisq <- fit.gamma(1/pr$variance[pr$n.spectra>3])

  c(prior.alpha=as.numeric(prior.phisq['alpha']),prior.beta=as.numeric(prior.phisq['beta']),
    prior.mu=as.numeric(prior.mu['mean']),prior.var=as.numeric(prior.mu['variance']))
}

shrink.ratios.gibbs <- function(pr,do.plot=0,S=1000,use.n=FALSE) {
  ## using the formulas developed in UW S509

  ## fit priors
  prior.mu <- fit.normal(pr$lratio)
  prior.phisq <- fit.gamma(1/pr$variance[pr$n.spectra>3])

  samples <- gibbs.sample.normalgamma(data.x=pr$lratio,data.var=pr$variance,n=pr$n.spectra,
                                      prior.alpha=prior.phisq['alpha'],prior.beta=prior.phisq['beta'],
                                      prior.mu=prior.mu['mean'],prior.var=prior.mu['variance'],S=S)
 
  pr$lratio <- colMeans(samples$mean)
  pr$variance <- colMeans(1/samples$phisq)

  pr
}

gibbs.sample.normalgamma <- function(data.x,data.var,n,prior.alpha,prior.beta,prior.mu,prior.var,S=1000,do.plot=0) {
  if (is.null(prior.alpha) || is.null(prior.beta) || is.null(prior.mu) || is.null(prior.var))
    stop("Prior parameters may not be NULL")

  dim.data <- length(data.x)

  samples.mean <- matrix(nrow=S,ncol=dim.data)
  samples.phisq <- matrix(nrow=S,ncol=dim.data)

  ### Starting values
  samples.mean[1,] <- data.x
  samples.phisq[1,] <- 1/data.var

  for (s in 2:S) {
    ## Generate a new posterior value of mu from f(mu | prev value of phi.sq, data)
    prev.phisq <- samples.phisq[s-1,]
    prev.mu <- samples.mean[s-1,]

    post.shape <- prior.alpha + n/2
    #post.scale <- prior.beta + 1/2*sum((pr$lratio-prev.mu)^2)  ### HERE should be the actual data such that it works as Gibbs sampler; ie \sum_{i=1}^n x_i - prev.mu for EACH protein!!
    post.scale <- prior.beta+1/2*data.var*(n-1)

    new.phisq <- rgamma(ncol(samples.phisq),post.shape,post.scale)

    ## Generate a new posterior value of phi.sq from f(phisq | prev value of mu, data)
    mu.post.phisq <- 1/prior.var + n*new.phisq
    mu.post.mu <- 1/mu.post.phisq * (prior.mu/prior.var + n*data.x*new.phisq)

    new.mu <- rnorm(ncol(samples.mean),mu.post.mu,1/sqrt(mu.post.phisq))
  
    samples.mean[s,] <- new.mu
    samples.phisq[s,] <- new.phisq
  }
  
  if (do.plot > 0) {
    par(mfrow=c(do.plot,2))

    ## variance
    if (do.plot > 1) {
    plot(density(data.var),main="variance prior")
    seqi <- seq(from=min(data.var),to=max(data.var),length.out=1000)
    lines(seqi,dgamma(1/seqi,
                      shape=fit.gamma$estimate['shape'],
                      scale=fit.gamma$estimate['rate']),
            col="red",type="l")
    }
    plot(data.var,1/colMeans(samples.phisq),
         cex=sqrt(n)/10,
         log="xy",col="#000000A0",main="sample vs posterior variance",
         xlab="sample variance", ylab="posterior variance")
    abline(0,1,col="red")

    ## means
    if (do.plot > 1) {
    plot(density(data.x),main="mean prior")
    seqi <- seq(from=-max(abs(data.x)),
                            to=max(abs(data.x)),length.out=10000)
    lines(seqi,dnorm(seqi,prior.mu,sqrt(prior.var)),
                col="red",type="l")
    }
    plot(data.x,colMeans(samples.mean),
         cex=sqrt(n)/10,
         log="xy",col="#000000A0",main="sample vs posterior mean",
         xlab="sample mean", ylab="posterior mean")
    abline(0,1,col="red")
  }

  return(list(mean=samples.mean,phisq=samples.phisq))
}

gibbs.sample.normalgamma.singleprotein <- function(data.x,data.var,prior.alpha,prior.beta,prior.mu,prior.var,S=1000,just.means=TRUE, burn.in = 0.1, prune=0.1,weighted=FALSE) {    
  if (length(data.x)==0 || length(data.var) == 0)
    return(c(mean=NA,phisq=NA))

  if (is.null(prior.alpha) || is.null(prior.beta) || is.null(prior.mu) || is.null(prior.var))
    stop("Prior parameters may not be NULL")

  if (is.null(S) || is.na(S))
    S <- 1000

  dim.data <- 1
  n <- length(data.x)

  samples.mean <- matrix(nrow=S,ncol=dim.data)
  samples.phisq <- matrix(nrow=S,ncol=dim.data)

  data.mean <- weightedMean(data.x,1/data.var)

  ### Starting values
  samples.mean[1,] <- data.mean
  if (length(data.var) == 1)
      samples.phisq[1,] <- 1/data.var
  else
      samples.phisq[1,] <- 1/weightedVariance(data.x,1/data.var,mean.estimate=data.mean)[1]

  for (s in 2:S) {
    ## Generate a new posterior value of mu from f(mu | prev value of phi.sq, data)
    prev.phisq <- samples.phisq[s-1,]
    prev.mu <- samples.mean[s-1,]

    post.shape <- prior.alpha + n/2
    if (weighted) 
      post.scale <- prior.beta + 1/2*n*biasedWeightedVariance(data.x,1/var.i,prev.mu)
    else
      post.scale <- prior.beta + 1/2*sum((data.x-prev.mu)^2)

    new.phisq <- rgamma(1,post.shape,post.scale)

    ## Generate a new posterior value of phi.sq from f(phisq | prev value of mu, data)
    mu.post.phisq <- 1/prior.var + n*new.phisq
    mu.post.mu <- 1/mu.post.phisq * (prior.mu/prior.var + n*data.mean*new.phisq)

    new.mu <- rnorm(1,mu.post.mu,1/sqrt(mu.post.phisq))
  
    samples.mean[s,] <- new.mu
    samples.phisq[s,] <- new.phisq
  }
 
  if (just.means) {
    return(c(mean=mean(samples.mean[,1]),phisq=mean(samples.phisq[,1]),mean.var=var(samples.mean[,1]),
             post.scale=as.numeric(post.scale),post.shape=as.numeric(post.shape)))
  } else {
    return(list(mean=samples.mean,phisq=samples.phisq))
  }
}



shrink.ratios <- function(pr,do.plot=FALSE,nrand=10000) {
  require(MASS)
  require(LaplacesDemon)

  ## fit inverse gamma distribution on sample variance (used as prior)
  x <- pr$variance[!is.na(pr$variance)&pr$n.spectra>2]
  (fit.gamma <- fitdistr(1/x,"gamma"))
  
  var.ks <- ks.test(1/x,pgamma,
          shape=fit.gamma$estimate['shape'],
          rate=fit.gamma$estimate['rate'])
  message("fitting inverse gamma on variance (two-sided ks D=%.2f, p-value=%.3f)",
          var.ks$statistic,var.ks$p.value)

  ## calculate posterior variance
  prior.alpha <- fit.gamma$estimate[1]
  prior.beta <- fit.gamma$estimate[2]
  sample.n <- pr$n.spectra
  sample.var <- pr$variance

  post.shape <- prior.alpha+sample.n/2-1/2
  post.scale <- prior.beta+1/2*sample.var*pr$n.spectra
  sel <- !is.na(pr$variance) 
  post.var <- sapply(seq_len(nrand),function(i)
       rinvgamma(sum(sel),post.shape[sel],post.scale[sel]))

  if (do.plot) {
    par(mfrow=c(2,2))
    plot(density(x),"variance prior")
    seqi <- seq(from=0,to=max(x),length.out=1000)
    lines(seqi,LaplacesDemon::dinvgamma(seqi,
                      shape=fit.gamma$estimate['shape'],
                      scale=fit.gamma$estimate['rate']),
            col="red",type="l")
    plot(pr$variance[sel],rowMeans(post.var),
         cex=sqrt(pr$n.spectra[sel])/10,
         log="xy",col="#000000A0",main="sample vs posterior variance")
    abline(0,1,col="red")
  }
  
  ## fit normal distribution on sample mu (used as prior)
  prior.mu <- median(pr$lratio[sel])
  prior.var <- mad(pr$lratio[sel])**2

  mu.ks <- ks.test(pr$lratio[sel],"dnorm",
                   prior.mu,sqrt(prior.var))

  message("fitting normal on mean (two-sided ks D=%.2f, p-value=%.3f)",
          mu.ks$statistic,mu.ks$p.value)

  ## calculate posterior mu
  post.mu <- sapply(seq_len(nrand),function(i) {
      mu.post.mu <- (prior.mu+pr$n.spectra[sel]*pr$lratio[sel])/
                     (pr$n.spectra[sel] + 1)
      mu.post.var <- post.var[,i]/(pr$n.spectra[sel]+1)
      rnorm(sum(sel),mu.post.mu,sqrt(mu.post.var))
  })


  if (do.plot) {
    plot(density(pr$lratio[sel]),main="mean prior")
    seqi <- seq(from=-max(abs(pr$lratio[sel])),
                            to=max(abs(pr$lratio[sel])),length.out=10000)
    lines(seqi,dnorm(seqi,prior.mu,sqrt(prior.var)),
                col="red",type="l")

    plot(pr$lratio[sel],rowMeans(post.mu),
         cex=sqrt(pr$n.spectra[sel])/10,
         log="xy",col="#000000A0",main="sample vs posterior mean")
    abline(0,1,col="red")
  }
 
  pr$lratio[sel] <- rowMeans(post.mu)
  pr$variance[sel] <- rowMeans(post.var)

  pr
}
