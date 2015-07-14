


calcProbXDiffNormals <- function(X,mu_Y,sd_Y,...,alternative=c("greater","less","two-sided"),progress=FALSE) {
  # calculates a vector of either P(X>=Y) or P(Y>=X)
  # (depending on 'alternative' option),
  # where Y ~ N(mu_Y,sd_Y)
  alternative <- match.arg(alternative)
  if (is.null(X)) {
	return(rep(NA,seq_along(mu_Y)))
  }

  require(distr)
  if (!is(X,"Distribution")) stop("X has to be of class Distribution (is ",class(X),")")

  if (length(mu_Y) != length(sd_Y)) 
    stop("mu_Y and sd_Y should have equal length")
 
  if (isTRUE(progress))
    pb <- txtProgressBar(min=1,max=length(mu_Y))

  X_mode <- distr::q(X)(0.5)
  p.values <- sapply(seq_along(mu_Y),function(i) {
    if (isTRUE(progress))
      setTxtProgressBar(pb,i)
    if (is.na(mu_Y[i]) || is.na(sd_Y[i])) {
      return(NA)
    }
    else {
      Y <- Norm(mu_Y[i],sd_Y[i])
      # FIXME P-value not adjusted for two-sided case
      if ( alternative == 'greater'
        || (alternative == 'two-sided' && X_mode < mu_Y[i])
      ){
        calcProbXGreaterThanY(X,Y,...)
      } else {
        calcProbXGreaterThanY(Y,X,...)
      }
    }
  })
  if (isTRUE(progress))
    close(pb)

  p.values
}

calcProbXGreaterThanY <- function(X, Y, rel.tol=.Machine$double.eps^0.25, subdivisions=100L)
{
  # numerically calculates P(X>=Y)
  require(distr)
  if (!is(X,"Distribution")) stop("X has to be of class Distribution")
  if (!is(Y,"Distribution")) stop("Y has to be of class Distribution")
#  return ( distr::p(Y-X)(0) ) # the easiest way, but it does not seem to give accurate results for Norm+Tlsd

  # FIXME use the distribution info to deduce integration boundaries
  return ( integrate( function(t) distr::d(X)(t)*distr::p(Y)(t),
                      lower=-Inf, upper=Inf, subdivisions = subdivisions,
                      rel.tol = rel.tol )$value )
}

.calcProbXGreaterThanY.orig <- function(X,Y,min.q=10^-6, subdivisions=100L) {
  # numerically calculates P(X>=Y)
  require(distr)
  if (!is(X,"Distribution")) stop("X has to be of class Distribution")
  if (!is(Y,"Distribution")) stop("Y has to be of class Distribution")

  steps.seq <- seq(from=min.q,to=1-min.q,length.out=subdivisions)
  steps <- sort(unique(c(q(X)(steps.seq),q(Y)(steps.seq))))
  nsubdivisions <- length(steps)
           
  dens.X <- d(X)(steps)
  dens.Y <- d(Y)(steps)

  dens.sum <- sum(dens.X*dens.Y)*.5 # P( X = Y )
  dens.total <- sum(dens.X)*sum(dens.Y)
  if (isTRUE(all.equal(dens.total,0))) return(0)

  for (i in 1:nsubdivisions) {
    dens.X <- dens.X[-1]
    dens.Y <- dens.Y[-length(dens.Y)]
    dens.sum <- dens.sum + sum(dens.X*dens.Y)
  }

  dens.sum/dens.total
}

calcCumulativeProbXGreaterThanY <- function(Xs, mu_Ys, sd_Ys,
                                            alternative=c("greater","less","two-sided"),
                                            rel.tol=.Machine$double.eps^0.25, subdivisions=100L)
{
  # numerically calculates P(X_cum>=Y_cum),
  # where X_cum variable has the PDF of joint X distribution induced to the diagonal (x_1=x_2=..x_n),
  # and Y_cum is the joint normal distribution induced to the diagonal as well
  # it's assumed that the mode of all Xs is very close to zero

  alternative <- match.arg(alternative)

  # cumulative Y distribution is Gaussian with explicitly calculated parameters 
  require(distr)
  inv_vars <- sd_Ys^(-2)
  var_Cum <- 1 / sum( inv_vars )
  mu_cum_Y <- sum( as.numeric( mu_Ys %*% inv_vars ) * var_Cum )
  # flip Norm around zero (the mode of X) for 'two-sided' or 'less' tests
  if ( alternative == 'less' || (alternative == 'two-sided' && mu_cum_Y < 0 ) ) {
    mu_cum_Y <- -mu_cum_Y
  }
  cum_Y <- distr::Norm( mu_cum_Y, sqrt( var_Cum ) )

  d_Xs <- lapply( Xs, distr::d )
  # multiplication of all X pdfs
  # note that pdf_cum_X is not a real PDF,
  # because it's not guaranteed to integrate to 1.0
  pdf_cum_X <- function( t ) {
    # apply each d_X to arguments
    pdf_Xs <- lapply( d_Xs, function(d_X) d_X(t) )
    # multiply all PDFs for each element in t
    sapply( seq_along(t), function(ix) prod(sapply(pdf_Xs, function(pdf_X) pdf_X[[ix]])) )
  } 

  pdf_cum_XmY <- function( t ) pdf_cum_X(t) * distr::p( cum_Y )(t)

  return ( integrate( pdf_cum_XmY,
                      lower=-Inf, upper=Inf, subdivisions = subdivisions,
                      rel.tol = rel.tol )$value /
           # correct the scaling of cum_X
           integrate( pdf_cum_X,
                      lower=-Inf, upper=Inf, subdivisions = subdivisions,
                      rel.tol = rel.tol )$value )
}

distrprint <- function(X,round.digits=5) {
  paste0(class(X)," [",
         paste(sapply(setdiff(slotNames(param(X)),'name'),
                      function(x) paste(x,round(slot(param(X),x),round.digits),
                                        sep=" ")),collapse="; "),
         "]")
}

twodistr.plot <- function(X,Y,n.steps=1000,min.q=10^-3) {
  require(distr)
  require(ggplot2)

  steps.seq <- seq(from=0,to=1,length.out=n.steps)
  steps <- sort(unique(c(q(X)(steps.seq),q(Y)(steps.seq))))
  ggplot(rbind(data.frame(x=steps,Distribution=paste0("X ~ ",distrprint(X)),density=d(X)(steps)),
               data.frame(x=steps,Distribution=paste0("Y ~ ",distrprint(Y)),density=d(Y)(steps)))) + 
    geom_line(aes_string(x="x",y="density",color="Distribution")) +
    ggtitle(paste0("P(X>=Y) = ",round(calcProbXGreaterThanY(X,Y),7)))
 
}


