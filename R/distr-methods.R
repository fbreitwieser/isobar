


calcProbXDiffNormals <- function(X,mu_Y,sd_Y,...,alternative="greater",progress=FALSE) {
  # calculates an estimate of min(P(X>Y),P(Y>X))
  #  Y ~ N(mu_Y,sd_Y)
  require(distr)

  if (length(mu_Y) != length(sd_Y)) 
    stop("mu_Y and sd_Y should have equal length")
 
  if (isTRUE(progress))
    pb <- txtProgressBar(min=1,max=length(mu_Y))

  p.values <- sapply(seq_along(mu_Y),function(i) {
    if (isTRUE(progress))
      setTxtProgressBar(pb,i)
    if(is.na(mu_Y[i]) || is.na(sd_Y[i])) return(NA)
    calcProbXGreaterThanY(X,Norm(mu_Y[i],sd_Y[i]),...)
  })
  if (isTRUE(progress))
    close(pb)

  switch(alternative,
         greater=p.values,
         less=1-p.values,
         "two-sided"=ifelse(p.values>0.5,1-p.values,p.values),
         stop("don't know alternative ",alternative))
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


