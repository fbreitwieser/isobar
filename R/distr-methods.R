


calcTwoSidedNormalProb <- function(X,mu_Y,sd_Y,...) {
  # calculates an estimate of min(P(X>Y),P(Y>X))
  #  Y ~ N(mu_Y,sd_Y)
  require(distr)

  if (length(mu_Y) != length(sd_Y)) 
    stop("mu_Y and sd_Y should have equal length")
 
  p.values <- sapply(seq_along(mu_Y),function(i) {
    if(is.na(mu_Y[i]) || is.na(sd_Y[i])) return(NA)
    calcProbXGreaterThanY(X,Norm(mu_Y[i],sd_Y[i]),...)
  })

  ifelse(p.values>0.5,1-p.values,p.values)
}

calcProbXGreaterThanY <- function(X,Y,n.steps=1000) {
  # calculates an estimate of P(X>Y) + 0.5 P(X=Y)
  #  a value of 0.5 represents complete overlap
  require(distr)
  if (!is(X,"Distribution")) stop("X has to be of class Distribution")
  if (!is(Y,"Distribution")) stop("Y has to be of class Distribution")

  steps.seq <- seq(from=0,to=1,length.out=n.steps)
  steps <- sort(unique(c(q(X)(steps.seq),q(Y)(steps.seq))))
  nn.steps <- length(steps)
           
  dens.X <- d(X)(steps)
  dens.Y <- d(Y)(steps)

  dens.sum <- sum(dens.X*dens.Y)*.5 # P( X = Y )
  dens.total <- sum(dens.X)*sum(dens.Y)
  if (isTRUE(all.equal(dens.total,0))) return(0)

  for (i in 1:nn.steps) {
    dens.X <- dens.X[-1]
    dens.Y <- dens.Y[-length(dens.Y)]
    dens.sum <- dens.sum + sum(dens.X*dens.Y)
  }

  dens.sum/dens.total
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
    geom_line(aes(x=x,y=density,color=Distribution)) +
    ggtitle(paste0(".5 P(X=Y) + P(X>Y) = ",round(calcProbXGreaterThanY(X,Y,n.steps=n.steps),7)))
 
}


