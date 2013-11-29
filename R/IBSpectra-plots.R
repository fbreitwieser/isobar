### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### plotting methods.
###

setGeneric("reporterMassPrecision", function(x,plot,...)
           standardGeneric("reporterMassPrecision"))
setGeneric("maplot",function(x,channel1,channel2,...) standardGeneric("maplot"))
setGeneric("plotRatio",function(x,channel1,channel2,protein,...)
           standardGeneric("plotRatio"))
setGeneric("protGgdata",function(x,relative.to,protein,...)
           standardGeneric("protGgdata"))
setGeneric("raplot",function(x,...) standardGeneric("raplot"))

setMethod("reporterMassPrecision",
          signature=c(x="IBSpectra",plot="missing"),
          function(x) reporterMassPrecision(x,TRUE))


setGeneric("reporterIntensityPlot",function(x) standardGeneric("reporterIntensityPlot"))

# Calculates and displays the deviation from the 'true' tag mass 
# - as specified in the IBSpectra object - of each channel.
setMethod("reporterMassPrecision", signature=c(x="IBSpectra",plot="logical"),
    function(x,plot=TRUE,nrow=NULL) {
      masses <- reporterMasses(x)
      melt.masses <- data.frame(reporter=rep(colnames(masses),each=nrow(masses)),
                                reporter.tag.mass=rep(reporterTagMasses(x),each=nrow(masses)),
                                observed.moz=as.numeric(masses),stringsAsFactors=FALSE)
      melt.masses$mass.difference <-
          melt.masses$observed.moz - rep(reporterTagMasses(x),each=nrow(masses))
      if (plot) {
        require(ggplot2)
        melt.masses$reporter <-
          factor(melt.masses$reporter,
                 levels=reporterTagNames(x),
                 labels=sprintf("%s: m/z %.2f",
                   reporterTagNames(x),reporterTagMasses(x)))
       
       if (compareVersion(packageDescription("ggplot2")$Version,"0.9.1") <= 0) {
        opts.f <- opts; text.f <- theme_text;
       } else {
        opts.f <- theme; text.f <- element_text;
       }
       if (is.null(nrow))
         nrow <- ceiling(ncol(masses)/5)

       ggplot(melt.masses,aes(x=mass.difference)) + geom_vline(xintercept=0,alpha=0.8) +
          geom_histogram(fill="white",aes(colour=factor(reporter)),alpha=0.8,
                         binwidth=1/20*(max(melt.masses$mass.difference,na.rm=TRUE)-min(melt.masses$mass.difference,na.rm=TRUE))) + 
            facet_wrap(~reporter,scales="fixed",nrow=nrow) + 
            theme_bw(base_size=10) + xlab("mass difference theoretical vs observed reporter tag mass") +
              opts.f(legend.position="none",
                   axis.text.x = text.f(angle=330,hjust=0,vjust=1,colour="grey50",size=7))
      } else {
        res <- ddply(melt.masses,'reporter',function(x) {
          c('true reporter mass'=x$reporter.tag.mass[1],
            'number of spectra'=sum(!is.na(x$observed.moz)),
            'mean of observed m/z'=mean(x$observed.moz,na.rm=TRUE),
            'sd of obseved m/z'=sd(x$observed.moz,na.rm=TRUE),
            'mean of mass difference * 10^3'=mean(x$mass.difference*10^3,na.rm=TRUE),
            'sd of mass difference * 10^3'=sd(x$mass.difference,na.rm=TRUE))
        })
        rownames(res) <- res$reporter
        res$reporter <- NULL
        res
      }
    }
)

.reshapeLong <- function(x,key="key",value="value") {
  res <- data.frame(key=rep(colnames(x),each=nrow(x)),
                    value=as.numeric(x),
                    stringsAsFactors=FALSE)
  colnames(res) <- c(key,value)
  res
}

setMethod("reporterIntensityPlot",
          signature=c(x="IBSpectra"),function(x) {
            require(ggplot2)
            intensities <- reporterIntensities(x)
            intensities.nn <- reporterData(x,element="ions_not_normalized") # null if not normalized
            
            melt.intensities <- data.frame(tag=rep(colnames(intensities),each=nrow(intensities)),
                       normalized=ifelse(is.null(intensities.nn),"no","2. after normalization"),
                       intensity=as.numeric(intensities),
                       stringsAsFactors=FALSE)

            if (!is.null(intensities.nn)) {
            	melt.intensities <- rbind(melt.intensities,
              	  data.frame(tag=rep(colnames(intensities.nn),each=nrow(intensities.nn)),
                             normalized="1. before normalization",
                             intensity=as.numeric(intensities.nn),
                             stringsAsFactors=FALSE))
            }
            
            ggplot(melt.intensities,aes(x=tag,y=intensity)) +
              geom_boxplot(aes(color=factor(normalized)),size=0.5,alpha=0.6,
                           outlier.size=0.5,position=position_dodge(width=0.25)) + 
              xlab("isobaric reporter tag") +
              scale_y_log10() + theme_bw(base_size=10) + scale_color_hue("") 
          }
)

# Ratio-Absolute intensity plot - will be deprecated by maplot
setMethod("raplot",signature(x="IBSpectra"),
    function(x,...) {
      ions <- reporterIntensities(x,na.rm=T)
      par(mfrow=c(ncol(ions)/2,ncol(ions)/2),mar=c(3,2,1,1))
      sums = rowSums(ions)
      for (i in colnames(ions)) {
        plot(sums,ions[,i]/sums*4,log="xy",...)
      }
    }
)

# plots ratio of one protein
setMethod("plotRatio",
    signature(x="IBSpectra",channel1="character",channel2="character",protein="character"),
    function(x,channel1,channel2,protein,noise.model=NULL,...) {
      ions <- reporterIntensities(x,na.rm=TRUE)
      R=log10(ions[,channel1])
      G=log10(ions[,channel2])
      M <- R - G
      A <- 0.5*(R + G)
      
      M <- M[order(A)]
      A <- sort(A)
      plot(A,M,pch=pch,...)
      #abline(h=0,col="red")
      if (!missing(noise.model) & !is.null(noise.model)) {
        if (is(noise.model,"NoiseModel"))
          noise.model <- c(noise.model)
        for (i in seq_along(noise.model)) {          
          lines(A,1.96*sqrt(variance(noise.model[[i]],A)),col=noise.model.col[i])
          lines(A,-1.96*sqrt(variance(noise.model[[i]],A)),col=noise.model.col[i])
        }
      }
      
      unspecific.spectra.sel    <- names(A) %in% spectrumSel(x,protein=protein,specificity=UNSPECIFIC)
      groupspecific.spectra.sel <- names(A) %in% spectrumSel(x,protein=protein,specificity=GROUPSPECIFIC)
      specific.spectra.sel      <- names(A) %in% spectrumSel(x,protein=protein,specificity=REPORTERSPECIFIC)
        
      points(A[unspecific.spectra.sel],M[unspecific.spectra.sel],pch=pch.p,col="yellow",...)
      points(A[groupspecific.spectra.sel],M[groupspecific.spectra.sel],pch=pch.p,col="orange",...)
      points(A[specific.spectra.sel],M[specific.spectra.sel],pch=pch.p,col="red",...)
        
    }
)

maplot.protein <- function(x,relative.to,protein,noise.model=NULL,
        channels=NULL,xlim=NULL,ylim=NULL,identify=FALSE,add=FALSE,pchs=NULL,log="xy",
        legend.pos="topright",names=NULL,legend.cex=0.8,cols=pchs,ltys=NULL,
        main=protein,xlab=NULL,ylab=NULL,type="ma",...) {
 	  if (is(x,"IBSpectra"))
      x <- c(x)
   
    if (!relative.to %in% reporterTagNames(x[[1]]))
      stop("relative.to must be one of the reporter tag names (",paste(reporterTagNames(x[[1]]),collapse=","),")")

    i.df <- data.frame()
    legend <- list(text=c(),lty=c(),pch=c())
    i <- 1
    if (is.null(channels)) channels <- setdiff(reporterTagNames(x[[1]]),relative.to)
    if (is.null(pchs)) pchs <- seq_along(channels)

    if (is.null(xlab) && type=="ma") xlab <- "average intensity"
    if (is.null(ylab) && type=="ma") ylab <- "ratio"

    if (is.null(xlab) && type!="ma") xlab <- paste("intensity",relative.to)
    if (is.null(ylab) && type!="ma") ylab <- paste("intensity",paste(channels,collapse=", "))
    

    for (ib in x) {
      ions <- reporterIntensities(ib,na.rm=FALSE,protein=protein)
      if (length(ions) == 0 || all(is.na(ions)))
        next;

      if (is.null(xlim))
        xlim <- range(ions,na.rm=TRUE)
      if (is.null(ylim) && type!="ma") ylim <- xlim
      if (any(!is.finite(xlim))) xlim  <- NULL
      if (type == "ma")
        ions[which(is.na(ions))] <- 0

      for (channel_i in seq_along(channels)) {
        channel <- channels[channel_i]
        channel.rt <- ifelse(length(relative.to) == 1,relative.to,relative.to[channel_i])
        if (length(ions[,channel])>0) {
          if (is.null(names))
            legend$text <- c(legend$text,sprintf("%s/%s",channel,channel.rt))
          else
            legend$text <- c(legend$text,names[i])
          legend$pch <- c(legend$pch,pchs[i])
          legend$lty <- c(legend$pch,ltys[i])
        }

        if (type=="ma") {
          div <- ions[,channel]/ions[,channel.rt]
          avg <- (ions[,channel]+ions[,channel.rt])/2
          div[div == Inf] = 100
          div[div == -Inf] = 0.01
          div[div > 100] = 90
          div[div < 0.01] = 1/90
          x <- avg
          y <- div
        } else {
          x <- ions[,channel.rt]
          y <- ions[,channel]
        }

        if (channel_i == 1 && add == FALSE) 
          plot(x,y,xlim=xlim,ylim=ylim,pch=pchs[i],log=log,main=main,col=cols[i],
               xlab=xlab,ylab=ylab,...)
        else
          points(x,y,pch=pchs[i],col=cols[i],...)

        add = TRUE
          if (identify)
            i.df <- rbind(i.df,data.frame(x=x,y=y,spectrum=rownames(ions),
                                          peptide=peptides(x=ib,spectrum=rownames(ions)),stringsAsFactors=T))

        if (!is.null(noise.model)) {
          ratio <- estimateRatio(ib,noise.model,channel.rt,channel,protein=protein)
          ratio.lm <- estimateRatio(ib,noise.model,channel.rt,channel,protein=protein,method="lm")
          #abline(h=10^ratio[1],lty=ltys[i],col=cols[i],...)

          if (type=="ma") {
            abline(h=10^ratio[1],lty=1,col=cols[i],...)
            abline(h=ratio.lm['ratio'],lty=2,col=cols[i],...)
          } else {
            #abline(0,10^ratio[1],lty=1,col=cols[i],...)
            abline(0,10^ratio[1],lty=1,col=cols[i],untf=TRUE,...)
            abline(0,ratio.lm['ratio'],lty=2,col=cols[i],untf=TRUE,...)
          }

          if (!is.na(ratio[1]))
            legend$text[length(legend$text)] <- 
              sprintf("%s = %.2f +/- %.2f",legend$text[length(legend$text)],10^ratio[1],10^sqrt(ratio[2]))
          #abline(h=10^(ratio[1]+1.96*sqrt(ratio[2])),lty=channel_i,lwd=0.5)
          #abline(h=10^(ratio[1]-1.96*sqrt(ratio[2])),lty=channel_i,lwd=0.5)
        }
        i <- i + 1
      }
    }
    if (length(legend$text) > 0 ) {
    legend(legend.pos,legend=legend$text,pch=legend$pch,lty=legend$lty,cex=legend.cex,col=cols)
    if (identify) {
      identify(x=i.df$x,y=i.df$y,labels=i.df$peptide)
    }}
}

.get.ri <- function(ri,ch) {
  if (ch == "ALL")
    rowSums(ri,na.rm=TRUE)
  else if (ch == "AVG")
    apply(ri,1,mean,na.rm=TRUE)
  else
    ri[,ch]
}

# MA plots
setMethod("maplot",
          signature(x="IBSpectra",channel1="character",channel2="character"),
          function(x,channel1,channel2,noise.model=NULL,colorize.protein=NULL,
              h=NULL,v=NULL,col.h="green",col.v="green",
              ylab=paste("ratio channel",channel1,"vs",channel2),
              pch=".",noise.model.col=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999"),
              pch.p=1,protein=NULL,peptide=NULL,smooth=FALSE,
              set.na.to=NULL,set.na.to.lim=NULL,na.rm=is.null(set.na.to),
              identify=FALSE,identify.column="spectrum",
              x.axis.labels=TRUE,y.axis.labels=TRUE,
              col="black",...) {
            if (is.null(protein)) 
              ions <- reporterIntensities(x,na.rm=FALSE)
            else
              ions <- reporterIntensities(x,na.rm=FALSE,protein=protein)

            fd <- fData(x)
            
            # remove points which are NA in both channels
            if (channel2=="ALL" | channel1=="ALL")
              sel <- apply(is.na(ions),1,all)
            else 
              sel <- is.na(ions[,channel1]) & is.na(ions[,channel2])
            ions <- ions[!sel,,drop=FALSE]
            fd <- fd[!sel,,drop=FALSE]
            
            if (na.rm) {
              if (channel2=="ALL" | channel1=="ALL")
                sel <- apply(is.na(ions),1,any)
              else 
                sel <- is.na(ions[,channel1])| is.na(ions[,channel2])
              if (any(sel))
                warning(sprintf("removing %s NA points",sum(sel)))
              ions <- ions[!sel,,drop=FALSE]
              fd <- fd[!sel,,drop=FALSE]
            }
            
            if (channel1=="ALL" && channel2=="ALL")
              stop("Cannot set both channels to ALL")

            if (channel1=="ALL") {
              R=as.numeric(ions[,setdiff(colnames(ions),channel2)])
              G=rep(ions[,channel2],ncol(ions)-1)
            } else if (channel2=="ALL") {
              R=rep(ions[,channel1],ncol(ions)-1)
              G=as.numeric(ions[,setdiff(colnames(ions),channel1)])
            } else {
              R=ions[,channel1]
              G=ions[,channel2]
            }
            
            M <- R/G
            A <- log10(sqrt(R*G))
            

            if (isTRUE(set.na.to)) {
              set.na.to <- ceiling(10^(max(abs(log10(M)), na.rm = TRUE)+0.2))
              set.na.to.lim <- ceiling(10^(max(abs(log10(M)), na.rm = TRUE)+0.1))
            }
            if (!is.null(set.na.to)) {
              sel <- is.na(R);
              M[sel] <- set.na.to + 1;
              A[sel] <- log10(G[sel])
              sel <- is.na(G); 
              M[sel] <- 1/set.na.to
              A[sel] <- log10(R[sel])
            }

            if (length(M)==0) {
              plot(1,1,cex=0,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
              text(1,1,"no datapoints")
              return()
            }
            
            if (smooth) {
#              smoothScatter(A,M,log="y",...)
              stop("option smooth deprecated!")
            } else {
              plot(A,M,pch=pch,log="y",
                  xlab=expression(paste(log[10]," average intensity")),
                  ylab=ylab,col=col,axes=FALSE,...)
              axis(side=1,labels=x.axis.labels)

              if (!is.null(set.na.to)) {
                y.axis <- axTicks(2,log=TRUE)
                y.axis <- y.axis[y.axis >= 1/set.na.to.lim & y.axis <= set.na.to.lim]
                y.axis <- y.axis[seq(from=2,to=length(y.axis)-1)]
                if (isTRUE(y.axis.labels))
                  y.axis.labels <- c(-Inf,y.axis,Inf)
                axis(side=2,at=c(1/set.na.to.lim,y.axis,set.na.to.lim),
                     labels= y.axis.labels,las=2)
                abline(h = y.axis, col = "#000000F0",lwd=0.15)
                abline(h = c(1/set.na.to.lim,set.na.to.lim), col = "blue", lty = 3)
              } else {
                axis(side=2,labels=y.axis.labels)
              }

            }
            if (!is.null(h)) {
              sel <- M > h | M < 1/h
              points(A[sel],M[sel],col=col.h,pch=pch,...)
            }
            if (!is.null(v)) {
              sel <- A < v
              points(A[sel],M[sel],col=col.v,pch=pch,...)
            }
            
            #abline(h=0,col="red")
            if (!missing(noise.model) & !is.null(noise.model)) {
              if (is(noise.model,"NoiseModel"))
                noise.model <- c(noise.model)

              ss <- seq(min(A,na.rm=TRUE),max(A,na.rm=TRUE),length.out=100)
              for (i in seq_along(noise.model)) {          
                lines(ss,10^(1.96*sqrt(variance(noise.model[[i]],ss))),col=noise.model.col[i],lwd=2)
                lines(ss,10^(-1.96*sqrt(variance(noise.model[[i]],ss))),col=noise.model.col[i],lwd=2)
              }
            }
            
            if (!is.null(colorize.protein)) {
              if (is.list(colorize.protein)) {
                cp.legend <- c()
                cp.legend.col <- c()
                for (i in names(colorize.protein)) {
                   data <- colorize.protein[[i]]
                   protein <- data$protein
                   p.col <- data$col
                   p.us.col <- data$us.col
                   cp.legend <- c(cp.legend,i)
                   cp.legend.col <- c(cp.legend.col,p.col)

                   protein.sel <- names(A) %in% spectrumSel(x,protein=protein,
                                                            specificity=c(GROUPSPECIFIC,REPORTERSPECIFIC),
                                                            spectrum.titles=TRUE)

                   protein.us.sel <- names(A) %in% spectrumSel(x,protein=protein,specificity=UNSPECIFIC,
                                                                spectrum.titles=TRUE)


                   points(A[protein.sel],M[protein.sel],pch=pch.p,col=p.col,...)
                   points(A[protein.us.sel],M[protein.us.sel],pch=pch.p,col=p.us.col,...)
                }
                legend("topright",legend=cp.legend,col=cp.legend.col,pch=pch.p)
              } else {
                unspecific.spectra.sel    <- names(A) %in% spectrumSel(x,protein=colorize.protein,specificity=UNSPECIFIC,spectrum.titles=T)
                groupspecific.spectra.sel <- names(A) %in% spectrumSel(x,protein=colorize.protein,specificity=GROUPSPECIFIC,spectrum.titles=T)
                specific.spectra.sel      <- names(A) %in% spectrumSel(x,protein=colorize.protein,specificity=REPORTERSPECIFIC,spectrum.titles=T)
                other.sel <- !unspecific.spectra.sel & !groupspecific.spectra.sel & !specific.spectra.sel
              
                points(A[unspecific.spectra.sel],M[unspecific.spectra.sel],pch=pch.p,col="yellow",...)
                points(A[groupspecific.spectra.sel],M[groupspecific.spectra.sel],pch=pch.p,col="orange",...)
                points(A[specific.spectra.sel],M[specific.spectra.sel],pch=pch.p,col="green",...)
  #points(A[other.sel],M[other.sel],pch=".",col="blue",...)
              }
              
            }
            if (identify)
              identify(x=A,y=M,labels=fd[,identify.column])
          }
)






setMethod("maplot",
    signature=c(x="missing",channel1="numeric",channel2="numeric"),
    function(channel1,channel2,noise.model=NULL,pch=".",noise.model.col=1:10,...) {
      sel <- !is.na(channel1) & !is.na(channel2)
      channel1 <- channel1[sel]
      channel2 <- channel2[sel]

      if (length(channel1)==0) {
        plot(1,1,cex=0,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        text(1,1,"no datapoints")
        return()
      }
      
      M <- log10(channel1) - log10(channel2)
      A <- 0.5*(log10(channel1) + log10(channel2))
      M <- M[order(A)]
      ch1 <- log10(channel1)[order(A)]
      ch2 <- log10(channel2)[order(A)]
      A <- sort(A)
      plot(A,M,pch=pch,xlab=expression(paste(log[10]," average intensity")),
          ylab="ratio",...)
      #abline(h=0,col="red")
      if (!missing(noise.model) & !is.null(noise.model)) {
        if (is(noise.model,"NoiseModel"))
          noise.model <- c(noise.model)
        ss <- seq(min(A),max(A),by=0.05)
        for (i in seq_along(noise.model)) {          
          lines(ss,1.96*sqrt(variance(noise.model[[i]],ss)),col=noise.model.col[i],lwd=1.5)
          lines(ss,-1.96*sqrt(variance(noise.model[[i]],ss)),col=noise.model.col[i],lwd=1.5)
          #lines(A,1.96*sqrt(variance(noise.model[[i]],channel1,channel2)),col=noise.model.col[i])
          #lines(A,-1.96*sqrt(variance(noise.model[[i]],channel1,channel2)),col=noise.model.col[i])
        }
      }
    }
)

setMethod("maplot",
          signature(x="IBSpectra",channel1="missing",channel2="missing"),
          function(x,noise.model=NULL,pairs=TRUE,xlim="fixed",ylim="fixed",...){
            ions <- reporterIntensities(x)
              set.na.to <- NULL
              set.na.to.lim <- NULL
              if (identical(xlim,"fixed")) 
                xlim <- log10(range(ions,na.rm=T)) 
              if (identical(ylim,"fixed")) {
                y.max <- max(apply(ions,1,function(x) {
                                   y <- x[!is.na(x)]
                                   if (length(y) < 2) NA
                                   else max(y)/min(y)
                                   }),na.rm=T)
                set.na.to <- ceiling(10^(log10(y.max)+0.2))
                set.na.to.lim <- ceiling(10^(log10(y.max)+0.1))
                y.max <- ceiling(10^(log10(y.max)+0.3))
                ylim <- c(1/y.max,y.max)

              }

            if (pairs) {
              histlimits=log10(range(ions,na.rm=TRUE))
              par(mfrow=c(ncol(ions),ncol(ions)),mar=c(2.5,2.25,0,0),
                  bty="n",fg="darkgray")
              for (i in colnames(ions)) {
                for (j in colnames(ions)) {
                  if (i == j) { 
                    if (all(is.na(ions[,i])))
                      plot(histlimits,c(1,1),type="n",bty="n", xaxt="n", yaxt="n", main="",ylim=c(0,1))
                    else
                      hist(log10(ions[,i]), col="#EEEEEE", freq=FALSE, 
                           xlim=histlimits,cex=0.5,breaks=20,
                           bty="n", xaxt="n", yaxt="n", main="",ylim=c(0,1))
                    text(sum(histlimits)/2,0.5,i,col="black",cex=1.2,font=2); 
                    #text(sum(histlimits)/2,0.4,sprintf("n = %s",sum(!is.na(ions[,i]))),col="black",cex=1,font=2); 
                    #legend("center",
                    #       legend=c(sprintf("%3s",i),
                    #                sprintf("n = %s",sum(!is.na(ions[,i])))),
                    #       bty="n",text.col="black",
                    #       cex=1,text.font=c(2,1))

                  } else if (i > j) {
                    plot.labels.x = is.null(xlim) || (i == tail(colnames(ions),1) && j == head(colnames(ions),1))
                    plot.labels.y = is.null(ylim) || (i == tail(colnames(ions),1) && j == head(colnames(ions),1))
                    maplot(x=x,channel1=as.character(i),channel2=as.character(j),
                           noise.model=noise.model,xlim=xlim,ylim=ylim,
                           set.na.to=set.na.to,set.na.to.lim=set.na.to.lim,
                           x.axis.labels = plot.labels.x,
                           y.axis.labels = plot.labels.y,
                           ...)
                    
                  }
                  else if (i < j) {
                    plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")
                  }
                        
                }
              }
            } else {
              par(mfrow=c(1,ncol(ions)),mar=c(2.5,2.25,2.25,0),
                  bty="n",fg="darkgray")
              for (i in colnames(ions)) {
                plot.labels.y = is.null(ylim) || (i == colnames(ions)[1])
                maplot(x=x,channel1=as.character(i),channel2="ALL",
                       main=substitute(paste(i," ",italic(vs)," ALL",sep="")),
                       noise.model=noise.model,xlim=xlim,ylim=ylim,
                       set.na.to=set.na.to,set.na.to.lim=set.na.to.lim,
                       x.axis.labels = TRUE,
                       y.axis.labels = plot.labels.y,
                       ...)

              }
            }
          }
)



# Should protGgdata be deprecated?
setMethod("protGgdata",
    signature(x="ANY",relative.to="character",protein="character"),
    function(x,relative.to,protein,noise.model=NULL,
        channels=setdiff(reporterTagNames(x),relative.to),
        names=NULL,legend.cex=0.8,...) {
      
      dfs <- data.frame()
      i <- 1
      for (ib in x) {
            ions <- reporterIntensities(ib,na.rm=FALSE,protein=protein)
            ions[which(is.na(ions))] <- 0
            
            for (channel_i in seq_along(channels)) {
              channel <- channels[channel_i]
              channel.rt <- ifelse(length(relative.to) == 1,
                                   relative.to,relative.to[channel_i])
              
              div <- ions[,channel]/ions[,channel.rt]
              avg <- (ions[,channel]+ions[,channel.rt])/2
              
              if (length(div)>0 & sum(avg,na.rm=T) > 0) {
                if (!is.null(noise.model)) {
                  ratio <- estimateRatio(ib,noise.model,channel.rt,channel,protein=protein)
                  
                  dfs <-
                    rbind(dfs,data.frame(channel=channel,channel.rt=channel.rt,
                                         average=avg,difference=div,ratio=10^ratio[1],
                                         name=ifelse(is.null(names),
                                          sprintf("%s/%s",channel,channel.rt),names[i]),
                                         lower=10^(ratio[1]-sqrt(ratio[2])),
                                         upper=10^(ratio[1]+sqrt(ratio[2]))))
                  
                }
              }
              i <- i + 1
            }
          }
          return(dfs)
          
    }
)


