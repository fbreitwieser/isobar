if(!isGeneric("as.data.frame")) setGeneric("as.data.frame", useAsDefault=as.data.frame)

# for R version pre 2.14
#if (!exists("paste0")) # will show a warning on newer R versions
  paste0 <- function(...,sep="") paste(...,sep=sep)

"%inrange%" <- function(a,b) {
  if (!is.numeric(a) || !is.numeric(b)) stop("Arguments must be numeric")
  if (length(b) != 2) stop("Second argument must have 2 elements")
  return(a >= b[1] & a <=b[2])
}

.cn <- function(x,y) {
  x[ , y ]
}

.gg_theme <- function(...) {
  if ( compareVersion(packageDescription("ggplot2")$Version,"0.9.1") > 0 )
    theme(...)
  else
    opts(...)
}

.gg_element_text <- function(...) {
  if ( compareVersion(packageDescription("ggplot2")$Version,"0.9.1") > 0 )
    element_text(...)
  else
    theme_text(...)
}

.abbrev <- function(strings,n=4,collapse=NULL) {
  ll <- length(strings)
  if (ll > n) {
    strings <- strings[seq_len(n-1)]
    strings[n] <- sprintf("... (%s in total)",ll)
  }
  if (!is.null(collapse))
    paste(strings,collapse=collapse)
  else
    strings
}


# are 'most' (default: 90%) of the values in x TRUE?
.most <- function(x,fraction=0.9) {
  if (is.null(dim(x)) || !is.logical(x)) stop(".most function works on logical vectors")
  sum(x)/length(x) > fraction
}

.check.isfunction <- function(f) {
  if (!is.function(f))
    stop(paste(deparse(substitute(f)),"must be a function!"))
  TRUE
}

.unique.or.collapse <- function(x,collapse=";") {
  if (is.null(x)) NA
  else
    ifelse(all(x==x[1]),x[1],paste0(x,collapse=collapse))
}

.paste_unique <- function(x,...,na.rm=TRUE) {
  x <- unique(x)
  x <- x[!is.na(x)]
  paste(x,...)
}

.grep_columns <- function(df,pattern,...,logical=TRUE) {
  if (logical)
    grepl(pattern,colnames(df),...)
  else
    grep(pattern,colnames(df),...)
}

#TODO: unify factor.as.character and factor.to.chr
.factor.as.character <- function(df) {
  for (col_i in seq_len(ncol(df))) {
    if (is.factor(df[,col_i]))
      df[,col_i] <- as.character(df[,col_i])
  }
  df
}

.names.as.vector <- function(x) {
  if (!is.null(names(x)))
    vec <- names(x)
  else 
    vec <- x
  setNames(vec,x)
}

# from Gavin Simpson [http://stackoverflow.com/questions/9788026/change-the-order-of-columns]
.moveToFirstCol <- function(df, colname) {
  cnams <- colnames(df)
  want <- which(colname == cnams)
  df[, c(cnams[want], cnams[-want])]
}

.factor.to.chr <- function(df) {
  for (col in colnames(df))
    if (is.factor(df[,col])) df[,col] <- as.character(df[,col])
  df
}

.all.duplicates <- function(x,n=2) {
  t <- table(x)
  names(t)[t>=n]
}

# get a number range. E.g. 1,2,3,5,6 -> 1-3,5,6
.string.number.ranges <- function(numbers) {
  n <- number.ranges(numbers)
  if (is.na(n)) return("")
  else return(sprintf("[%s]",n))
}

number.ranges <- function(numbers) {
  if (all(is.na(numbers))) { return(NA) }
  numb=c()
  numb_string=c()

  for (i in sort(unique(as.numeric(numbers)))) {
    if (length(numb) == 0) {
      numb = c(i)
    } else {
      if (i == numb[length(numb)]+1) {
        numb = c(numb,i)
      } else {
        if (length(numb) <= 2) {
          numb_string = c(numb_string,numb)
        } else {
          numb_string = c(numb_string,paste(min(numb),max(numb),sep="-"))
        }
        numb = c(i)
      }
    }
  }
  if (length(numb) <= 2) {
        numb_string = c(numb_string,numb)
  } else {
       numb_string = c(numb_string,paste(min(numb),max(numb),sep="-"))
  }
  return(paste(numb_string,collapse=","))
}

.return.equal.or.na <- function(df) {
  apply(df, 1, function(x) {

    if (all(is.na(x))) return(NA)

    y <- x[!is.na(x)]
    if (!all(y == y[1])) return(NA)
    return(y[1])
  })
}

.all.duplicate.rows <- function(df,column,n=2) {
  t <- table(df[,column])
  res <- df[df[,column] %in% names(t)[t>=n],]
  res[order(res[,column]),]
}

.vector.as.data.frame <- function(vect,colnames=NULL,stringsAsFactors=FALSE) {
	as.data.frame(.as.matrix(vect,colnames),stringsAsFactors=stringsAsFactors)
}

.as.matrix <- function(vect,colnames=NULL) {
  mat <- matrix(c(names(vect),vect),ncol=2,byrow=FALSE)
  colnames(mat) <- colnames
  mat
}

.as.vect <- function(my.matrix,col.data=2,col.names=1) {
  setNames(my.matrix[,col.data],my.matrix[,col.names])
}

.stopifnot <- function(cond,...) if (!cond) stop(...)

.stopiflengthnotequal <- function(x,y,...) {
  if (length(x) != length(y)) {
    stop(..., " length(x) = ",length(x),", length(y) = ",length(y),"")
  }
}


.stopifna  <- function(data,...) if (any(is.na(data))) stop(...)


.expand.w.vector <- function(v,m,by.m,v.name='v',...) {
  df <- do.call(rbind,lapply(names(v),function(x) {
    pos <- which(m[,by.m] %in% v[x])
    data.frame(x,m[pos,],row.names=NULL,stringsAsFactors=FALSE,...)
  }))
  colnames(df)[1] <- v.name
  df
}

.trim <- function(x,side="both") {
  if (side=="left" | side=="both")
    x <- sub("^\\s*","",x)
  if (side=="right" | side=="both")
    x <- sub("\\s*$","",x)
  x
}

.bnds <- function(x,bnd=NULL,min.x=-bnd,max.x=bnd) {
  max(
      min(x,max.x),
      min.x)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### weightedVariance and weightedMean generics and functions
### replace by Hmisc functions?

setGeneric("weightedMean", function(data, weights,trim=0) standardGeneric("weightedMean"))
setGeneric("weightedVariance", function(data, weights, mean.estimate,trim=0) standardGeneric("weightedVariance"))

setMethod("weightedVariance",
    signature(data = "numeric", weights = "numeric", mean.estimate = "missing"),
    function(data, weights, trim=0) {
      mean.estimate <- weightedMean(data,weights,trim=trim)
      weightedVariance(data,weights=weights,mean.estimate=mean.estimate,trim=trim)
    }
)

setMethod("weightedVariance",
    signature(data = "numeric", weights = "numeric", mean.estimate = "numeric"),
    function(data, weights, mean.estimate, trim=0) {
      # Validation? if (length(data)==1) { return(0); }
      # Should we rescale the weights? weights=1/sum(weights)
      
      # weighted sample variance - for not normally distributed data
      # see http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html and
      # http://en.wikipedia.org/wiki/Weighted_mean#Weighted_sample_variance
      sel <- !is.na(data) & !is.na(weights)
      weights <- weights[sel]
      data <- data[sel]

      if (trim < 0 || trim > 0.5)
        stop("trim has to be between 0 and 0,5")

      if (trim > 0) {
        sel <- data > quantile(data,trim) & data < quantile(data,1-trim)
        weights <- weights[sel]
        data <- data[sel]
      }
      
      V1 <- sum(weights)
      V2 <- sum(weights**2)
      variance <- ( V1 / (V1**2 - V2) ) * sum(weights*(data-mean.estimate)**2)
      
      return(variance)
    }
)

setMethod("weightedMean",
    signature(data = "numeric", weights = "numeric"),
    function(data, weights, trim = 0) {
      sel <- !is.na(data) & !is.na(weights)
      weights <- weights[sel]
      data <- data[sel]

      if (trim < 0 | trim > 0.5)
        stop("trim has to be between 0 and 0,5")

      if (trim > 0) {
        sel <- data > quantile(data,trim) & data < quantile(data,1-trim)
        weights <- weights[sel]
        data <- data[sel]
      }
 
      return(
          sum(data * weights) / sum(weights)
      )
    })


# split string into named vector
.strsplit_vector <- function(x,pattern) {
  pos <- sapply(x,regexpr,pattern="=",fixed=TRUE)
  res <- substring(x,pos+1)
  names(res) <- substring(x,rep(1,length(x)),pos-1)
  res
}

.sum.bool.na  <- function(x) 
  c('TRUE'=sum(x,na.rm=TRUE),'FALSE'=sum(!x,na.rm=TRUE),'TRUE %'=round(sum(x,na.rm=TRUE)/length(x)*100,1),length=length(x),n.na=sum(is.na(x)))

# provide summary of boolean vector values
.sum.bool  <- function(x) 
  c('TRUE'=sum(x),'FALSE'=sum(!x),'TRUE %'=round(sum(x)/length(x)*100,1),length=length(x))

# provide character summary of boolean vectors
.sum.bool.c  <- function(x) 
  paste('TRUE: ',sum(x),'; FALSE: ',sum(!x),'; TRUE %: ',round(sum(x)/length(x)*100,1),"; length: ",length(x),sep="")

.combine.fisher <- function(p.values,signs) {
  if (length(signs) != length(p.values))
    stop("lratios and pvalues must have equal length!")
  
  sel.notna <- !is.na(p.values)
  p.values <- p.values[sel.notna]
  signs <- signs[sel.notna]
  k <- length(p.values)

  if (length(p.values) == 1) return(p.values)

  ## require that the direction is the same for all p-values
  if (!all(signs == signs[1])) return(1)

  return(pchisq(-2*sum(log(p.values)),2*k,lower.tail=FALSE))
}

.combine.fisher.tblwide <- function(my.df) {
  lr.cols <- grepl("^lratio.",colnames(my.df))
  p.cols <- grepl("^p.value.rat.",colnames(my.df)) & !grepl("adj",colnames(my.df))

  if (sum(lr.cols) != sum(p.cols))
    stop("unequal number of '^lratio.' and '^p.value.rat' columns")

  combined.p <- rep(1,nrow(my.df))
  signs.equal <- apply(sign(my.df[,lr.cols]),1,function(x) { y=x[!is.na(x)]; all(y==y[1])})
  ks <- apply(!is.na(my.df[,p.cols]),1,sum)
  logsums <- rowSums(log(my.df[,p.cols]),na.rm=TRUE)
  sel <-  ks > 1 & signs.equal
  combined.p[sel] <- pchisq(-2*logsums[sel],2*ks[sel],lower.tail=FALSE)

  combined.p[ks==1] <- as.numeric(apply(my.df[ks==1,p.cols],1,function(x) { y=x[!is.na(x)]; ifelse(length(y)==0,NA,y) } ))
  return(combined.p)
}

.call.cmd <- function(cmd,stdout.to=NULL) {
    if (is.null(stdout.to)) {
      message("  calling system command [",cmd,"]",appendLF=FALSE)
      if (system(cmd) != 0) stop("\nError executing [",cmd,"]")
      message(" finished.")
    } else {
      message("  calling system command [",cmd," > ",stdout.to,"]",appendLF=FALSE)
      if (system(paste(cmd,">",stdout.to)) != 0) 
        stop("\nError executing [",cmd,"]: \n\n ...\n",
             paste(tail(readLines(stdout.to),n=10),collapse="\n"))
      message(" finished.")
    }
}


.get.cmbn <- function(combn,tags,cl) {
  if (!all(unlist(combn) %in% cl))
    stop("incorrect combn specification")

  res <- c()
  for (cc in combn)
    for (tag1 in tags[cl==cc[1]&!is.na(cl)]) 
      for (tag2 in tags[cl==cc[2]&!is.na(cl)])
        res <- cbind(res,c(r1=tag1,r2=tag2,
                     class1=cc[1],
                     class2=cc[2]))
  res
}


.sanitize.sh <- function(str) {
  gsub("[^a-zA-Z\\.0-9_\\-]","", str)
}

.weighted.cor <- function( x, y, w = rep(1,length(x)),use='complete.obs') {
  ## (c) Heather Turner, Vincent Zoonekynd at http://stackoverflow.com/questions/9460664/weighted-pearsons-correlation
  stopifnot(length(x) == dim(y)[2] )
  if (use=='complete.obs') {
    sel <- !is.na(x) & !is.na(y) & !is.na(w)
    x <- x[sel]
    y <- y[sel]
    w <- w[sel]
  }

  w <- w / sum(w)
  # Center x and y, using the weighted means
  x <- x - sum(x * w)
  ty <- t(y - colSums(t(y) * w))
  # Compute the variance
  vx <- sum(w * x * x)
  vy <- colSums(w * ty * ty)
  # Compute the covariance
  vxy <- colSums(ty * x * w)
  # Compute the correlation
  vxy / sqrt(vx * vy)
}


.concensus.il.peptide <- function(peptide) { 
    pep <- do.call(rbind,strsplit(peptide,""))
    paste0(apply(pep,2,function(x) {
      xu <- unique(x)
      if (length(xu)==1) xu
      else "L"
    }),collapse="")
  }

.fix.il.peptide <- function(from,sub.il=TRUE) {
   if (.PEPTIDE.COLS['REALPEPTIDE'] %in% colnames(from))
     return(from)

   from[,.PEPTIDE.COLS['REALPEPTIDE']] <- from[,.SPECTRUM.COLS['PEPTIDE']]
   l.peptide <- gsub("I","L",from[,.SPECTRUM.COLS['PEPTIDE']])
   if (sub.il) {
     from$peptide <- l.peptide
   } else {
     from$peptide <- as.vector( tapply(from[,.SPECTRUM.COLS['PEPTIDE']],l.peptide,function(x) {
           if (all(x == x[1])) x[1]
           else 
              .concensus.il.peptide(unique(x))
     }))[l.peptide]
   }
   return(from)
}


