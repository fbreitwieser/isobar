if(!isGeneric("as.data.frame")) setGeneric("as.data.frame", useAsDefault=as.data.frame)

#TODO: unify factor.as.character and factor.to.chr
.factor.as.character <- function(df) {
  for (col_i in seq_len(ncol(df))) {
    if (is.factor(df[,col_i]))
      df[,col_i] <- as.character(df[,col_i])
  }
  df
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

.as.vect <- function(matrix,col.data=2,col.names=1) {
  vect <- matrix[,col.data]
  names(vect) <- matrix[,col.names]
  vect
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

      if (trim < 0 | trim > 0.5)
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


