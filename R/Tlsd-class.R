## Class: TlsParameter
setClass("TlsParameter",
          representation = representation(df = "numeric",
                                          location = "numeric",
                                          scale = "numeric"),
          prototype = prototype(df = 1, location = 0, scale = 1, name =
                      gettext("Parameter of a Student/T location scale distribution")
                      ),
          contains = "Parameter"
          )

setClass("Tlsd",
          prototype = prototype(
                      r = function(n) { rt(n, df)*scale + location },
                      d = function(x, log = FALSE){ 
                            1/scale * dt((x - location)/scale, df, log = log) },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){
                            pt((x - location)/scale, df,lower.tail=lower.tail,log.p=log.p)},
                      q = function(q, lower.tail = TRUE, log.p = FALSE ){
                            qt(p, df,lower.tail=lower.tail,log.p=log.p)*scale + location},
                      param = new("TlsParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE,
                     Symmetry = new("SphericalSymmetry",
                                     type = "univariate symmetric distribution",
                                     SymmCenter = 0)
                      ),
          contains = "AbscontDistribution"
          )

## TODO: log=TRUE does not work correctly
## Class: TlsDistribution
setMethod("initialize", "Tlsd",
          function(.Object, df = 1, location = 0, scale = 1) {
            .Object@img <- new("Reals")
            .Object@param <- new("TlsParameter", df = df, location = location,
                                  scale = scale)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){}
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){}
            
            rtls = function(n,df,location,scale) rt(n, df)*scale + location 
            dtls = function(x,df,location,scale,log = FALSE)
                     1/scale * dt((x - location)/scale, df, log = log) 
            ptls = function(q,df,location,scale,lower.tail = TRUE, log.p = FALSE )
                     pt((q - location)/scale, df,lower.tail=lower.tail,log.p=log.p) 
                     
            qtls = function(p,df,location,scale, lower.tail = TRUE, log.p = FALSE )
                            qt(p, df,lower.tail=lower.tail,log.p=log.p)*scale + location
                            
            body(.Object@r) <- substitute(
                           { rtls(n, df=dfSub,location=locationSub,scale=scaleSub) },
                             list(dfSub = df,
                                  locationSub = location,
                                  scaleSub = scale)
                                          )
            body(.Object@d) <- substitute(
                           { dtls(x, df=dfSub, location = locationSub,
                                     scale = scaleSub, log = log) },
                             list(dfSub = df,locationSub = location, scaleSub = scale)
                                          )
            body(.Object@p) <- substitute(
                           { ptls(q, df=dfSub, location = locationSub,
                                     scale = scaleSub, lower.tail = lower.tail,
                                     log.p = log.p) },
                             list(dfSub = df,locationSub = location, scaleSub = scale)
                                          )
            body(.Object@q) <- substitute(
                           { qtls(p, df=dfSub, location = locationSub,
                                     scale = scaleSub, lower.tail = lower.tail,
                                     log.p = log.p) },
                             list(dfSub = df,locationSub = location, scaleSub = scale)
                                          )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- FALSE
            .Object@Symmetry <- SphericalSymmetry(location)
            .Object
          })


## Access Methods
setMethod("df", "TlsParameter", function(x, ...) x@df)
setMethod("location", "TlsParameter", function(object) object@location)
setMethod("scale", "TlsParameter",
           function(x, center = TRUE, scale = TRUE) x@scale)
## Replace Methods
setReplaceMethod("df", "TlsParameter",
                  function(object, value){ object@df <- value; object})
setReplaceMethod("location", "TlsParameter",
                  function(object, value){ object@location <- value; object})
setReplaceMethod("scale", "TlsParameter",
                  function(object, value){ object@scale <- value; object})

setValidity("TlsParameter", function(object){
  if(length(df(object)) != 1)
    stop("df has to be a numeric of length 1")
  if(df(object) <= 0)
    stop("df has to be positive")
  # TODO: more tests
  else return(TRUE)
})

################################
##
## Class: Student T Location Scale distribution
##
################################

Tlsd <- function(df = 1, location = 0, scale = 1) new("Tlsd", df = df, location = location, scale = scale)

## wrapped access methods
setMethod("df", "Tlsd", function(x, ...) df(param(x)))
setMethod("location", "Tlsd", function(object) location(param(object)))
setMethod("scale", "Tlsd",
           function(x, center = TRUE, scale = TRUE) scale(param(x)))
           ### odd arg-list due to existing function in base package 


## wrapped replace methods
#setMethod("df<-", "Tlsd",
#           function(object, value) new("Tlsd", df = value, ncp = ncp(object)))
