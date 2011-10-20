# ==========================================================================
# isobar package initialization
# ==========================================================================
.onLoad <- function(libname, pkgname) {
    ## need contents to load at library attach - not at build time
    #.initContents() ## in environment.R
    .buildIsobarOpts() ## in environment.R
}

.buildIsobarOpts <- function() {
    BioC <- getOption("BioC")
	 #TODO: assert that BioC is not null

	 ## add isobar specific options
    Isobar <- BioC$Isobar
    if (is.null(Isobar)) {
        Isobar <- list()
        class(Isobar) <- "BioCPkg"
    }

	if (is.null(Isobar$algo.foo)) Isobar$algo.foo = "bla"

    BioC$Isobar <- Isobar
    options("BioC"=BioC)

}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("\nWelcome to Isobar package\n",
                "Vignettes contain introductory material. To view, type",
                "'openVignette(\"isobar\")'. To cite Isobar, see",
                "'citation(\"isobar\")'.\n", sep="\n  "))
   addVigs2WinMenu("isobar") 
}

.onUnload <- function( libpath ) {
  #library.dynam.unload( "isobar", libpath )
}
