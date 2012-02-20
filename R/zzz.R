# ==========================================================================
# isobar package initialization
# ==========================================================================

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("Welcome to isobar (v ",packageDescription("isobar")$Version,")\n",
                "   'openVignette(\"isobar\")' and '?isobar' provide help on usage.\n", sep=""))
   addVigs2WinMenu("isobar") 
}
