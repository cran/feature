#.First.lib <- function(lib, pkg=feature) {
#  library.dynam("feature", pkg, lib)
#  x <- installed.packages()
#  cat(pkg, x[x[,1]==pkg,3], "installed\n")
#  cat("Copyright T. Duong & M. P. Wand (2006)\n")
#  rm(x)
#  require(rgl)
#  require(misc3d)
#}

.onLoad <- function(libname=NULL, pkgname=feature)
{
  #x <- installed.packages()
  #cat(pkgname, x[x[,1]==pkgname,3], "(2006)\n")
  #rm(x)
  cat("feature 1.1-5 (2006)\n") 
}
