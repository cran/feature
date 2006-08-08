########## R-function: linbin ##########

# For application of linear binning to a 
# univariate data set.

# Last changed: 16 JUNE 1995

linbin <- function(X,gpoints,truncate=TRUE)

{
   n <- length(X)
   M <- length(gpoints)  
   trun <- 0
   if (truncate) trun <- 1
   a <- gpoints[1]
   b <- gpoints[M]
   out <- .Fortran("linbin",as.double(X),as.integer(n),
           as.double(a),as.double(b),as.integer(M),
           as.integer(trun),double(M),PACKAGE="feature")
   return(out[[7]])
}

########## End of linbin ##########
