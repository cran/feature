######### R-function: linbin3D #########
 
# Creates the grid counts from a trivariate data set X 
# over an equally-spaced set of grid points
# contained in "gpoints" using the linear 
# binning strategy. Note that the FORTRAN subroutine
# "lbthrd" is called. 

# Last changed: 27 JUL 2005

linbin3D <- function(X,gpoints1,gpoints2,gpoints3)
{
   n <- nrow(X)
   X <- c(X[,1],X[,2],X[,3]) 
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3) 
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   out <- .Fortran("lbthrd",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(a3),as.double(b1),
           as.double(b2),as.double(b3),as.integer(M1),as.integer(M2),
           as.integer(M3),double(M1*M2*M3),PACKAGE="feature")
   return(array(out[[12]],c(M1,M2,M3)))
}

########## End of linbin3D ##########

