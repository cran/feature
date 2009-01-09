######### R-function: linbin4D #########
 
# Creates the grid counts from a quadrivariate data set X 
# over an equally-spaced set of grid points
# contained in "gpoints" using the linear 
# binning strategy. Note that the FORTRAN subroutine
# "lbfoud" is called. 

# Last changed: 31 AUG 2005

linbin4D <- function(X,gpoints1,gpoints2,gpoints3,gpoints4)
{
   n <- nrow(X)
   X <- c(X[,1],X[,2],X[,3],X[,4]) 
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3)
   M4 <- length(gpoints4)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   a4 <- gpoints4[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   b4 <- gpoints4[M4]
   out <- .Fortran("lbfoud",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(a3),as.double(a4),
           as.double(b1),as.double(b2),as.double(b3),as.double(b4),
           as.integer(M1),as.integer(M2),as.integer(M3),as.integer(M4),
           double(M1*M2*M3*M4),PACKAGE="feature")
   return(array(out[[15]],c(M1,M2,M3,M4)))
}

########## End of linbin4D ##########

