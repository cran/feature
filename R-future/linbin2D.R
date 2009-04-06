######### R-function: linbin2D #########
 
# Creates the grid counts from a bivariate data set X 
# over an equally-spaced set of grid points
# contained in "gpoints" using the linear 
# binning strategy. Note that the FORTRAN subroutine
# "lbtwod" is called. 

# Last changed: 25 AUG 1995

linbin2D <- function(X,gpoints1,gpoints2)
{
   n <- nrow(X)
   X <- c(X[,1],X[,2]) 
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   out <- .Fortran("lbtwod",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(b1),as.double(b2),
           as.integer(M1),as.integer(M2),double(M1*M2),
           PACKAGE="feature")
   return(matrix(out[[9]],M1,M2))
}

########## End of linbin2D ##########
