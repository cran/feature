\name{drvkde}
\alias{drvkde}
\title{Kernel density derivative estimation}

\description{
  Compute kernel density derivative
  estimates and standard errors for multivariate data. 
}

\usage{
drvkde(x, drv, bandwidth, gridsize, range.x, binned=FALSE, se=TRUE)
}

\arguments{
  \item{x}{data matrix or matrix of binning counts} 
  \item{drv}{vector of derivative indices}
  \item{bandwidth}{vector of bandwidths}
  \item{gridsize}{vector of grid sizes}
  \item{range.x}{list of vector of ranges for x}
  \item{binned}{flag to indicate: TRUE = x is binned counts or FALSE = x is
    data matrix. Default is TRUE} 
  \item{se}{flag for computing standard error of kernel
    estimate. Default is TRUE}    
}

\value{
  Returns a list with fields
  
  \code{x.grid}  - grid points \cr
  \code{est} - kernel estimate of partial derivative of density function
  indicated by \code{drv}\cr
  \code{se} -  estimate of standard error of \code{est} (if \code{se=TRUE}).

}

\details{The estimates and standard errors are computed over a grid of binned counts
  \code{x.grid}. If the binned counts are not supplied then they are computed
  inside this function.

  If \code{gridsize} and \code{range.x} are not supplied, they are
  computed inside this function.
  
  This function doesn't need to be used directly as it called from
  \code{\link{featureSignif}}.
}


\examples{
## univariate
data(earthquake)
eq3 <- -log10(-earthquake[,3])
fhat <- drvkde(x=eq3, drv=0, bandwidth=0.1)    ## KDE of f
fhat1 <- drvkde(x=eq3, drv=1, bandwidth=0.1)   ## KDE of df/dx

## trivariate
data(earthquake)
earthquake[,3] <- -log10(-earthquake[,3])
fhat <- drvkde(x=earthquake, drv=c(0,0,0), bandwidth=c(0.04,0.04,0.05),
               se=FALSE)     ## KDE of f
                                              
fhat212 <- drvkde(x=earthquake, drv=c(2,1,2), bandwidth=c(0.04,0.04,0.05),
                  se=FALSE)  ## KDE of d^3 f/ dx^2 dy dz^2
}

\references{
     Wand, M.P. and Jones, M.C. (1995) \emph{Kernel Smoothing} Chapman and Hall.
   }
   
\author{M.P. Wand  \email{wand@maths.unsw.edu.au}}

\seealso{
\code{\link{featureSignif}}
}

\keyword{smooth} 

