\name{drvkde}
\alias{drvkde}
\title{Kernel density derivative estimation}

\description{
  Compute kernel density derivative
  estimates and standard errors for multivariate data. 
}

\usage{
drvkde(x,drv,bandwidth,gridsize,range.x,binned=FALSE,se=TRUE)
}

\arguments{
  \item{x}{data matrix or matrix of binning counts} 
  \item{drv}{vector of derivative indices}
  \item{bandwidth}{vector of bandwidths}
  \item{gridsize}{vector of grid sizes}
  \item{range.x}{vector of ranges for x}
  \item{binned}{TRUE if x is binned counts or FALSE (default) if x is
    data matrix} 
  \item{se}{TRUE (default) computes standard error (SE) for kernel
    estimate or FALSE skips computing SE}    
}

\value{
  If \code{se=TRUE}, it returns a list with fields
  
  \code{x.grid}  - grid points \cr
  \code{est} - kernel estimate of partial derivative of density function
  indicated by \code{drv}\cr
  \code{se} -  estimate of standard error of \code{est}.

  Otherwise if
  \code{se=FALSE}, only \code{x.grid} and \code{est} are returned.  
}

\details{The estimates and standard errors are computed over a grid of binned counts
  \code{x.grid}. If the binned counts are not supplied then they are computed
  inside this function.

  This function doesn't need to be used directly as it called from
  \code{\link{featureSignif}}.
}


\references{
     Wand, M.P. and Jones, M.C. (1995) \emph{Kernel Smoothing} Chapman and Hall.
   }
   
\author{M.P. Wand  \email{wand@maths.unsw.edu.au}}

\seealso{
\code{\link{featureSignif}}
}

\keyword{smooth} 

