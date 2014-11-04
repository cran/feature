\name{featureSignif}
\alias{featureSignif}
\title{Feature significance for kernel density estimation}

\description{
  Identify significant features of kernel density estimates of 1- to
  4-dimensional data. 
}

\usage{
featureSignif(x, bw, gridsize, scaleData=FALSE, addSignifGrad=TRUE,
   addSignifCurv=TRUE, signifLevel=0.05)  
}
		  

\arguments{
  \item{x}{data matrix}
  \item{bw}{vector of bandwidth(s)}
  \item{gridsize}{vector of estimation grid sizes}
  \item{scaleData}{flag for scaling the data i.e. transforming to unit
    variance for each dimension.}
  \item{addSignifGrad}{flag for computing significant gradient regions}
  \item{addSignifCurv}{flag for computing significant curvature regions}
  \item{signifLevel}{significance level}
}

\value{
  Returns an object of class \code{fs} which is a list with the following fields
  \describe{
  \item{\code{x}}{data matrix}
  \item{\code{names}}{name labels used for plotting}
  \item{\code{bw}}{vector of bandwidths}
  \item{\code{fhat}}{kernel density estimate on a grid}
  \item{\code{grad}}{logical grid for significant gradient}
  \item{\code{curv}}{logical grid for significant curvature}
  \item{\code{gradData}}{logical vector for significant gradient data points}
  \item{\code{gradDataPoints}}{significant gradient data points}
  \item{\code{curvData}}{logical vector for significant curvature data points}
  \item{\code{curvDataPoints}}{significant curvature data points}
  } 
}

\details{
  Feature significance is based on significance testing of the gradient
  (first derivative) and curvature (second derivative) of a kernel
  density estimate.
  This was developed for 1-d data by Chaudhuri & Marron (1995), for
  2-d data by Godtliebsen, Marron & Chaudhuri (1999), and for 3-d and 4-d
  data by Duong, Cowling, Koch & Wand (2007).

  The test statistic for gradient testing is at a point \eqn{\mathbf{x}}{x} is
  \deqn{W(\mathbf{x}) = \Vert \widehat{\nabla f} (\mathbf{x}; \mathbf{H}) \Vert^2}{%
    W(x) = || hat{grad f}(x; H)||^2} where
  \eqn{\widehat{\nabla f} (\mathbf{x};\mathbf{H})}{hat{grad f}(x; H)}
  is kernel estimate of the gradient
  of \eqn{f(\mathbf{x})}{f(x)} with bandwidth \eqn{\mathbf{H}}{H}, and
  \eqn{\Vert\cdot\Vert}{||.||} is the Euclidean norm. \eqn{W(\mathbf{x})}{W(x)} is
  approximately chi-squared distributed with \eqn{d} degrees of freedom
  where \eqn{d} is the dimension of the data.

  The analogous test statistic for curvature is  
  \deqn{W^{(2)}(\mathbf{x}) = \Vert \mathrm{vech}
    \widehat{\nabla^{(2)}f} (\mathbf{x}; \mathbf{H})\Vert ^2}{%
    W2(x) = ||vech hat{curv f}(x; H)||^2}
  where \eqn{\widehat{\nabla^{(2)} f} (\mathbf{x};\mathbf{H})}{hat{curv f}(x;
    H)} is the kernel estimate of the curvature of
  \eqn{f(\mathbf{x})}{f(x)}, and vech is the vector-half operator.
  \eqn{W^{(2)}(\mathbf{x})}{W2(x)} is
  approximately chi-squared distributed with \eqn{d(d+1)/2}{d(d+1)/2} degrees of freedom.

  Since this is a situation with many dependent hypothesis tests, we
  use a multiple comparison or simultaneous test to control the
  overall level of significance. We use a Hochberg-type procedure. See
  Hochberg (1988) and Duong, Cowling, Koch & Wand (2007).
}


\references{
  Chaudhuri, P. & Marron, J.S. (1999) 
  SiZer for exploration of structures in curves.
  \emph{Journal of the American Statistical Association},
  \bold{94},  807-823.

  Duong, T., Cowling, A., Koch, I. & Wand, M.P. (2008)
  Feature significance for multivariate kernel density estimation.
  \emph{Computational Statistics and Data Analysis}, \bold{52}, 4225-4242. 

  Hochberg, Y. (1988) A sharper Bonferroni procedure for multiple tests
  of significance. \emph{Biometrika}, \bold{75}, 800-802.
  
  Godtliebsen, F., Marron, J.S. & Chaudhuri, P. (2002) 
  Significance in scale space for bivariate density estimation.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{11}, 1-22.

  Wand, M.P. & Jones, M.C. (1995)
  \emph{Kernel Smoothing.} Chapman & Hall/CRC, London.
}



\examples{
## Univariate example
data(earthquake)
eq3 <- -log10(-earthquake[,3])
fs <- featureSignif(eq3, bw=0.1)
plot(fs, addSignifGradRegion=TRUE)

## Bivariate example
library(MASS)
data(geyser)
fs <- featureSignif(geyser)
plot(fs, addKDE=FALSE, addData=TRUE)  ## data only
plot(fs, addKDE=TRUE)                 ## KDE plot only
plot(fs, addSignifGradRegion=TRUE)    
plot(fs, addKDE=FALSE, addSignifCurvRegion=TRUE)
plot(fs, addSignifCurvData=TRUE, curvCol="cyan")
}

\seealso{\code{\link{featureSignifGUI}}, \code{\link{plot.fs}}}

\keyword{smooth} 

