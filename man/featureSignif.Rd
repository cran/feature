\name{featureSignif}
\alias{featureSignif}
\title{Feature significance for kernel density estimation}

\description{
  Identify significant features of kernel density estimates of 1- to
  4-dimensional data. This feature has an interactive and a
  non-interactive mode.
}

\usage{
featureSignif(x, bw, xlab, ylab, zlab, xlim, ylim, zlim,
   addData=FALSE, scaleData=FALSE, addDataNum=1000,
   addKDE=TRUE, jitterRug=TRUE, signifLevel=0.05, 
   addSignifGradRegion=FALSE, addSignifGradData=FALSE,
   addSignifCurvRegion=FALSE, addSignifCurvData=FALSE,
   plotSiZer=FALSE, logbwSiZer=TRUE, addAxes3d=TRUE,
   densCol, dataCol="black", gradCol="green4", curvCol="blue",
   axisCol="black", bgCol="white", dataAlpha=0.1, gradDataAlpha=0.3,
   gradRegionAlpha=0.2, curvDataAlpha=0.3, curvRegionAlpha=0.3,
   gridsize, gridsizeSiZer)
}
		  

\arguments{
  \item{x}{data matrix}
  \item{bw}{bandwidth(s) - see below for details on how to specify bandwidths}
  \item{xlim, ylim, zlim}{x-, y-, z-axis limits}
  \item{xlab, ylab, zlab}{x-, y-, z-axis labels}
  \item{scaleData}{flag for scaling the data i.e. transforming to unit
    variance for each dimension. Default is FALSE.}
  \item{addData}{flag for display of the data. Default is FALSE.}
  \item{addDataNum}{maximum number of data points plotted in
    displays. Default is 1000.}
  \item{addKDE}{flag for display of kernel density estimates. Default
    is TRUE. Not available for 4-d data.}
  \item{jitterRug}{flag for jittering of rug-plot for univariate 
    data display. Default is TRUE.}
  \item{addSignifGradRegion}{flag for display of significant gradient regions.
    Default is FALSE. Not available for 4-d data.}
  \item{addSignifGradData}{flag for display of significant gradient
    data points. Default is FALSE.}
  \item{addSignifCurvRegion}{flag for display of significant curvature regions.
    Default is FALSE. Not available for 4-d data.}
  \item{addSignifCurvData}{flag for display of significant curvature
    data points. Default is FALSE.}
  \item{plotSiZer}{flag for displaying 1-d gradient SiZer map. Default is FALSE.}
  \item{logbwSiZer}{flag for displaying log bandwidths in SiZer
    map. Default is TRUE.}
  \item{addAxes3d}{flag for displaying axes in 3-d displays. Default is
    TRUE.}
  \item{signifLevel}{significance level. Default is 0.05.}
  \item{densCol}{colour of density estimate curve. Default for 1-d data
    "DarkOrange", for 2-d data is heat.colors(1000), for
    3-d data is heat.colors(5).}
  \item{dataCol}{colour of data points. Default is
    "black".}
  \item{gradCol}{colour of significant gradient regions/points. Default
    is "green4".}
  \item{curvCol}{colour of significant curvature regions/points.
    Default is "blue".}
  \item{axisCol}{colour of axes. Default is "black".}
  \item{bgCol}{colour of background. Default is "white".}
  \item{dataAlpha}{alpha-blending transparency value for data points.}
  \item{gradDataAlpha}{alpha-blending transparency value for
    significant gradient data points.}
  \item{gradRegionAlpha}{alpha-blending transparency value for
    significant gradient regions.}
  \item{curvDataAlpha}{alpha-blending transparency value for
    significant curvature data points.}
  \item{curvRegionAlpha}{alpha-blending transparency value for
    significant curvature regions.}
  \item{gridsize}{vector of the number of grid points in each
    direction.}
  \item{gridsizeSiZer}{number of x-axis grid points for SiZer
    map. Default is 101.} 
}

\value{
  If \code{bw} is not specified, then a range of possible bandwidths is
  automatically calculated. 
  
  For univariate data, \code{bw} can be either a scalar or a
  vector. With the former, a 
  KDE is computed with this scalar bandwidth. The latter is interpreted
  as a range of bandwidths.

  For multivariate data, \code{bw} can either be a vector or a matrix.
  With the former, a KDE is computed with this vector bandwidth. The
  latter is interpreted as a range of bandwidths with the first row are
  the minimum values and the second row the maximum values. 

  If a range of bandwidths is supplied, it goes into interactive mode. 
  If a single bandwidth is supplied, it goes into non-interactive mode.

  Returns a list with the following fields\cr
  \code{x} - data matrix \cr
  \code{bw} - vector of bandwidths \cr
  \code{fhat} - kernel density estimate on a grid (output from
  \code{\link{drvkde}}) \cr
  \code{grad} - logical array indicating significant gradient (if
  \code{addSignigGradRegion=TRUE})\cr 
  \code{curv} - logical array indicating significant curvature (if
  \code{addSignigCurvRegion=TRUE})\cr 
  \code{gradData} - logical vector indicating significant gradient (if
  \code{addSignigGradData=TRUE})\cr 
  \code{curvData} - logical vector indicating significant curvature (if
  \code{addSignigCurvData=TRUE}) 

  In the interactive case, the return values are
  based on the last bandwidths chosen before the interactive session
  was ended. In the non-interactive case, the return values are based on the
  specified bandwidth.
  
  For 1-d data, the gradient SiZer map  of Chaudhuri \& Marron (1999) is
  implemented. If this option is selected, it automatically goes into
  non-interactive mode. The horizontal axis is the data axis, the
  vertical axis are the bandwidths. It returns a list with the following
  fields \cr 
  \code{x.grid} - vector of grid points \cr
  \code{bw} - vector of bandwidths at grid points\cr
  \code{SiZer} - matrix (rows = grid points, columns = bandwidths) for
  SiZer map:  3 = decreasing gradient (red), 
  2 = increasing gradient (blue), 1 = zero gradient (purple), 0 =
  sparse region (grey).
  
}

\details{
  Feature significance is based on significance testing of the gradient
  (first derivative) and curvature (second derivative) of a kernel
  density estimate.
  This was developed for 1-d data by Chaudhuri \& Marron (1995), for
  2-d data by Godtliebsen, Marron \& Chaudhuri (1999), and for 3-d and 4-d
  data by Duong, Cowling, Koch \& Wand (2006).

  The test statistic for gradient testing is at a point \eqn{\mathbf{x}}{x} is
  \deqn{W(\mathbf{x}) = \Vert \widehat{\nabla f} (\mathbf{x}; \mathbf{H}) \Vert^2}{%
    W(x) = || hat{grad f}(x; H)||^2} where
  \eqn{\widehat{\nabla f} (\mathbf{x};\mathbf{H})}{hat{grad f}(x; H)}
  is kernel estimate of the gradient
  of \eqn{f(\mathbf{x})}{f(x)} with bandwidth \eqn{\mathbf{H}}{H}, and
  \eqn{\Vert\cdot\Vert}{||.||} is the Euclidean norm. \eqn{W(\mathbf{x})}{W(x)} is
  approximately chi-squared distributed with \eqn{d} degrees of freedom
  where \eqn{d} is the dimension of the data.

  The test statistic for curvature is analogous to that for gradient
  testing:
  \deqn{W^{(2)}(\mathbf{x}) = \Vert \mathrm{vech}
    \widehat{\nabla^{(2)}f} (\mathbf{x}; \mathbf{H})\Vert ^2}{%
    W2(x) = ||vech hat{curv f}(x; H)||^2}
  where \eqn{\widehat{\nabla^{(2)} f} (\mathbf{x};\mathbf{H})}{hat{curv f}(x;
    H)} is the kernel estimate of the curvature of
  \eqn{f(\mathbf{x})}{f(x)}, and vech is the vector-half operator.
  \eqn{W^{(2)}(\mathbf{x})}{W2(x)} is
  approx. chi-squared distributed with \eqn{d(d+1)/2} degrees of freedom.

  Since this is a situation with many dependent hypothesis tests, we
  use a multiple comparison or simultaneous test to control the
  overall level of significance. We use a Hochberg-type procedure. See
  Hochberg (1988) and Duong, Cowling, Koch \& Wand (2006).
}
  


\references{
  Chaudhuri, P. and Marron, J.S. (1999) 
  SiZer for exploration of structures in curves.
  \emph{Journal of the American Statistical Association},
  \bold{94},  807-823.

  Duong, T., Cowling, A., Koch, I., Wand, M.P. (2006)
  Feature significance for multivariate kernel density estimation.
  \emph{Submitted.}

  Hochberg, Y. (1988) A sharper Bonferroni procedure for multiple tests
  of significance. \emph{Biometrika}, \bold{75}, 800-802.
  
  Godtliebsen, F., Marron, J.S. and Chaudhuri, P. (2002) 
  Significance in scale space for bivariate density estimation.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{11}, 1-22.

  Wand, M.P. and Jones, M.C. (1995)
  \emph{Kernel Smoothing} Chapman and Hall.
}


\seealso{
\code{\link[KernSmooth]{bkde}} (in package `KernSmooth'),
\code{\link[KernSmooth]{bkde2D}} (in package `KernSmooth'),
\code{\link[stats]{density}}
}
\examples{
## Non-interactive examples

## Univariate example
data(earthquake)
eq3 <- -log10(-earthquake[,3])

featureSignif(eq3, addSignifGradRegion=TRUE,xlab="-log(-depth)", bw=0.1)
featureSignif(eq3, plotSiZer=TRUE, xlab="-log(-depth)")

## Bivariate example

library(MASS)
data(geyser)

fs <- featureSignif(geyser, addSignifGradRegion=TRUE,
     addSignifCurvRegion=TRUE, bw=c(4.5, 0.37))
names(fs)

## Trivariate example

data(earthquake)
earthquake[,3] <- -log10(-earthquake[,3])

featureSignif(earthquake, scaleData=TRUE, addData=TRUE,
    bw=c(0.0381, 0.0381, 0.0442))

featureSignif(earthquake, addKDE=FALSE, scaleData=TRUE,
   addSignifGradRegion=TRUE, addSignifCurvRegion=TRUE,
   bw=c(0.0381, 0.0381, 0.0442),
   xlim=c(0.4,0.5), ylim=c(0.4,0.5), zlim=c(0.8,0.9))


\dontrun{
## Interactive examples

library(MASS)
data(geyser)
duration <- geyser$duration 

## Univariate example

featureSignif(duration)
featureSignif(duration, addData=TRUE)
featureSignif(duration, addSignifGradRegion=TRUE, addSignifGradData=TRUE)
featureSignif(duration, addSignifCurvRegion=TRUE, addSignifCurvData=TRUE)

## Bivariate example

featureSignif(geyser, addData=TRUE, addSignifGradRegion=TRUE,
   addSignifGradData=TRUE, bw=rbind(c(1, 0.1), c(5, 0.9)))
   ## bandwidths ranges: h1 in c(1, 5), h2 in c(0.1, 0.9)

## Trivariate example

data(earthquake)
earthquake$depth <- -log10(-earthquake$depth)
featureSignif(earthquake, addSignifGradRegion=TRUE, scaleData=TRUE)

}
}

\keyword{smooth} 

