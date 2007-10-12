\name{plot.fs}
\alias{plot.fs}

\title{Feature signficance plot for 1- to 3-dimensional data}
\description{
  Feature signficance plot for 1- to 3-dimensional data.
}

\synopsis{\method{plot}{fs}(x, ...)}

\usage{
\method{plot}{fs}(x, xlab, ylab, zlab, xlim, ylim, zlim,
   addData=FALSE, scaleData=FALSE, addDataNum=1000,
   addKDE=TRUE, jitterRug=TRUE,  
   addSignifGradRegion=FALSE, addSignifGradData=FALSE,
   addSignifCurvRegion=FALSE, addSignifCurvData=FALSE,
   addAxes3d=TRUE, densCol, dataCol="black", gradCol="green4",
   curvCol="blue", axisCol="black", bgCol="white",
   dataAlpha=0.1, gradDataAlpha=0.3, gradRegionAlpha=0.2,
   curvDataAlpha=0.3, curvRegionAlpha=0.3)
}

\arguments{
  \item{x}{an object of class \code{fs} (output from
    \code{\link{featureSignif}} function)}
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
  \item{addAxes3d}{flag for displaying axes in 3-d displays. Default is
    TRUE.}
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
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of 1-d and 2-d kernel density estimates are sent to graphics
  window. Plot for 3-d is sent to RGL window. 
}

\references{
  Chaudhuri, P. and Marron, J.S. (1999) 
  SiZer for exploration of structures in curves.
  \emph{Journal of the American Statistical Association},
  \bold{94},  807-823.

  Duong, T., Cowling, A., Koch, I., Wand, M.P. (2007)
  Feature significance for multivariate kernel density estimation.
  \emph{Submitted.}
  
  Godtliebsen, F., Marron, J.S. and Chaudhuri, P. (2002) 
  Significance in scale space for bivariate density estimation.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{11}, 1-22.
}

\details{
  This is a plot method for \code{fs} objects created from
  \code{featureSignif}. This plotting function is called
  automatically from inside
  \code{featureSignif} except when \code{plotFS=FALSE}.
  If the user creates an \code{fs} object with the
  significant gradient and curvature, then it is easier and more
  efficient to modify the graphical display just using \code{plot.fs}. See
  examples below.    
  
}
 

\seealso{\code{\link{featureSignif}}}

\examples{
library(MASS)
data(geyser)

fs <- featureSignif(geyser, addSignifGradRegion=TRUE,
     addSignifCurvRegion=TRUE, bw=c(4.5, 0.37), plotFS=FALSE)

plot(fs, addKDE=FALSE, addData=TRUE)  ## data only
plot(fs, addKDE=TRUE)                 ## KDE plot only
plot(fs, addSignifGradRegion=TRUE)    
plot(fs, addKDE=FALSE, addSignifCurvRegion=TRUE)
plot(fs, addSignifCurvData=TRUE, curvCol="cyan")
}

\keyword{ hplot}
