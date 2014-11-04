\name{plot.fs}
\alias{plot.fs}

\title{Feature signficance plot for 1- to 3-dimensional data}
\description{
  Feature signficance plot for 1- to 3-dimensional data.
}

\usage{
\method{plot}{fs}(x, ..., xlab, ylab, zlab, xlim, ylim, zlim, add=FALSE, addData=FALSE,
   scaleData=FALSE, addDataNum=1000, addKDE=TRUE,jitterRug=TRUE,
   addSignifGradRegion=FALSE, addSignifGradData=FALSE,
   addSignifCurvRegion=FALSE, addSignifCurvData=FALSE, addAxes3d=TRUE,
   densCol, dataCol="black", gradCol="green4", curvCol="blue",
   axisCol="black", bgCol="white", dataAlpha=0.1, gradDataAlpha=0.3,
   gradRegionAlpha=0.2, curvDataAlpha=0.3, curvRegionAlpha=0.3)
}

\arguments{
  \item{x}{object of class \code{fs} (output from
    \code{\link{featureSignif}} function)}
  \item{xlim,ylim,zlim}{x-, y-, z-axis limits}
  \item{xlab,ylab,zlab}{x-, y-, z-axis labels}
  \item{scaleData}{flag for scaling the data i.e. transforming to unit
    variance for each dimension}
  \item{add}{flag for adding to an existing plot}
  \item{addData}{flag for display of the data}
  \item{addDataNum}{maximum number of data points plotted in
    displays}
  \item{addKDE}{flag for display of kernel density estimates}
  \item{jitterRug}{flag for jittering of rug-plot for univariate 
    data display}
  \item{addSignifGradRegion,addSignifGradData}{flag for display of significant gradient regions/data points}
  \item{addSignifCurvRegion,addSignifCurvData}{flag for display of significant curvature regions/data points}
  \item{addAxes3d}{flag for displaying axes in 3-d displays}
  \item{densCol}{colour of density estimate curve}
  \item{dataCol}{colour of data points}
  \item{gradCol}{colour of significant gradient regions/data points}
  \item{curvCol}{colour of significant curvature regions/data points}
  \item{axisCol}{colour of axes}
  \item{bgCol}{colour of background}
  \item{dataAlpha}{transparency of data points}
  \item{gradRegionAlpha,gradDataAlpha}{transparency of significant gradient regions/data points}
  \item{curvRegionAlpha,curvDataAlpha}{transparency of significant curvature regions/data points}
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of 1-d and 2-d kernel density estimates are sent to graphics
  window. Plot for 3-d is sent to RGL window. 
}
 

\seealso{\code{\link{featureSignif}}}

\examples{
## See ? featureSignif for uni- and bivariate examples
## Trivariate example
data(earthquake)
earthquake[,3] <- -log10(-earthquake[,3])
fs <- featureSignif(earthquake, scaleData=TRUE, bw=c(0.06, 0.06, 0.05))
plot(fs, addKDE=FALSE, addSignifCurvRegion=TRUE)
}

\keyword{ hplot}
