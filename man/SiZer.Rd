\name{SiZer}
\alias{SiZer}

\title{SiZer  plot for 1-dimensional data}
\description{
  SiZer (\bold{Si}gnificant \bold{Zer}o crossings) plot for 1-dimensional data.
}

\usage{
SiZer(x, bw, gridsize, scaleData=FALSE, signifLevel=0.05,
      plotSiZer=TRUE, logbw=TRUE, xlim, xlab, 
      addLegend=TRUE, posLegend="topright") 
}

\arguments{
  \item{x}{data vector}
  \item{bw}{vector of range of bandwidths}
  \item{gridsize}{number of x- and y-axis grid points}
  \item{scaleData}{flag for scaling the data i.e. transforming to unit
    variance for each dimension.}
  \item{signifLevel}{significance level}
  \item{plotSiZer}{flag for displaying SiZer map}
  \item{logbw}{flag for displaying log bandwidths on y-axis}
  \item{xlim}{x-axis limits}
  \item{xlab}{x-axis label}
  \item{addLegend}{flag for legend display}
  \item{posLegend}{legend position}
}
  

\value{
  SiZer plot sent to graphics window.
}

\details{
  The gradient SiZer map of Chaudhuri & Marron (1999) is
  implemented. The horizontal axis is the data axis, the
  vertical axis are the bandwidths.
  %It returns a list with the following fields \cr
  %\code{x.grid} - vector of grid points \cr
  %\code{bw} - vector of bandwidths at grid points\cr
  %\code{SiZer} - matrix (rows = grid points, columns = bandwidths) for
  %SiZer map:  3 = decreasing gradient (red), 
  %2 = increasing gradient (blue), 1 = zero gradient (purple), 0 =
  %sparse region (grey).     
}
 

\seealso{\code{\link{featureSignif}}}

\references{
  Chaudhuri, P. & Marron, J.S. (1999) 
  SiZer for exploration of structures in curves.
  \emph{Journal of the American Statistical Association},
  \bold{94},  807-823.
}

\examples{
data(earthquake)
eq3 <- -log10(-earthquake[,3])
SiZer(eq3, xlab="-log(-depth)")
SiZer(eq3, xlab="-log(-depth)", bw=c(0.05, 1.3), gridsize=401)
}

\keyword{hplot}
