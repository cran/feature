\name{SiZer,SiCon}
\alias{SiZer}
\alias{SiCon}

\title{SiZer and SiCon plots for 1-dimensional data}
\description{
  SiZer (\bold{Si}gnificant \bold{Zer}o crossings) and SiCon (\bold{Si}gnificant \bold{Con}vexity) plots for 1-dimensional data.
}

\usage{
SiZer(x, bw, gridsize, scaleData=FALSE, signifLevel=0.05, plotSiZer=TRUE,
   logbw=TRUE, xlim, xlab, addLegend=TRUE, posLegend="bottomright") 

SiCon(x, bw, gridsize, scaleData=FALSE, signifLevel=0.05, plotSiCon=TRUE,
   logbw=TRUE, xlim, xlab, addLegend=TRUE, posLegend="bottomright")
}

\arguments{
  \item{x}{data vector}
  \item{bw}{vector of range of bandwidths}
  \item{gridsize}{number of x- and y-axis grid points}
  \item{scaleData}{flag for scaling the data i.e. transforming to unit
    variance for each dimension.}
  \item{signifLevel}{significance level}
  \item{plotSiZer,plotSiCon}{flag for displaying SiZer/SiCon map}
  \item{logbw}{flag for displaying log bandwidths on y-axis}
  \item{xlim}{x-axis limits}
  \item{xlab}{x-axis label}
  \item{addLegend}{flag for legend display}
  \item{posLegend}{legend position}
}
  

\value{
  SiZer/SiCon plot sent to graphics window.
}

\details{
  The gradient SiZer and curvature SiCon maps of Chaudhuri & Marron (1999) are
  implemented. The horizontal axis is the data axis, the
  vertical axis are the bandwidths. The colour scheme for the SiZer map is 
  red: negative gradient, blue: positive gradient, purple: zero gradient
  and grey: sparse regions. For the SiCon map, 
  orange: negative curvature (concave), blue: positive curvature (convex), 
  green: zero curvature and grey: sparse regions.
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
SiZer(eq3)
SiCon(eq3)
}

\keyword{hplot}
