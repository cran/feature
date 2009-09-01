\name{feature}
\alias{feature}
\docType{package}
\title{
feature
}
\description{
  Package for feature significance for multivariate kernel density estimation. 
  
}
\details{
  The \pkg{feature} package contains functions 
  to display and compute kernel density estimates, significant gradient
  and significant curvature regions. Significant gradient and/or
  curvature regions often correspond to significant features (e.g.
  local modes). 

  There are two main functions in this package. 
  \code{\link{featureSignifGUI}} is the interactive function where
  the user can select bandwidths from a pre-defined range. This
  mode is useful for initial exploratory
  data analysis. \code{\link{featureSignif}} is the non-interactive function.
  This is useful when the user has a more
  definite idea of suitable values for the bandwidths.     
  For a more detailed example for 1-d and 2-d data, see
  \code{vignette("feature")}. 
}


\author{
Tarn Duong <tarn.duong@gmail.com> & Matt Wand <wand@uow.edu.au>
}
\references{
  Chaudhuri, P. and Marron, J.S. (1999) 
  SiZer for exploration of structures in curves.
  \emph{Journal of the American Statistical Association},
  \bold{94},  807-823.

  Duong, T., Cowling, A., Koch, I., Wand, M.P. (2008)
  Feature significance for multivariate kernel density estimation.
  \emph{Computational Statistics and Data Analysis}, \bold{52}, 4225-4242. 

  Godtliebsen, F., Marron, J.S. and Chaudhuri, P. (2002) 
  Significance in scale space for bivariate density estimation.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{11}, 1-22.
  
  Scott, D.W. (1992) \emph{Multivariate Density Estimation: Theory,
    Practice, and Visualization}. John Wiley \& Sons. New York.

  Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman \&
  Hall/CRC. London.
}

\keyword{ package }

\seealso{
\code{sm}, \code{KernSmooth}
}

