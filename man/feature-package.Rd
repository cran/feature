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

\keyword{ package }

\seealso{
\code{ks}, \code{sm}, \code{KernSmooth}
}

