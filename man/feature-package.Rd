\name{feature}
\alias{feature}
\docType{package}
\title{
feature
}
\description{
  Package for feature significance for multivariate kernel density estimation (for
  1- to 4-dimensional data). 
  
}
\details{
  There is one main function in this package,
  \code{\link{featureSignif}}. It has a range of options which allow
  the user to display compute and kernel density estimates, significant gradient
  and significant curvature regions. Significant gradient and/or
  curvature regions often correspond to significant features (e.g.
  local modes). 

  It's available in an interactive mode where the user selects
  bandwidths in the graphics window and the significant regions are computed in
  real-time. This interactive mode is useful for initial exploratory
  data analysis.

  Otherwise there is the non-interactive mode where the user can specify
  a single set of bandwidths. This is useful when the user has a more
  definite idea of suitable bandwidths.     

  For a for comprehensive example for 1-d and 2-d data, see
  \code{vignette("feature")}. 
}


\author{
Tarn Duong <tduong@maths.unsw.edu.au> & Matt Wand <wand@maths.unsw.edu.au>
}
\references{
  Chaudhuri, P. and Marron, J.S. (1999) 
  SiZer for exploration of structures in curves.
  \emph{Journal of the American Statistical Association},
  \bold{94},  807-823.

  Duong, T., Cowling, A., Koch, I., Wand, M.P. (2006)
  Feature significance for multivariate kernel density estimation.
  \emph{Submitted.}

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
\code{\link[sm:sm-package]{sm}}, \code{\link[KernSmooth]{KernSmooth}}
}

}
