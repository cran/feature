\name{featureSignifGUI}
\alias{featureSignifGUI}
\title{GUI for feature significance for kernel density estimation}

\description{
  GUI for feature significance for kernel density estimation.
}

\usage{
featureSignifGUI(x, scaleData=FALSE)
}
		  

\arguments{
  \item{x}{data matrix}
  \item{scaleData}{flag for scaling the data to the unit interval in
	each dimension}
}

%\value{
%  The result is the output from \link{featureSignif}, an object of class
%  \code{fs}. 
%}

\details{
  In the first column are the sliders for selecting the bandwidths (one
  for each dimension). Move
  the slider buttons to change the value of the bandwidths. 
  The text field is for the grid size which specifies the number of
  points in each dimension of the kernel estimation binning grid. Press the `Compute significant features' button to begin the
  computation. This creates a plot of the kernel density estimate (KDE)
  from the data with the specified bandwidths by calling
  \code{\link{featureSignif}}. Once this complete, a
   pop-up window will appear. 
  
  In the second column are the axis limits and labels. The last text
  field is for the (maximum) number of data points used in the display.
  Press the `Reset plot (except KDE)' button to clear the plot of all
  added features except for the KDE itself. 
  
  In the third column are 5 buttons which can be used to add to the KDE plot
  such as the data points, significant gradient points/regions and
  significant curvature points/regions. 
  For 1-d data, the button in the third column is `Compute SiZer
  map'. Press this button to compute a gradient SiZer plot using the
  \code{\link{SiZer}} function. Once this complete, a pop-up window will appear. 
  For 2- and 3-d data, the button in the third column is 'Reset plot'. This
  will clear the plot of all features as well as the KDE. This is useful
  for showing only the significant features when the KDE
  may interfere with their display.

  For 3-d data, there is an extra fourth column of options: these are
  sliders for the transparency values for the 
  features. Move the slider button along to the desired value (between 0
  and 1) and then press the `Add ...' button to the left. 
  Repeatedly pressing the `Add ...' button will cause the transparency
  of the features to decrease. In this case, press the one of the 'Reset plot'
  buttons to clear the plot window, and replot the significant feature
  with the desired transparency. 
}
  

\examples{
\dontrun{
library(MASS)
data(geyser)
duration <- geyser$duration 
featureSignifGUI(duration)  ## univariate example
featureSignifGUI(geyser)    ## bivariate example 

data(earthquake)            ## trivariate example
earthquake$depth <- -log10(-earthquake$depth)
featureSignifGUI(earthquake,  scaleData=TRUE)
}
}

\keyword{smooth} 

