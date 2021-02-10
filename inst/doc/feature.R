## ---- echo=FALSE, message=FALSE-----------------------------------------------
knitr::opts_chunk$set(global.par=TRUE, collapse=TRUE, comment="#>", fig.width=5, fig.height=5, fig.align="center")
options(tibble.print_min=4L, tibble.print_max=4L)

## ---- fig.width=7-------------------------------------------------------------
library(feature)
data(earthquake)
eq3 <- log10(-earthquake[,3])
eq3.fs <- featureSignif(eq3, bw=0.1)
plot(eq3.fs, xlab="-log(-depth)", addSignifGradRegion=TRUE, addData=TRUE)
xlim <- par()$usr[1:2]  ## save x-axis limits to align following SiZer plot

## ---- fig.width=7-------------------------------------------------------------
eq3.SiZer <- SiZer(eq3, xlim=xlim, bw=c(0.05, 0.5), xlab="-log(-depth)")
abline(h=log(0.1))

## -----------------------------------------------------------------------------
library(MASS)
data(geyser)
geyser.fs <- featureSignif(geyser, bw=c(4.5, 0.37))
plot(geyser.fs, addSignifCurvRegion=TRUE)

## -----------------------------------------------------------------------------
plot(geyser.fs, addSignifCurvData=TRUE)

## -----------------------------------------------------------------------------
data(earthquake)
earthquake[,3] <- -log10(-earthquake[,3])
earthquake.fs <- featureSignif(earthquake, scaleData=TRUE, bw=c(0.06, 0.06, 0.05))
plot(earthquake.fs, addKDE=FALSE, addSignifCurvRegion=TRUE)

## -----------------------------------------------------------------------------
names(earthquake.fs)

