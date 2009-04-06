

featureSignifGUI <- function(x, scaleData=FALSE)
{
  fscreate <- function(panel)
  {
    xlim <- as.numeric(c(panel$xlim1, panel$xlim2))
    if (d>=2) ylim <- as.numeric(c(panel$ylim1, panel$ylim2))
    if (d>=3) zlim <- as.numeric(c(panel$zlim1, panel$zlim2))
    bw <- c(panel$bw1, panel$bw2, panel$bw3)
   
    gs <- rep(as.numeric(panel$gridsize), d)
    panel$fs <- featureSignif(x, bw=bw, addSignifGrad=TRUE, addSignifCurv=TRUE, gridsize=gs)
    if (d==1) plot(panel$fs, addKDE=TRUE, add=FALSE, xlim=xlim, xlab=panel$xlab)
    if (d==2) plot(panel$fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, xlab=panel$xlab, ylab=panel$ylab)
    if (d==3) plot(panel$fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim, xlab=panel$xlab, ylab=panel$ylab, zlab=panel$zlab)
    ##assign("fs.gui", panel$fs, envir=globalenv())
    ##rp.messagebox("Feature significance computations complete and saved as `fs.gui' in the R workspace.")
    rp.messagebox("Feature significance computations complete.")
    return(panel)
  }
  
  fssizer <- function(panel)
  {
    xlim <- as.numeric(c(panel$xlim1, panel$xlim2))
    ##bw <- c(panel$bw1, panel$bw2, panel$bw3)
    xlab <- panel$xlab
    gs <- rep(as.numeric(panel$gridsize), d)
    panel$fs.SiZer <- SiZer(x, bw=bw.range[[1]], gridsize=gs, xlim=xlim, xlab=xlab, plotSiZer=TRUE)
    ##assign("fs.SiZer.gui", panel$fs.SiZer, envir=globalenv())
    ##rp.messagebox("SiZer map complete and saved as `fs.SiZer.gui' in the R workspace.")
    rp.messagebox("SiZer map complete.")
    return(panel)
  }

  fsdata <- function(panel)
  {
    addDataNum <- as.numeric(panel$addDataNum)
    if (!is.null(panel$fs))
      if (d<3)
        plot(panel$fs, addKDE=FALSE, addData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(panel$fs, addKDE=FALSE, addData=TRUE, add=TRUE, dataAlpha=panel$dataAlpha, addDataNum=addDataNum)
    return(panel)
  }

  fsgraddata <- function(panel)
  {
    addDataNum <- as.numeric(panel$addDataNum)
    if (!is.null(panel$fs))
      if (d<3)
        plot(panel$fs, addKDE=FALSE, addSignifGradData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(panel$fs, addKDE=FALSE, addSignifGradData=TRUE, add=TRUE, gradDataAlpha=panel$gradDataAlpha, addDataNum=addDataNum)
    return(panel)
  }

  fsgradregion <- function(panel)
  {
    if (!is.null(panel$fs))
      if (d<3)
        plot(panel$fs, addKDE=FALSE, addSignifGradRegion=TRUE, add=TRUE)
      else
        plot(panel$fs, addKDE=FALSE, addSignifGradRegion=TRUE, add=TRUE, gradRegionAlpha=panel$gradRegionAlpha)
    return(panel)
  }

  fscurvdata <- function(panel)
  {
    addDataNum <- as.numeric(panel$addDataNum)
    if (!is.null(panel$fs))
      if (d <3)
        plot(panel$fs, addKDE=FALSE, addSignifCurvData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(panel$fs, addKDE=FALSE, addSignifCurvData=TRUE, add=TRUE, curvDataAlpha=panel$curvDataAlpha, addDataNum=addDataNum)
    return(panel)
  }

  fscurvregion <- function(panel)
  {
    if (!is.null(panel$fs))
      if (d <3)
         plot(panel$fs, addKDE=FALSE, addSignifCurvRegion=TRUE, add=TRUE)
      else
        plot(panel$fs, addKDE=FALSE, addSignifCurvRegion=TRUE, add=TRUE, curvRegionAlpha=panel$curvRegionAlpha)
    return(panel)
  }

  fsclearall1 <- function(panel)
  {
    xlim <- as.numeric(c(panel$xlim1, panel$xlim2))
    if (d>=2) ylim <- as.numeric(c(panel$ylim1, panel$ylim2))
    if (d>=3) zlim <- as.numeric(c(panel$zlim1, panel$zlim2))
   
    if (!is.null(panel$fs))
    {
      if (d==1) plot(panel$fs, addKDE=TRUE, add=FALSE, xlim=xlim, xlab=panel$xlab)
      if (d==2) plot(panel$fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, xlab=panel$xlab, ylab=panel$ylab)
      if (d==3) plot(panel$fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim,xlab=panel$xlab, ylab=panel$ylab, zlab=panel$zlab)
    }
    return(panel)
  }
  fsclearall2 <- function(panel)
  {
    xlim <- as.numeric(c(panel$xlim1, panel$xlim2))
    if (d>=2) ylim <- as.numeric(c(panel$ylim1, panel$ylim2))
    if (d>=3) zlim <- as.numeric(c(panel$zlim1, panel$zlim2))

    if (!is.null(panel$fs))
    {
      if (d==1) plot(panel$fs, addKDE=FALSE, add=FALSE, xlim=xlim, xlab=panel$xlab)
      if (d==2) plot(panel$fs, addKDE=FALSE, add=FALSE, xlim=xlim, ylim=ylim, xlab=panel$xlab, ylab=panel$ylab)
      if (d==3) plot(panel$fs, addKDE=FALSE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim, xlab=panel$xlab, ylab=panel$ylab, zlab=panel$zlab)
    }
    return(panel)
  }

  ## Create GUI panel 
  require(rpanel)
  options(guiStyle = "Rcmdr")

  ## Set some defaults

  if (is.vector(x))
  { d <- 1; n <- length(x); names.x <- deparse(substitute(x))}
  else
  { d <- ncol(x); n <- nrow(x); names.x <- colnames(x)
    if (is.null(names.x))
    {
       names.xx <- deparse(substitute(x))
       names.xx <- strsplit(names.xx, "\\[")[[1]][1]
       names.x <- paste(names.xx, "[,", 1:d, "]", sep="")
    }
  }
  if (d>4)
    stop("Feature significance only available for 1- to 4-d data")

  tau <- 5

  if (d==1) gridsize <- 401
  if (d==2) gridsize <- 151
  if (d==3) gridsize <- 31
  if (d==4) gridsize <- 21
  bw.range <- dfltBWrange(x, tau)
  bw <- matrix(unlist(bw.range), nrow=2, byrow=FALSE)
  h.low <- bw[1,]
  h.upp <- bw[2,]
  hmix.prop <- 1/4
  h <- h.low^(hmix.prop)*h.upp^(1-hmix.prop)

  if (d==1)
  {  
    xlim <- c(min(x)-h[1],max(x)+h[1])
    dfltCounts.out <- dfltCounts(x,gridsize, apply(bw, 2, max))
    gcounts <- dfltCounts.out$counts
    range.x <- dfltCounts.out$range.x  
    dest <- drvkde(gcounts, rep(0,d), bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE)
    ylim <- c(0,1.5)*max(dest$est)
  }
  else
  {
    xlim <- c(min(x[,1])-h[1],max(x[,1])+h[1])
    ylim <- c(min(x[,2])-h[2],max(x[,2])+h[2])
    if (d>2)
      zlim <- c(min(x[,3])-h[3],max(x[,3])+h[3])
  }


  ## panel dimensions
  if (d <3)
  {  
    fspanel.width <- 740
    fspanel.height <- 410
    fspanel.button.height <- 50
    fspanel.button.width <- (fspanel.width-40)/3
  }
  else
  {
    fspanel.width <- 950
    fspanel.height <- 410
    fspanel.button.height <- 50
    fspanel.button.width <- (fspanel.width-50)/4
  }
  fspanel <- rp.control(panelname="fsgui", title="featureSignif GUI: feature significance for multivariate kernel density estimation", size=c(fspanel.width, fspanel.height))


  ## bandwidth sliders and axes limits/labels
  if (d>=1)
  {
    rp.slider(fspanel, bw1, h.low[1], h.upp[1], showvalue=TRUE, title=paste("Bandwidth for", names.x[1]), initval=h[1], pos=c(10,10,fspanel.button.width, 2*fspanel.button.height))
    rp.textentry(fspanel, xlim1, title="Axis lower limit", init=xlim[1],4, pos=c(1*fspanel.button.width+20,10,fspanel.button.width, fspanel.button.height/2))
    rp.textentry(fspanel, xlim2, title="Axis upper limit", init=xlim[2],4, pos=c(1*fspanel.button.width+20, 0.5*fspanel.button.height+10,fspanel.button.width, fspanel.button.height/2))
    rp.textentry(fspanel, xlab, title="Axis label", init=names.x[1], pos=c(1*fspanel.button.width+20, 1*fspanel.button.height+10,fspanel.button.width, 0.5*fspanel.button.height))
  }
  
  if (d>=2)
  {
    rp.slider(fspanel, bw2, h.low[2], h.upp[2], showvalue=TRUE, title=paste("Bandwidth for", names.x[2]), initval=h[2], pos=c(10,2*fspanel.button.height+10,fspanel.button.width, 2*fspanel.button.height))
    rp.textentry(fspanel, ylim1, title="Axis lower limit", init=ylim[1],4, pos=c(1*fspanel.button.width+20,2*fspanel.button.height+10,fspanel.button.width, 0.5*fspanel.button.height))
    rp.textentry(fspanel, ylim2, title="Axis upper limit", init=ylim[2],4, pos=c(1*fspanel.button.width+20, 2.5*fspanel.button.height+10,fspanel.button.width, 0.5*fspanel.button.height))
    rp.textentry(fspanel, ylab, title="Axis label", init=names.x[2], pos=c(1*fspanel.button.width+20, 3*fspanel.button.height+10,fspanel.button.width, 0.5*fspanel.button.height))
  }
  if (d>=3)
  {
    rp.slider(fspanel, bw3, h.low[3], h.upp[3], showvalue=TRUE, title=paste("Bandwidth for", names.x[3]), initval=h[3], pos=c(10,4*fspanel.button.height+10,fspanel.button.width, 2*fspanel.button.height))
 
    rp.textentry(fspanel, zlim1, title="Axis lower limit", init=zlim[1],4, pos=c(1*fspanel.button.width+20,4*fspanel.button.height+10,fspanel.button.width, 0.5*fspanel.button.height))
    rp.textentry(fspanel, zlim2, title="Axis upper limit", init=zlim[2],4, pos=c(1*fspanel.button.width+20, 4.5*fspanel.button.height+10,fspanel.button.width, 0.5*fspanel.button.height))
    rp.textentry(fspanel, zlab, title="Axis label", init=names.x[3], pos=c(1*fspanel.button.width+20, 5*fspanel.button.height+10,fspanel.button.width, 0.5*fspanel.button.height))
  }

  ## Gridsize and data number text fields

  rp.textentry(fspanel, gridsize, title=paste("Grid size"), initval=gridsize, pos=c(10,6*fspanel.button.height,fspanel.button.width, fspanel.button.height-10))

  rp.textentry(fspanel, addDataNum, title="No. of display data points", init=1000, pos=c(1*fspanel.button.width+20,6*fspanel.button.height,fspanel.button.width, fspanel.button.height-10))
  
  
  ## Computation buttons
  
  rp.button(fspanel, title="Compute significant features", action=fscreate, pos=c(10,7*fspanel.button.height+10,fspanel.button.width, fspanel.button.height-10))

  rp.button(fspanel, clearAll, title="Reset plot (except KDE)", action=fsclearall1, pos=c(1*fspanel.button.width+20,7*fspanel.button.height+10,fspanel.button.width, fspanel.button.height-10))
  
  if (d==1)  rp.button(fspanel, fsSiZer, title="Compute SiZer map", action=fssizer, pos=c(2*fspanel.button.width+30,7*fspanel.button.height+10,fspanel.button.width, fspanel.button.height-10))
  else
    rp.button(fspanel, clearAll, title="Reset plot", action=fsclearall2, pos=c(2*fspanel.button.width+30,7*fspanel.button.height+10,fspanel.button.width, fspanel.button.height-10))

  ## Buttons for feature signif. plots options
  
  rp.button(fspanel, addData, title="Add data points", action=fsdata, pos=c(2*fspanel.button.width+30,10,fspanel.button.width, fspanel.button.height-10))
  rp.button(fspanel, addSignifGradRegion, title="Add significant gradient regions", action=fsgradregion, pos=c(2*fspanel.button.width+30,1*fspanel.button.height+10,fspanel.button.width, fspanel.button.height-10))
  rp.button(fspanel, addSignifGradData, title="Add significant gradient data", action=fsgraddata, pos=c(2*fspanel.button.width+30,2*fspanel.button.height+10,fspanel.button.width, fspanel.button.height-10))
  rp.button(fspanel, addSignifCurvRegion, title="Add significant curvature regions", action=fscurvregion, pos=c(2*fspanel.button.width+30,3*fspanel.button.height+10,fspanel.button.width, fspanel.button.height-10))
  rp.button(fspanel, addSignifCurvData, title="Add significant curvature data", action=fscurvdata, pos=c(2*fspanel.button.width+30,4*fspanel.button.height+10,fspanel.button.width, fspanel.button.height-10))
  if (d==3)
  {
    rp.slider(fspanel, dataAlpha, 0, 1, title="Transparency of data points", initval=0.1, pos=c(3*fspanel.button.width+40,10,fspanel.button.width, fspanel.button.height))
    rp.slider(fspanel, gradRegionAlpha, 0, 1, title="Transparency of signif. grad regions", initval=0.2, pos=c(3*fspanel.button.width+40,fspanel.button.height+10,fspanel.button.width, fspanel.button.height))
    rp.slider(fspanel, gradDataAlpha, 0, 1, title="Transparency of signif. grad data", initval=0.3, pos=c(3*fspanel.button.width+40,2*fspanel.button.height+10,fspanel.button.width, fspanel.button.height))
    rp.slider(fspanel, curvRegionAlpha, 0, 1, title="Transparency of signif. curv regions", initval=0.3, pos=c(3*fspanel.button.width+40,3*fspanel.button.height+10,fspanel.button.width, fspanel.button.height))
    rp.slider(fspanel, curvDataAlpha, 0, 1, title="Transparency of signif. curv data", initval=0.3, pos=c(3*fspanel.button.width+40,4*fspanel.button.height+10,fspanel.button.width, fspanel.button.height))
  }
  ## to prevent 'no visible binding for global variable' errors in R CMD check 
  bw1 <- NULL; bw2 <- NULL; bw3 <- NULL;
  xlim1 <- NULL; xlim2 <- NULL; xlim3 <- NULL; xlab <- NULL;
  ylim1 <- NULL; ylim2 <- NULL; ylim3 <- NULL; ylab <- NULL;
  zlim1 <- NULL; zlim2 <- NULL; zlim3 <- NULL; zlab <- NULL;
  clearAll <- NULL; fsSiZer <- NULL;
  addSignifGradRegion <- NULL; addSignifGradData <- NULL;
  addSignifCurvRegion <- NULL; addSignifCurvData <- NULL;
  gradRegionAlpha <- NULL; gradDataAlpha  <- NULL;
  curvRegionAlpha <- NULL; curvDataAlpha  <- NULL;
  addData <- NULL; addDataNum <- NULL; dataAlpha <- NULL
}
