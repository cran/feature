featureSignifGUI <- function(x, scaleData=FALSE)
{
  fscreate.tcl <- function()
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    bw <- as.numeric(tclvalue(bw1.tcl))
    if (d>=2)
    {
      ylim <- as.numeric(c(tclvalue(ylim1.tcl), tclvalue(ylim2.tcl)))
      bw <- as.numeric(c(tclvalue(bw1.tcl), tclvalue(bw2.tcl)))
    }
    if (d>=3)
    {
      zlim <- as.numeric(c(tclvalue(zlim1.tcl), tclvalue(zlim2.tcl)))
      bw <- as.numeric(c(tclvalue(bw1.tcl), tclvalue(bw2.tcl), tclvalue(bw3.tcl)))
    }
    gs <- rep(as.numeric(tclvalue(gridsize.tcl)), d)
    fs <- featureSignif(x, bw=bw, addSignifGrad=TRUE, addSignifCurv=TRUE, gridsize=gs)
    assign("fs", fs, envir=fs.env)

    if (d==1) plot(fs, addKDE=TRUE, add=FALSE, xlim=xlim, xlab=tclvalue(xlab.tcl))
    if (d==2) plot(fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl))
    if (d==3) plot(fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl), zlab=tclvalue(zlab.tcl), addAxes3d=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    tkmessageBox(title="featureSignifGUI", message="Feature significance computations complete.", type="ok")
    return()
  }
  
  fssizer.tcl <- function()
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    xlab <- tclvalue(xlab.tcl)
    gs <- rep(as.numeric(tclvalue(gridsize.tcl)), d)
    fs.env$fs.SiZer <- SiZer(x, bw=bw.range[[1]], gridsize=gs, xlim=xlim, xlab=xlab, plotSiZer=TRUE)
    tkmessageBox(title="featureSignifGUI", message="SiZer map computed.", type="ok")
    return()
  }

  fsdata.tcl <- function()
  {
    addDataNum <- as.numeric(tclvalue(addDataNum.tcl))
    if (d>=3) dataAlpha <- as.numeric(tclvalue(dataAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d<3)
        plot(fs.env$fs, addKDE=FALSE, addData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(fs.env$fs, addKDE=FALSE, addData=TRUE, add=TRUE, dataAlpha=dataAlpha, addDataNum=addDataNum)
    return()
  }

  fsgraddata.tcl <- function()
  {
    addDataNum <- as.numeric(tclvalue(addDataNum.tcl))
    if (d>=3) gradDataAlpha <- as.numeric(tclvalue(gradDataAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d<3)
        plot(fs.env$fs, addKDE=FALSE, addSignifGradData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(fs.env$fs, addKDE=FALSE, addSignifGradData=TRUE, add=TRUE, gradDataAlpha=gradDataAlpha, addDataNum=addDataNum)
    return()
  }

  fsgradregion.tcl <- function(panel)
  {
    if (d>=3) gradRegionAlpha <- as.numeric(tclvalue(gradRegionAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d<3)
        plot(fs.env$fs, addKDE=FALSE, addSignifGradRegion=TRUE, add=TRUE)
      else
        plot(fs.env$fs, addKDE=FALSE, addSignifGradRegion=TRUE, add=TRUE, gradRegionAlpha=gradRegionAlpha, addAxes=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    return()
  }

  fscurvdata.tcl <- function(panel)
  {
    addDataNum <- as.numeric(tclvalue(addDataNum.tcl))
    if (d>=3) curvDataAlpha <- as.numeric(tclvalue(curvDataAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d <3)
        plot(fs.env$fs, addKDE=FALSE, addSignifCurvData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(fs.env$fs, addKDE=FALSE, addSignifCurvData=TRUE, add=TRUE, curvDataAlpha=curvDataAlpha, addDataNum=addDataNum)
    return()
  }

  fscurvregion.tcl <- function(panel)
  {
    if (d>=3) curvRegionAlpha <- as.numeric(tclvalue(curvRegionAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d <3)
         plot(fs.env$fs, addKDE=FALSE, addSignifCurvRegion=TRUE, add=TRUE)
      else
        plot(fs.env$fs, addKDE=FALSE, addSignifCurvRegion=TRUE, add=TRUE, curvRegionAlpha=curvRegionAlpha, addAxes=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    return(panel)
  }

  fsclearall1.tcl <- function()
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    if (d>=2) ylim <- as.numeric(c(tclvalue(ylim1.tcl), tclvalue(ylim2.tcl)))
    if (d>=3) zlim <- as.numeric(c(tclvalue(zlim1.tcl), tclvalue(zlim2.tcl)))
           
    if (!is.null(fs.env$fs))
    {
      if (d==1) plot(fs.env$fs, addKDE=TRUE, add=FALSE, xlim=xlim, xlab=tclvalue(xlab.tcl))
      if (d==2) plot(fs.env$fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl))
      if (d==3) plot(fs.env$fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl), zlab=tclvalue(zlab.tcl),  addAxes3d=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    }
    return()
  }
  fsclearall2.tcl <- function(panel)
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    if (d>=2) ylim <- as.numeric(c(tclvalue(ylim1.tcl), tclvalue(ylim2.tcl)))
    if (d>=3) zlim <- as.numeric(c(tclvalue(zlim1.tcl), tclvalue(zlim2.tcl)))
           
    if (!is.null(fs.env$fs))
    {
      if (d==1) plot(fs.env$fs, addKDE=FALSE, add=FALSE, xlim=xlim, xlab=tclvalue(xlab.tcl))
      if (d==2) plot(fs.env$fs, addKDE=FALSE, add=FALSE, xlim=xlim, ylim=ylim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl))
      if (d==3) plot(fs.env$fs, addKDE=FALSE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl), zlab=tclvalue(zlab.tcl),  addAxes3d=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    }
    return()
  }
  
  ######################################################################################################
  ## GUI
  #####################################################################################################

  fs.env <- new.env()
    
  button.col <- "lightblue"
  bg.col <- "white"
  fg.col <- "black"
  heading.col <- "blue"
  heading.font <- tkfont.create(family="helvetica", weight="bold")
  space.font <- tkfont.create(size=10)
  textwidth <- 10
  #space.label <- paste(rep(" ", scalebar.width), collapse="")
  tcl("tk_setPalette", "background", bg.col, "selectColor", "grey75") 

  ## set defaults
  
  if (is.vector(x))
  { d <- 1; n <- length(x); names.x <- deparse(substitute(x));  if (scaleData)  x <- (x-min(x))/(max(x) - min(x))}
  else
  { d <- ncol(x); n <- nrow(x); names.x <- colnames(x)
    if (is.null(names.x))
    {
       names.xx <- deparse(substitute(x))
       names.xx <- strsplit(names.xx, "\\[")[[1]][1]
       names.x <- paste(names.xx, "[,", 1:d, "]", sep="")
    }
    if (scaleData)
      for (i in 1:d)
        x[,i] <- (x[,i]-min(x[,i]))/(max(x[,i]) - min(x[,i]))
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
  
  tt <- tktoplevel(width=fspanel.width, height=fspanel.height)
  tt.name <- "featureSignif" 
  tktitle(tt) <- paste("featureSignif GUI")
  tkgrid(tklabel(tt,text="featureSignif GUI: feature significance for multivariate kernel density estimation", font=heading.font, foreground=heading.col), row=1, columnspan=7+2*(d>=3))
  tkgrid(tklabel(tt,text="", font=space.font), row=3, column=1)
  tkgrid(tklabel(tt,text="", font=space.font), row=3, column=3)
  tkgrid(tklabel(tt,text="", font=space.font), row=3, column=5)
  tkgrid(tklabel(tt,text="", font=space.font), row=3, column=7)
  if (d==3)
  {
    tkgrid(tklabel(tt,text="", font=space.font), row=3, column=7)
    tkgrid(tklabel(tt,text="", font=space.font), row=3, column=9)
    tkgrid(tklabel(tt,text=""), row=2, columnspan=11)
    tkgrid(tklabel(tt,text=""), row=4, columnspan=11)
    tkgrid(tklabel(tt,text=""), row=6, columnspan=11)  
  }
  else
  {
    tkgrid(tklabel(tt,text=""), row=2, columnspan=9)
    tkgrid(tklabel(tt,text=""), row=4, columnspan=9)
    tkgrid(tklabel(tt,text=""), row=6, columnspan=9)
  }

  ## Bandwidths
  
  col1.frame <- tkframe(tt)
  col1b.frame <- tkframe(col1.frame)
  bw1.tcl <- tclVar(h[1])
  bw2.tcl <- tclVar(h[2])
  bw3.tcl <- tclVar(h[3])

  if (d>=1)
  {
    bw1.s <- tkscale(col1.frame, from=h.low[1], to=h.upp[1], showvalue=TRUE, resolution=(h.upp[1]-h.low[1])/100, length=fspanel.button.width, variable=bw1.tcl, orient="horizontal", label=paste("Bandwidth for", names.x[1]))
    tkconfigure(bw1.s, variable=bw1.tcl) 
    tkgrid(bw1.s, sticky="we")
  }
  if (d>=2)
  {
    bw2.s <- tkscale(col1.frame, from=h.low[2], to=h.upp[2], showvalue=TRUE,  resolution=(h.upp[2]-h.low[2])/100, length=fspanel.button.width, variable=bw2.tcl, orient="horizontal", label=paste("Bandwidth for", names.x[2]))
    tkconfigure(bw2.s, variable=bw2.tcl) 
    tkgrid(bw2.s, sticky="we")
  }
  if (d>=3)
  {
    bw3.s <- tkscale(col1.frame, from=h.low[3], to=h.upp[3], showvalue=TRUE,  resolution=(h.upp[3]-h.low[3])/100, length=fspanel.button.width, variable=bw3.tcl, orient="horizontal", label=paste("Bandwidth for", names.x[3]))
    tkconfigure(bw3.s, variable=bw3.tcl) 
    tkgrid(bw3.s, sticky="we")
  }
  tkgrid(col1b.frame, sticky="w")
  tkgrid(tklabel(col1.frame, text=" "), sticky="w")
  
  ## Grid size

  col1b.frame <- tkframe(col1.frame)
  gridsize.tcl <- tclVar(gridsize)
  grid.e <- tkentry(col1b.frame, textvariable=gridsize.tcl, width=textwidth)
  tkgrid(tklabel(col1b.frame, text="Grid size "), grid.e)
  tkgrid(col1b.frame, sticky="w")
  tkgrid(col1.frame, row=3, column=2, sticky="nw")


  ## Axis limits
  col2.frame <- tkframe(tt)

  if (d>=1)
  {  
    col2b.frame <- tkframe(col2.frame)
    xlim1.tcl <- tclVar(xlim[1])
    xlim1.e <- tkentry(col2b.frame, textvariable=xlim1.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis lower limit "), xlim1.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    xlim2.tcl <- tclVar(xlim[2])
    xlim2.e <- tkentry(col2b.frame, textvariable=xlim2.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis upper limit "), xlim2.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    xlab.tcl <- tclVar(names.x[1])
    xlab.e <- tkentry(col2b.frame, textvariable=xlab.tcl)
    tkgrid(tklabel(col2b.frame, text="Axis label "), xlab.e)
    tkgrid(col2b.frame, sticky="w")
    tkgrid(tklabel(col2.frame, text=" "), sticky="w")
    tkgrid(col2.frame, row=3, column=4, sticky="ne")
    
  }
  if (d>=2)
  {  
    col2b.frame <- tkframe(col2.frame)
    ylim1.tcl <- tclVar(ylim[1])
    ylim1.e <- tkentry(col2b.frame, textvariable=ylim1.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis lower limit "), ylim1.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    ylim2.tcl <- tclVar(ylim[2])
    ylim2.e <- tkentry(col2b.frame, textvariable=ylim2.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis upper limit "), ylim2.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    ylab.tcl <- tclVar(names.x[2])
    ylab.e <- tkentry(col2b.frame, textvariable=ylab.tcl)
    tkgrid(tklabel(col2b.frame, text="Axis label "), ylab.e)
    tkgrid(col2b.frame, sticky="w")
    tkgrid(tklabel(col2.frame, text=" "), sticky="w")
  }
  if (d>=3)
  {  
    col2b.frame <- tkframe(col2.frame)
    zlim1.tcl <- tclVar(zlim[1])
    zlim1.e <- tkentry(col2b.frame, textvariable=zlim1.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis lower limit "), zlim1.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    zlim2.tcl <- tclVar(zlim[2])
    zlim2.e <- tkentry(col2b.frame, textvariable=zlim2.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis upper limit "), zlim2.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    zlab.tcl <- tclVar(names.x[3])
    zlab.e <- tkentry(col2b.frame, textvariable=zlab.tcl)
    tkgrid(tklabel(col2b.frame, text="Axis label "), zlab.e)
    tkgrid(col2b.frame, sticky="w")
    tkgrid(tklabel(col2.frame, text=" "), sticky="w")
  }
  
  col2b.frame <- tkframe(col2.frame)
  addDataNum.tcl <- tclVar(1000)
  addDataNum.e <- tkentry(col2b.frame, textvariable=addDataNum.tcl, width=textwidth)
  tkgrid(tklabel(col2b.frame, text="No. of display data points "), addDataNum.e)
  tkgrid(col2b.frame, sticky="w")
  tkgrid(col2.frame, row=3, column=4, sticky="nw")

  ## Plot options
  
  col3.frame <- tkframe(tt)
  fsdata.b <- tkbutton(col3.frame, text="Add data\npoints", command=fsdata.tcl, background=button.col)
  fsgradregion.b <- tkbutton(col3.frame, text="Add significant\ngradient region", command=fsgradregion.tcl, background=button.col)
  fsgraddata.b <- tkbutton(col3.frame, text="Add significant\ngradient data", command=fsgraddata.tcl, background=button.col)
  fscurvregion.b <- tkbutton(col3.frame, text="Add significant\ncurvature region", command=fscurvregion.tcl, background=button.col)
  fscurvdata.b <- tkbutton(col3.frame, text="Add significant\ncurvature data", command=fscurvdata.tcl, background=button.col)
  tkgrid(fsdata.b)
  tkgrid(fsgradregion.b)
  tkgrid(fsgraddata.b)
  tkgrid(fscurvregion.b)
  tkgrid(fscurvdata.b)
  if (d==3) tkgrid(col3.frame, row=3, column=8, sticky="nw")
  else tkgrid(col3.frame, row=3, column=6, sticky="nw")

  ## Transparency options

  if (d==3)
  {
    col4.frame <- tkframe(tt)

    dataAlpha.tcl <- tclVar(0.1)
    dataAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1, length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of data points")
    tkconfigure(dataAlpha.s, variable=dataAlpha.tcl) 
    tkgrid(dataAlpha.s, sticky="we")

    gradRegionAlpha.tcl <- tclVar(0.2)
    gradRegionAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1, length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of signif. grad. regions")
    tkconfigure(gradRegionAlpha.s, variable=gradRegionAlpha.tcl) 
    tkgrid(gradRegionAlpha.s, sticky="we")

    gradDataAlpha.tcl <- tclVar(0.2)
    gradDataAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1,length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of signif. grad. data")
    tkconfigure(gradDataAlpha.s, variable=gradDataAlpha.tcl) 
    tkgrid(gradDataAlpha.s, sticky="we")

    curvRegionAlpha.tcl <- tclVar(0.3)
    curvRegionAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1, length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of signif. curv. regions")
    tkconfigure(curvRegionAlpha.s, variable=curvRegionAlpha.tcl) 
    tkgrid(curvRegionAlpha.s, sticky="we")

    curvDataAlpha.tcl <- tclVar(0.3)
    curvDataAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1, length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of signif. curv. data")
    tkconfigure(curvDataAlpha.s, variable=curvDataAlpha.tcl) 
    tkgrid(curvDataAlpha.s, sticky="we")

    col4b.frame <- tkframe(col4.frame)
    addAxes3d.tcl <- tclVar(TRUE)
    addAxes3d.cb <- tkcheckbutton(col4b.frame, variable=addAxes3d.tcl)
    tkgrid(addAxes3d.cb, tklabel(col4b.frame, text="Add 3D axes "))
    tkgrid(col4b.frame, sticky="nw")
    tkgrid(col4.frame, row=3, column=6, sticky="nw")
  }

  ## Computation buttons
  
  col1.frame <- tkframe(tt)
  fscreate.b <- tkbutton(col1.frame, text="Compute significant\nfeatures", command=fscreate.tcl, background=button.col)
  tkgrid(fscreate.b)
  tkgrid(col1.frame, row=5, column=2, sticky="n")

  col2.frame <- tkframe(tt)
  fsclearall1.b <- tkbutton(col2.frame, text="Reset plot\n(except KDE)", command=fsclearall1.tcl, background=button.col)
  tkgrid(fsclearall1.b)
  tkgrid(col2.frame, row=5, column=4, sticky="n")

  col3.frame <- tkframe(tt)
  if (d==1)
  {
    fssizer.b <- tkbutton(col3.frame, text="Compute SiZer\nmap", command=fssizer.tcl, background=button.col)
    tkgrid(fssizer.b)
  }
  else
  {
    fsclearall2.b <- tkbutton(col3.frame, text="Reset plot\n", command=fsclearall2.tcl, background=button.col)
    tkgrid(fsclearall2.b)
  }
  tkgrid(col3.frame, row=5, column=6, sticky="n")
}


