########## R function: featureSignif ##########

# For dynamic kernel density estimation; with
# feature significance options.

# Last changed: 19 JUL 2006
  
featureSignif <-
  function(x, bw, xlab, ylab, zlab, xlim, ylim, zlim,
           addData=FALSE, scaleData=FALSE, addDataNum=1000,
           addKDE=TRUE, jitterRug=TRUE, signifLevel=0.05, 
           addSignifGradRegion=FALSE, addSignifGradData=FALSE,
           addSignifCurvRegion=FALSE, addSignifCurvData=FALSE,
           plotSiZer=FALSE, logbwSiZer=TRUE, addAxes3d=TRUE, 
           densCol, dataCol="black", gradCol="green4", curvCol="blue",
           axisCol="black", bgCol="white",
           dataAlpha=0.1, gradDataAlpha=0.3,
           gradRegionAlpha=0.2, curvDataAlpha=0.3, curvRegionAlpha=0.3,
           gridsize, gridsizeSiZer)
                    
{
  options(locatorBell=FALSE) 
  names.x <- names(x)

  ## Determine appropriate value of parameter
  ## denoted by "tau" (effective kernel support).
  
  if (addSignifGradRegion | addSignifCurvRegion |
      addSignifGradData | addSignifCurvData)
    tau <- 5
  else
    tau <- 4
  
  ## Determine the dimension of the data

  x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)

  if (d>4)
    stop("Feature significance currently only available for 1- to 4-dimensional data")
  if (scaleData)
    for (i in 1:d)
      x[,i] <- (x[,i]-min(x[,i]))/(max(x[,i]) - min(x[,i]))
  
  if (missing(gridsize)) 
  {
    if (d==1) gridsize <- 401
    if (d==2) gridsize <- rep(151,2)
    if (d==3) gridsize <- rep(51,3)
    if (d==4) gridsize <- rep(21,4)
  }
  if (missing(gridsizeSiZer))
    gridsizeSiZer <- 101
  
  ## Set some defaults

  if (missing(bw))           ## b/w not specified -> interactive 
  {
    bw.range <- dfltBWrange(x,gridsize,tau)
    bw <- matrix(unlist(bw.range), nrow=2, byrow=FALSE)
    dfltCounts.out <- dfltCounts(x,gridsize, apply(bw, 2, max))
    interactive <- TRUE 
  }
  else
    if (d==1)
    {
      if ((length(bw)==1))   ## scalar b/w -> non-interactive
      {
        dfltCounts.out <- dfltCounts(x,gridsize, bw)
        interactive <- FALSE
      }
      else if (length(bw)>1) ## vector of b/w -> interactive
      {
        bw <- matrix(bw, nrow=2, byrow=FALSE)
        dfltCounts.out <- dfltCounts(x,gridsize, apply(bw, 2, max))
        interactive <- TRUE 
      }
    } 
    else if (d > 1)
    {
      if (is.vector(bw))     ## vector of b/w -> non-interactive
      {
        dfltCounts.out <- dfltCounts(x,gridsize, bw)
        interactive <- FALSE
      }
      else                   ## matrix of b/w -> interactive
      {
        dfltCounts.out <- dfltCounts(x,gridsize, apply(bw, 2, max))
        interactive <- TRUE 
      }
    }
       
  ## defaults for interactive plots
  if (interactive)
  {
    h.low <- bw[1,]
    h.upp <- bw[2,]
    h.init <- sqrt(h.low*h.upp)
    logh.low <- logb(h.low,10)
    logh.upp <- logb(h.upp,10)
    h <- h.init
  }
  else
    h <- bw

  gcounts <- dfltCounts.out$counts
  range.x <- dfltCounts.out$range.x  

  dest <- drvkde(gcounts, rep(0,d), bandwidth=h, binned=TRUE,
                 range.x=range.x, se=FALSE)
  
  ## adjust values of density estimate to be positive
  dest$est[dest$est<0] <- 0 #1e-6*min(abs(dest$est[dest$est>0]))
  
  ## random sample of data points used for display
  n <- nrow(x)
  nsamp <- min(addDataNum, n)
  
  if (nsamp < n)
  {
    rand.inds <- sort(sample(1:n, nsamp, replace=FALSE))
    x.rand <- as.matrix(x[rand.inds,])
  }
  else
    x.rand <- x
  
  ## Determine default xlim, ylim, zlim
  if (missing(xlim))
     if (d==1)
       xlim <- c(min(x)-h[1],max(x)+h[1])
     else
       xlim <- c(min(x[,1])-h[1],max(x[,1])+h[1])

  if (missing(ylim))
    if (d==1)
      ylim <- c(0,1.5)*max(dest$est)
    else if (d>1)
      ylim <- c(min(x[,2])-h[2],max(x[,2])+h[2])
  
  if (missing(zlim) & d>2)
    zlim <- c(min(x[,3])-h[3],max(x[,3])+h[3])

  if (d==1)
    lims <- list(xlim)
  if (d==2)
    lims <- list(xlim, ylim)
  if (d==3)
    lims <- list(xlim, ylim, zlim)
  if (d==4)
    lims <- list(xlim, ylim, zlim, c(min(x[,4])-h[4],max(x[,4])+h[4]))

  ## Determine default axis labels.

  if (missing(xlab)) xlab <- NULL
  if (missing(ylab)) ylab <- NULL
  if (missing(zlab)) zlab <- NULL
  labs <- dfltLabs(d,names.x,xlab,ylab,zlab)
  xlab <- labs$xlab ; ylab <- labs$ylab ; zlab <- labs$zlab
  
  ## Determine plotting indices.
  
  plot.inds <- list()
  for (id in 1:d)
  {
    plot.inds.l <- (1:gridsize[id])[dest$x.grid[[id]]>=lims[[id]][1]]
    plot.inds.u <- (1:gridsize[id])[dest$x.grid[[id]]<=lims[[id]][2]]
    plot.inds[[id]] <- intersect(plot.inds.l,plot.inds.u)
  }


  ## Do initial plot.

  if (missing(densCol))
    if (d==1)
      densCol <- "DarkOrange"
    else if (d==2)
      densCol <- rev(heat.colors(1000))
    else if (d==3)
      densCol <- rev(heat.colors(5))
 
  if (d<3)
  {
  #  op <- par(plt=par("plt"),xpd=par("xpd"),bg=par("bg"))
  #  plt.dflt <- par("plt"); plt.new <-  plt.dflt
  #  plt.new[3] <- 2.5*plt.new[3] ; par(plt=plt.new)
  #  par(xpd=TRUE)
    par(bg=bgCol)
    old.mar <- c(5.1, 4.1, 4.1, 4.1) 
  }
  if (d >=3)
  {
    open3d()
    clear3d()
    screen1 <- rgl.cur()
    if (interactive)
      rgl.viewpoint(theta=0, phi=0)
    else
      rgl.viewpoint(theta=0, phi=-90)
    rgl.bg(col=bgCol)
    pop3d(type="lights")
    light3d(theta=0, phi=30)
    material3d(alpha=1)
    material3d(back="fill")
  }
  if (d==1)
  {
    clk.bar.msg1 <- "Click in: green bar to vary bandwidth or red box to stop"
    clk.bar.msg2 <- "          green bar to vary bandwidth                   "
    clk.bar.msg3 <- "                                         red box to stop"
  }
  if (d==2)
  {
    clk.bar.msg1 <- "Click in: green bar to vary bandwidths or red box to stop"
    clk.bar.msg2 <- "          green bar to vary bandwidths                   "
    clk.bar.msg3 <- "                                          red box to stop"
  }
    
  
  if (plotSiZer)
    interactive <- FALSE

  if (d==1)
    bw.text <- paste("bandwidth =", toString(signif(h,3)))
  else if (d>1)
    bw.text <- paste("bandwidths = (", toString(signif(h,3)), ")", sep="")
  
  
  if (interactive)
  {
    if (d <3)
    {
      split.screen(matrix(c(0,1,0,0.2, 0,1,0.2,1), byrow=TRUE, ncol=4))
      screen(1, new=TRUE)
      par(mar=c(1,1,1,1))
      plot(0, xlim=c(0,1), ylim=c(0,1), type="n", axes=FALSE)

      if (d==2)
      {  
        ## draw wait symbol
        wait.x <- 0.1
        wait.y <- 0.4
        wait.r <- 0.05 
        symbols(x=wait.x, y=wait.y, circle=wait.r, add=TRUE, inches=FALSE,
              fg="blue", lwd=3)
        text(wait.x, wait.y, "WAIT", col="blue")
      }
      
      ## draw bandwidth bar 
      bw.bar <- c(0.2, 0.8, 0.3, 0.5) ## (x1, x2, y1, y2)
      rect(bw.bar[1], bw.bar[3], bw.bar[2], bw.bar[4], col="green", border=FALSE) 

      ## draw stop box
      stop.box <- c(0.9, 0.95, 0.3, 0.5) 
      rect(stop.box[1],stop.box[3], stop.box[2], stop.box[4], col="red",
           border=FALSE) 
      text((stop.box[1]+stop.box[2])/2, 0.15, "STOP", col="red")
      
      ## draw bandwidths notch
      sca.fac <- (logh.upp-logh.low)/(bw.bar[2]-bw.bar[1])
      x.ck <- (log(h, 10) - logh.low)/sca.fac + bw.bar[1]
      notch.y <- c(0.90*bw.bar[3]+0.10*bw.bar[4],0.10*bw.bar[3]+0.90*bw.bar[4])
      lines(rep(x.ck[1],2), notch.y,lwd=3,col="darkgreen")

      ## draw bandwidths
      if (addKDE | addSignifGradRegion | addSignifGradData |
          addSignifCurvRegion | addSignifCurvData)
        text((bw.bar[1]+bw.bar[2])/2, bw.bar[4]+0.2, bw.text,col="purple4")
    }
    else if (d >=3)
    {
      bw.bar <- c(0.2, 0.8, 0.4, 0.5) ## (x1, x2, y1, y2)
      stop.box <- c(1, 1.1, 0.4, 0.5) 
      wait.x <- 0

      ## draw bandwidth bar
      quads3d(c(bw.bar[1], bw.bar[1], bw.bar[2], bw.bar[2]),
              c(bw.bar[3], bw.bar[4], bw.bar[4], bw.bar[3]),c(0,0,0,0),
              col="green")
      
      ## draw stop box
      quads3d(c(stop.box[1], stop.box[1], stop.box[2], stop.box[2]),
                   c(stop.box[3], stop.box[4], stop.box[4], stop.box[3]),
                   c(0,0,0,0), col="red")
      texts3d(c(stop.box[1]+stop.box[2])/2,stop.box[3]-0.1, 0,
                   "STOP", col="red", adj=0.3)
      
      ## draw bandwidth notch
      sca.fac <- (logh.upp-logh.low)/(bw.bar[2]-bw.bar[1])
      x.ck <- (log(h, 10) - logh.low)/sca.fac + bw.bar[1]
      lines3d(c(x.ck[1], x.ck[1]),
              c(0.9*bw.bar[3]+0.1*bw.bar[4], 0.1*bw.bar[3]+0.9*bw.bar[4]),
              c(0.05,0.05), col="darkgreen", size=3)

      ## draw bandwidths 
      
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.2, 0,
                   "bandwidths = ", col="purple4", adj=0.3)
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.1, 0,  
                   paste("(", toString(signif(h,3)), ")", sep=""),
                   col="purple4", adj=0.3)

      ## draw wait symbol
      spheres3d(wait.x, (bw.bar[3]+bw.bar[4])/2, 0, r=0.05,
                     col="blue", alpha=0.5)
      texts3d(wait.x,bw.bar[3]-0.1,0, "WAIT", col="blue", adj=0.3) 
    }
  }
  
  ## draw main plot
  
  if (d>=3)
  {
    if (interactive)
    {
       open3d()
       screen2 <- rgl.cur() 
      clear3d()
      rgl.viewpoint(theta=0, phi=-90)
      rgl.bg(col=bgCol)
    }
  }
 
    
  if (d==1)
  {
    if (!plotSiZer)
    {
      if (interactive) screen(2)
      par(mar=old.mar)
      plot(dest$x.grid[[1]][plot.inds[[1]]], dest$est[plot.inds[[1]]],
           type="n",bty="l" ,col=densCol, lwd=2, xlim=xlim, ylim=ylim,
           xlab=xlab,ylab="kernel density estimate")
      
      lines(dest$x.grid[[1]][plot.inds[[1]]],dest$est[plot.inds[[1]]],
            bty="l",col=densCol,lwd=2)
    }
  }
  
  if (d==2)
  {
     x.grid.1 <- dest$x.grid[[1]] ; x.grid.2 <- dest$x.grid[[2]]

     if (interactive) screen(2)
     par(mar=old.mar)
     image(x.grid.1[plot.inds[[1]]],x.grid.2[plot.inds[[2]]],
           dest$est[plot.inds[[1]],plot.inds[[2]]],col=densCol,
           xlim=xlim, ylim=ylim, xlab=xlab,ylab=ylab,bty="n")
     
     ## Add a border around the image

     x1.bor.low <- min(x.grid.1[plot.inds[[1]]]) 
     x1.bor.upp <- max(x.grid.1[plot.inds[[1]]]) 
     x2.bor.low <- min(x.grid.2[plot.inds[[2]]]) 
     x2.bor.upp <- max(x.grid.2[plot.inds[[2]]]) 

     lines(c(x1.bor.low,x1.bor.upp),rep(x2.bor.low,2),lwd=2,col="black")
     lines(c(x1.bor.low,x1.bor.upp),rep(x2.bor.upp,2),lwd=2,col="black")
     lines(rep(x1.bor.low,2),c(x2.bor.low,x2.bor.upp),lwd=2,col="black")
     lines(rep(x1.bor.upp,2),c(x2.bor.low,x2.bor.upp),lwd=2,col="black")
   }

  if (d==3) 
  {
    x.gd.1 <- dest$x.grid[[1]] ; x.gd.2 <- dest$x.grid[[2]]
    x.gd.3 <- dest$x.grid[[3]]
    num.levs <- length(densCol)
    
    if (addKDE)
    { 
      alph <- seq(0.1,0.5,length=num.levs)
      lev.vals <- seq(0, max(dest$est), length=num.levs+2)[-c(1, num.levs+2)]

      for (il in 1:num.levs)
        contour3d(dest$est,level=lev.vals[il],
                  x=x.gd.1,y=x.gd.2,z=x.gd.3,color=densCol[il],alpha=alph[il],
                  add=(il!=1))   
    }
  }
  
  if (d==4) 
  { 
    x.gd.1 <- dest$x.grid[[1]] ; x.gd.2 <- dest$x.grid[[2]]
    x.gd.3 <- dest$x.grid[[3]] ; x.gd.4 <- dest$x.grid[[4]]

    addKDE <- FALSE
    if (addSignifGradRegion)
      cat("\nDisplay of signif. gradient regions for 4-d data not available\n")
    addSignifGradRegion <- FALSE
  } 

  ## significant features sub-function

  addSignifFeature <- function(h, dest)
  {
    if (d>=3 & addAxes3d)
    {
      lines3d(xlim, rep(ylim[1],2), rep(zlim[1],2), size=3, color=axisCol, alpha=1)
      lines3d(rep(xlim[1],2), ylim, rep(zlim[1],2), size=3, color=axisCol, alpha=1)
      lines3d(rep(xlim[1],2), rep(ylim[1],2), zlim, size=3, color=axisCol, alpha=1)
  
      texts3d(xlim[2],ylim[1],zlim[1],xlab,size=3, color=axisCol, adj=0, alpha=1)
      texts3d(xlim[1],ylim[2],zlim[1],ylab,size=3, color=axisCol, adj=1, alpha=1)
      texts3d(xlim[1],ylim[1],zlim[2],zlab,size=3, color=axisCol, adj=1, alpha=1)

      ## something wrong with the stack of rgl objects - objects being deleted
      ## when swapping between rgl windows
      ## so need to repeat plotting of axes and axes labels
      lines3d(xlim, rep(ylim[1],2), rep(zlim[1],2), size=3, color=axisCol, alpha=1)
      lines3d(rep(xlim[1],2), ylim, rep(zlim[1],2), size=3, color=axisCol, alpha=1)
      lines3d(rep(xlim[1],2), rep(ylim[1],2), zlim, size=3, color=axisCol, alpha=1)
      texts3d(xlim[2],ylim[1],zlim[1],xlab,size=3, color=axisCol, adj=0, alpha=1)
      texts3d(xlim[1],ylim[2],zlim[1],ylab,size=3, color=axisCol, adj=1, alpha=1)
      texts3d(xlim[1],ylim[1],zlim[2],zlab,size=3, color=axisCol, adj=1, alpha=1)

      lines3d(xlim, rep(ylim[1],2), rep(zlim[1],2), size=3, color=axisCol, alpha=1)
      lines3d(rep(xlim[1],2), ylim, rep(zlim[1],2), size=3, color=axisCol, alpha=1)
      lines3d(rep(xlim[1],2), rep(ylim[1],2), zlim, size=3, color=axisCol, alpha=1)
      texts3d(xlim[2],ylim[1],zlim[1],xlab,size=3, color=axisCol, adj=0, alpha=1)
      texts3d(xlim[1],ylim[2],zlim[1],ylab,size=3, color=axisCol, adj=1, alpha=1)
      texts3d(xlim[1],ylim[1],zlim[2],zlab,size=3, color=axisCol, adj=1, alpha=1)
    }

    SignifFeatureRegion.mat <-
      SignifFeatureRegion(n,d,gcounts,gridsize,dest,
                          h,signifLevel,range.x,
                          grad=(addSignifGradRegion | addSignifGradData),
                          curv=(addSignifCurvRegion | addSignifCurvData))

    ESS <- n*dest$est*prod(h)*(sqrt(2*pi)^d)
    SigESS <- ESS >= 5

    if (addSignifGradRegion | addSignifGradData)
      SignifGradRegion.mat <- SignifFeatureRegion.mat$grad
    
    if (addSignifGradData)  
      SignifGradData.mat <- SignifFeatureData(x.rand, d, dest,SignifGradRegion.mat)
    
    if (addSignifCurvRegion | addSignifCurvData)
      SignifCurvRegion.mat <- SignifFeatureRegion.mat$curv
    
    if (addSignifCurvData)
      SignifCurvData.mat <- SignifFeatureData(x.rand, d, dest,SignifCurvRegion.mat)
    
    if (addSignifGradRegion)
      if (d<3)  
        addSignifFeatureRegion(d,gridsize,SignifGradRegion.mat,plot.inds,gradCol,
                               dest,lims)
      else if (d==3)
        addSignifFeatureRegion(d,gridsize,SignifGradRegion.mat,plot.inds,gradCol,
                              dest,lims,trans.alpha=gradRegionAlpha)
                               
    if (addSignifCurvRegion)
      if (d<3)
        addSignifFeatureRegion(d,gridsize,SignifCurvRegion.mat,plot.inds,
                               curvCol,dest,lims)
      else if (d==3)
        addSignifFeatureRegion(d,gridsize,SignifCurvRegion.mat,plot.inds,curvCol,
                               dest,lims,trans.alpha=curvRegionAlpha)
      else if (d==4)
        addSignifFeatureRegion(d,gridsize,SignifCurvRegion.mat,plot.inds,curvCol,
                               dest,lims,trans.alpha=c(0.1,0.4))
   
    if (addSignifGradData)
      addSignifFeatureData(x.rand,SignifGradData.mat,gradCol, trans.alpha=gradDataAlpha)
    
    if (addSignifCurvData)
      addSignifFeatureData(x.rand,SignifCurvData.mat,curvCol, trans.alpha=curvDataAlpha)

    if (addData)
      if (d==1)
      {
        if (jitterRug)
          x.rug <- jitter(x.rand)
        else
          x.rug <- x.rand
        rug(x.rug)
      }
      else if (d==2)
        points(x.rand, col=dataCol)
      else if (d>=3)
        points3d(x.rand[,1],x.rand[,2],x.rand[,3],size=3,col=dataCol, alpha=dataAlpha)

    
    if (!addSignifGradData & ! addSignifCurvData)
      return (SignifFeatureRegion.mat)
    else if (!addSignifGradData & addSignifCurvData)
      return (c(SignifFeatureRegion.mat, list(curvData=SignifCurvData.mat)))
    else if (addSignifGradData & !addSignifCurvData)
      return (c(SignifFeatureRegion.mat, list(gradData=SignifGradData.mat)))
    else if (addSignifGradData & addSignifCurvData)
      return (c(SignifFeatureRegion.mat, list(gradData=SignifGradData.mat, curvData=SignifCurvData.mat)))
  }

  if (!plotSiZer)                  ## draw feature significance plot
    feat <- addSignifFeature(h=h, dest=dest)
  else                             ## draw SiZer plot
  {
    gs.SiZer <- gridsizeSiZer
    ##if ((length(bw)==1))   ## scalar b/w -> non-interactive
    ##{
      bw.range.SiZer  <- dfltBWrange(x,gridsize=gs.SiZer,tau)
      bw.SiZer <- matrix(unlist(bw.range.SiZer), nrow=2, byrow=FALSE)
    ##}
    ##else
    ##  bw.SiZer <- bw

    
    dfltCounts.out.SiZer  <- dfltCounts(x,gridsize=gs.SiZer, apply(bw.SiZer, 2, max))
    range.x.SiZer <-dfltCounts.out.SiZer$range.x
    gcounts.SiZer <- dfltCounts.out.SiZer$counts
    x.SiZer  <- seq(range.x.SiZer[[1]][1], range.x.SiZer[[1]][2], length=gs.SiZer) 
    bw.SiZer  <- seq(log(bw.SiZer[1,1]), log(bw.SiZer[2,1]), length=101)
    SiZer.map <- matrix(0, ncol=length(bw.SiZer), nrow=length(x.SiZer))

    i <- 0
    for (logh in bw.SiZer) 
    {
      h <- exp(logh)
      i <- i + 1
      
      est.dens <- drvkde(gcounts.SiZer,drv=0,bandwidth=h, binned=TRUE,
                         range.x=range.x.SiZer, se=FALSE)
      est.dens$est[est.dens$est<0] <- 0
      ESS <- n*est.dens$est*prod(h)*(sqrt(2*pi)^d)
      sig.ESS <- ESS >= 5
      
      sig.grad <- SignifFeatureRegion(n,d,gcounts.SiZer,gridsize=gs.SiZer,
                                      est.dens, h,signifLevel,
                                      range.x.SiZer, grad=TRUE, curv=FALSE)$grad

      est.grad <- drvkde(gcounts.SiZer, drv=1, bandwidth=h, binned=TRUE,
                         range.x=range.x.SiZer, se=FALSE)$est
          
      ## Gradient SiZer map colours
      ## 0 = grey   = sparse data     
      ## 1 = purple = zero grad
      ## 2 = blue   = +ve grad
      ## 3 = red    = -ve grad
      SiZer.col <- rep(0, length(ESS))
      SiZer.col[sig.ESS] <- 1
      SiZer.col[sig.ESS & sig.grad & est.grad >0] <- 2 
      SiZer.col[sig.ESS & sig.grad & est.grad <0] <- 3
      SiZer.map[,i] <- SiZer.col
    }
   
    if (logbwSiZer)
      image(x.SiZer, bw.SiZer, SiZer.map, breaks=c(-1,0,1,2,3),
            col=c("grey", "purple", "blue", "red"), ylab="log(bandwidth)", xlab=xlab,
            xlim=xlim)
    else
      image(x.SiZer, exp(bw.SiZer), SiZer.map, breaks=c(-1,0,1,2,3),
            col=c("grey", "purple", "blue", "red"), ylab="bandwidth", xlab=xlab,
            xlim=xlim)
  }

  if (d < 3)
    if (interactive)
    {
      screen(1, new=FALSE)
      plot(0, xlim=c(0,1), ylim=c(0,1), type="n", axes=FALSE)
      if (d==2)
      {  
        ## remove wait symbol
        symbols(x=wait.x, y=wait.y, circle=1.05*wait.r, add=TRUE, inches=FALSE,
                fg=bgCol, bg=bgCol, lwd=3)
      }
      
      ## draw bandwidth bar message
      text((bw.bar[1]+bw.bar[2])/2, 0.15, clk.bar.msg1, col="black")
      #text((bw.bar[1]+bw.bar[2])/2, 0.15, clk.bar.msg2, col="dark green")
      #text((bw.bar[1]+bw.bar[2])/2, 0.15, clk.bar.msg3, col="red")
    }
    else
      if ((addKDE | addSignifGradRegion | addSignifGradData |
          addSignifCurvRegion | addSignifCurvData) & !plotSiZer)
        title(sub=bw.text,col.sub="purple4")

  if (d>=3) 
  {
    
    if (interactive)
    {
      rgl.set(screen1)

      ## something wrong with stack of RGL objects when swapping between RGL
      ## windows - instead of just pop3d-ing the wait symbol, have to redraw
      ## all objects
      
      ## erase wait symbol
      ## pop3d();pop3d()

      clear3d(type="shapes")
      material3d(alpha=1)
      
      ## draw bandwidth bar
      quads3d(c(bw.bar[1], bw.bar[1], bw.bar[2], bw.bar[2]),
              c(bw.bar[3], bw.bar[4], bw.bar[4], bw.bar[3]),c(0,0,0,0),
              col="green")
      
      ## draw stop box
      quads3d(c(stop.box[1], stop.box[1], stop.box[2], stop.box[2]),
              c(stop.box[3], stop.box[4], stop.box[4], stop.box[3]),
              c(0,0,0,0), col="red")
      texts3d(c(stop.box[1]+stop.box[2])/2,stop.box[3]-0.1, 0,
                "STOP", col="red", adj=0.3)

      ## draw bandwidths notch
      x.bw.ck <- (log(h, 10) - logh.low)/sca.fac + bw.bar[1]    
      lines3d(c(x.bw.ck[1], x.bw.ck[1]),
              c(0.9*bw.bar[3]+0.1*bw.bar[4], 0.1*bw.bar[3]+0.9*bw.bar[4]),
              c(0.05,0.05), color="darkgreen", size=3)
      
      ## draw bandwidths 
      bw.text <- paste("bandwidths = (", toString(signif(h,3)), ")", sep="")
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.2, 0,
                   "bandwidths = ", col="purple4", adj=0.3)
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.1, 0,  
                   paste("(", toString(signif(h,3)), ")", sep=""),
                   col="purple4", adj=0.3)
     
      ## draw bandwidth bar message
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[3]-0.15, 0,
              "Select rectangular region",col="black",adj=0.3)
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[3]-0.25, 0,
              "using mouse button in:",col="black",adj=0.3)
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[3]-0.35, 0,
              "green bar to vary bandwidths", col="darkgreen", adj=0.3)
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[3]-0.45, 0,
              "red box to stop", col="red", adj=0.3)
      
     
      
    }
    else if (addKDE | addSignifGradRegion | addSignifGradData |
        addSignifCurvRegion | addSignifCurvData)
      texts3d(xlim[2], ylim[2], zlim[2], bw.text, color="purple4", alpha=1)
  }


  #################################################################
  ## Facilitate dynamic kernel density estimation.
  #################################################################
 
  finished <- FALSE
  if (plotSiZer) finished <- TRUE
  legal <- TRUE
  prev.dest <- dest  
  prev.click.vec <- NULL; 
  
  while (!finished & interactive & d < 3)
  {
    prev.x.ck <- x.ck 
    sca.fac <- (logh.upp-logh.low)/(bw.bar[2]-bw.bar[1])
    screen(1, new=FALSE)

    click.vec <- locator(1)
    
    x.ck <- click.vec$x
    y.ck <- click.vec$y
    
    ## Check for finish
    
    if ((x.ck>=stop.box[1])&(x.ck<=stop.box[2]) &
        (y.ck>=stop.box[3])&(y.ck<=stop.box[4]))
      finished <- TRUE

    
    if (!finished)
    {
      ## Check for legality
      legal <- (x.ck>=bw.bar[1]) & (x.ck<=bw.bar[2]) & 
               (y.ck>=bw.bar[3]) & (y.ck<=bw.bar[4])

      if (legal)
      {
        screen(1, new=FALSE)

        ## re-draw bandwidth bar to erase bandwidths notch 
        bw.bar <- c(0.2, 0.8, 0.3, 0.5) ## (x1, x2, y1, y2)
        rect(bw.bar[1], bw.bar[3], bw.bar[2], bw.bar[4], col="green", border=FALSE) 

        ## erase bandwidth bar message
        text((bw.bar[1]+bw.bar[2])/2, 0.15, clk.bar.msg1, col=bgCol)

        prev.click.vec <- click.vec
        curr.log10h <- logh.low + sca.fac*(x.ck-bw.bar[1])    
        h <- 10^curr.log10h

        ## draw bandwidths notch
        notch.y <- c(0.90*bw.bar[3]+0.10*bw.bar[4],0.10*bw.bar[3]+0.90*bw.bar[4])
        lines(rep(x.ck[1],2), notch.y,lwd=3,col="darkgreen")

        ## draw bandwidth
        text((bw.bar[1]+bw.bar[2])/2, bw.bar[4]+0.2, bw.text,col=bgCol)
        if (d==1)
          bw.text <- paste("bandwidth =", toString(signif(h,3)))
        else if (d>1)
          bw.text <- paste("bandwidths = (", toString(signif(h,3)), ")", sep="")
        
        text((bw.bar[1]+bw.bar[2])/2, bw.bar[4]+0.2, bw.text,col="purple4")

        if (d==2)
        {  
          ## draw wait symbol
          symbols(x=wait.x, y=wait.y, circle=wait.r, add=TRUE, inches=FALSE,
                  fg="blue", lwd=3)
          text(wait.x, wait.y, "WAIT", col="blue")
        }
        
        dest <- drvkde(gcounts,rep(0,d),bandwidth=h,binned=TRUE,
                       range.x=range.x, se=FALSE)
        dest$est[dest$est<0] <- 0
        prev.dest <- dest

        screen(2, new=TRUE)
        par(mar=old.mar)
        
        if (d==1)
        {
          plot(dest$x.grid[[1]][plot.inds[[1]]],dest$est[plot.inds[[1]]],
               type="n",bty="l",col=densCol,lwd=2, xlim=xlim, ylim=ylim,
               xlab=xlab, ylab="kernel density estimate")
          
          lines(dest$x.grid[[1]][plot.inds[[1]]],dest$est[plot.inds[[1]]],
            bty="l",col=densCol,lwd=2)
        }
        
        if (d==2)
        {
          image(x.grid.1[plot.inds[[1]]],x.grid.2[plot.inds[[2]]],
                dest$est[plot.inds[[1]],plot.inds[[2]]], xlab=xlab,
                ylab=ylab, bty="n", col=densCol, xlim=xlim, ylim=ylim)
          
          lines(c(x1.bor.low,x1.bor.upp),rep(x2.bor.low,2),
                lwd=2,col="black")
          lines(c(x1.bor.low,x1.bor.upp),rep(x2.bor.upp,2),
                lwd=2,col="black")
          lines(rep(x1.bor.low,2),c(x2.bor.low,x2.bor.upp),
                lwd=2,col="black")
          lines(rep(x1.bor.upp,2),c(x2.bor.low,x2.bor.upp),
                lwd=2,col="black")
        }

        feat <- addSignifFeature(h=h, dest=dest)

        screen(1, new=FALSE)
        plot(0, xlim=c(0,1), ylim=c(0,1), type="n", axes=FALSE)
        if (d==2)
        {  
          ## erase wait symbol
          symbols(x=wait.x, y=wait.y, circle=1.05*wait.r, add=TRUE, inches=FALSE,
                  fg=bgCol, bg=bgCol, lwd=3)
        }
        
        ## draw bandwidth bar message
        text((bw.bar[1]+bw.bar[2])/2, 0.15, clk.bar.msg1, col="black")
        #text((bw.bar[1]+bw.bar[2])/2, 0.15, clk.bar.msg2, col="dark green")
        #text((bw.bar[1]+bw.bar[2])/2, 0.15, clk.bar.msg3, col="red")
      }
    }                 
  }

  
  while (!finished & interactive & d >= 3)
  {
    
    rgl.set(screen1)
    rgl.viewpoint(theta=0, phi=0)
    
    sca.fac <- (logh.upp-logh.low)/(bw.bar[2]-bw.bar[1])

    ## RGL analgoue to locator()
    f <- select3d()
    bw.bar.xy <- expand.grid(seq(bw.bar[1], bw.bar[2], length=101),
                             seq(bw.bar[3], bw.bar[4], length=101))
    stop.box.xy <- expand.grid(seq(stop.box[1], stop.box[2], length=101),
                               seq(stop.box[3], stop.box[4], length=101))
   
    click.bw <- f(bw.bar.xy[,1], bw.bar.xy[,2], rep(0, length(bw.bar.xy[,1])))
    click.stop <- f(stop.box.xy[,1], stop.box.xy[,2],rep(0, length(stop.box.xy[,1])))
    
    ## mean of selected region is 'click'
    click.bw.vec <- apply(bw.bar.xy[click.bw,], 2, mean)
    click.stop.vec <- apply(stop.box.xy[click.stop,], 2, mean)
    x.bw.ck <- click.bw.vec[1]
    y.bw.ck <- click.bw.vec[2]
    x.stop.ck <- click.stop.vec[1]
    y.stop.ck <- click.stop.vec[2]
    
    ## Check for finish

    finished <- !(is.na(x.stop.ck) & is.na(y.stop.ck))

    if (!finished)
    {
      ## Check for legality
      legal <- !(is.na(x.bw.ck) & is.na(y.bw.ck))
      
      if (legal)
      {  
        prev.click.bw.vec <- click.bw.vec
        curr.log10h <- logh.low + sca.fac*(x.bw.ck-bw.bar[1])    
        h <- 10^curr.log10h
        
        rgl.set(screen1)

        ## remove old bandwidths, bandwidths bar message, bandwidths notch
        pop3d();  pop3d();  pop3d();  pop3d();  pop3d();
        pop3d();  pop3d();  

        ## draw bandwidths notch
        x.bw.ck <- (log(h, 10) - logh.low)/sca.fac + bw.bar[1]    
        lines3d(c(x.bw.ck[1], x.bw.ck[1]),
                c(0.9*bw.bar[3]+0.1*bw.bar[4], 0.1*bw.bar[3]+0.9*bw.bar[4]),
                c(0.05,0.05), color="darkgreen", size=3)

        ## draw bandwidths 
        texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.2, 0,
                "bandwidths = ", color="purple4", adj=0.3)
        texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.1, 0,  
                paste("(", toString(signif(h,3)), ")", sep=""),
                color="purple4", adj=0.3)
        
        ## draw wait symbol
        spheres3d(wait.x, (bw.bar[3]+bw.bar[4])/2, 0, r=0.05,
                  color="blue", alpha=0.5)
        texts3d(wait.x,bw.bar[3]-0.1,0, "WAIT", color="blue", adj=0.3)    
        
        ## draw KDE & feature sig plot
        
        rgl.set(screen2)
        clear3d()
        rgl.viewpoint(theta=0, phi=-90)
        dest <- drvkde(gcounts,rep(0,d),bandwidth=h,binned=TRUE,
                       range.x=range.x, se=FALSE)
        dest$est[dest$est<0] <- 0
        prev.dest <- dest
         
        if (d==3)
        {        
          x.gd.1 <- dest$x.grid[[1]]; x.gd.2 <- dest$x.grid[[2]]
          x.gd.3 <- dest$x.grid[[3]]

          alph <- seq(0.1,0.5,length=num.levs)
          lev.vals <- seq(0, max(dest$est), length=num.levs+2)[-c(1, num.levs+2)]
          if (addKDE)
            for (il in 1:num.levs)
              contour3d(dest$est, level=lev.vals[il],x=x.gd.1,y=x.gd.2,z=x.gd.3,
                        color=densCol[il],alpha=alph[il],add=(il!=1))
         
        }

        feat <- addSignifFeature(h=h, dest=dest)
        
        rgl.set(screen1)
             
        ## something wrong with stack of RGL objects when swapping between RGL
        ## windows - instead of just pop3d-ing the wait symbol, have to redraw
        ## all objects

        ## erase wait symbol
        ##pop3d(); pop3d(); 

        
        clear3d(type="shapes")
        material3d(alpha=1)
        
        ## draw bandwidth bar
        quads3d(c(bw.bar[1], bw.bar[1], bw.bar[2], bw.bar[2]),
                c(bw.bar[3], bw.bar[4], bw.bar[4], bw.bar[3]),c(0,0,0,0),
                col="green")
        
        ## draw stop box
        quads3d(c(stop.box[1], stop.box[1], stop.box[2], stop.box[2]),
                c(stop.box[3], stop.box[4], stop.box[4], stop.box[3]),
                   c(0,0,0,0), col="red")
        texts3d(c(stop.box[1]+stop.box[2])/2,stop.box[3]-0.1, 0,
                "STOP", col="red", adj=0.3)

        ## draw bandwidths notch
        x.bw.ck <- (log(h, 10) - logh.low)/sca.fac + bw.bar[1]    
        lines3d(c(x.bw.ck[1], x.bw.ck[1]),
                c(0.9*bw.bar[3]+0.1*bw.bar[4], 0.1*bw.bar[3]+0.9*bw.bar[4]),
                c(0.05,0.05), color="darkgreen", size=3)
      
        ## draw bandwidths 
        texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.2, 0,
                "bandwidths = ", color="purple4", adj=0.3)
        texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.1, 0,  
                paste("(", toString(signif(h,3)), ")", sep=""),
                color="purple4", adj=0.3)

        ## draw bandwidth bar message
        texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[3]-0.15, 0,
                "Select rectangular region",col="black",adj=0.3)
        texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[3]-0.25, 0,
                "using mouse button in:",col="black",adj=0.3)
        texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[3]-0.35, 0,
                "green bar to vary bandwidths", col="darkgreen", adj=0.3)
        texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[3]-0.45, 0,
                "red box to stop", col="red", adj=0.3)
      }                
    }
  }

  if (interactive)
  {
    if (d < 3)
    {  
      #text((bw.bar[1]+bw.bar[2])/2, bw.bar[4]+0.4, col=bgCol,
      #     " WARNING: clicks must be within the green bar.\n")
      text((bw.bar[1]+bw.bar[2])/2, bw.bar[4]+0.5 , col="deeppink",
           "\n End of interactive session")
      par(bg="transparent", mar=old.mar)
      close.screen(all=TRUE)
    }
    if (d >=3)
    {
      #rgl.set(screen2)
      #texts3d(xlim[2], ylim[2], zlim[2],
      #        paste("bandwidths = (", toString(signif(h,3)), ")",
      #              sep=""), color="purple4", alpha=1)
      rgl.set(screen1)
      texts3d((bw.bar[1]+bw.bar[2])/2,bw.bar[4]+0.4,0,
              "End of interactive session", color="deeppink", adj=0.3)
    }
  }

  if (plotSiZer)
    feat.temp <- list(x.grid=x.SiZer,bw=exp(bw.SiZer), SiZer=SiZer.map)
  else
  {
    feat.temp <- list(x=x, bw=h, fhat=dest)
    feat.temp <- c(feat.temp, feat)

    class(feat.temp) <- "fs"
  }
  invisible(feat.temp)
}

########## End of featureSignif ##########

