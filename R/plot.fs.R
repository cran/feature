plot.fs <- function(x, ...)
{
  plotfs(x, ...)
  invisible() 
}

plotfs <-  function(fs, xlab, ylab, zlab, xlim, ylim, zlim,
           addData=FALSE, scaleData=FALSE, addDataNum=1000,
           addKDE=TRUE, jitterRug=TRUE,  
           addSignifGradRegion=FALSE, addSignifGradData=FALSE,
           addSignifCurvRegion=FALSE, addSignifCurvData=FALSE,
           addAxes3d=TRUE, 
           densCol, dataCol="black", gradCol="green4", curvCol="blue",
           axisCol="black", bgCol="white",
           dataAlpha=0.1, gradDataAlpha=0.3,
           gradRegionAlpha=0.2, curvDataAlpha=0.3, curvRegionAlpha=0.3, gridsize) 
{
  names.x <- colnames(fs$x)
  
  x <- as.matrix(fs$x)
  d <- ncol(x)
  n <- nrow(x)
  h <- fs$bw
  

  ## Determine default axis labels.

  if (missing(xlab)) xlab <- NULL
  if (missing(ylab)) ylab <- NULL
  if (missing(zlab)) zlab <- NULL
  labs <- dfltLabs(d,names.x,xlab,ylab,zlab)
  xlab <- labs$xlab ; ylab <- labs$ylab ; zlab <- labs$zlab
  
  if (missing(gridsize)) 
  {
    if (d==1) gridsize <- 401
    if (d==2) gridsize <- rep(151,2)
    if (d==3) gridsize <- rep(51,3)
    if (d==4) gridsize <- rep(21,4)
  }

 
  dfltCounts.out <- dfltCounts(x, gridsize, fs$bw)
  gcounts <- dfltCounts.out$counts
  range.x <- dfltCounts.out$range.x  

  dest <- drvkde(gcounts, rep(0,d), bandwidth=h, binned=TRUE,
                 range.x=range.x, se=FALSE)
  
  
  ## random sample of data points used for display
  nsamp <- min(addDataNum, n)
  
  if (nsamp < n)
  {
    rand.inds <- sort(sample(1:n, nsamp, replace=FALSE))
    x.rand <- as.matrix(x[rand.inds,])
  }
  else
    x.rand <- x
  
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

  
  plot.inds <- list()
  for (id in 1:d)
  {
    plot.inds.l <- (1:gridsize[id])[dest$x.grid[[id]]>=lims[[id]][1]]
    plot.inds.u <- (1:gridsize[id])[dest$x.grid[[id]]<=lims[[id]][2]]
    plot.inds[[id]] <- intersect(plot.inds.l,plot.inds.u)
  }
  
  if (missing(densCol))
    if (d==1)
      densCol <- "DarkOrange"
    else if (d==2)
      densCol <- rev(heat.colors(1000))
    else if (d==3)
      densCol <- rev(heat.colors(5))

  if (d==1)
  {
    par(bg=bgCol)
    plot(dest$x.grid[[1]][plot.inds[[1]]], dest$est[plot.inds[[1]]],
         type="n",bty="l" ,col=densCol, lwd=2, xlim=xlim, ylim=ylim,
         xlab=xlab,ylab="kernel density estimate")
    
    lines(dest$x.grid[[1]][plot.inds[[1]]],dest$est[plot.inds[[1]]],
          bty="l",col=densCol,lwd=2)
 
    SignifGradRegion.mat <- fs$grad
    SignifCurvRegion.mat <- fs$curv
    
    if (addSignifGradData)  
      SignifGradData.mat <- SignifFeatureData(x.rand, d, dest,SignifGradRegion.mat)
    
    if (addSignifCurvData)
      SignifCurvData.mat <- SignifFeatureData(x.rand, d, dest,SignifCurvRegion.mat)
    
    if (addSignifGradRegion)
      addSignifFeatureRegion(d,gridsize,SignifGradRegion.mat,plot.inds,gradCol,
                             dest,lims)
      
    if (addSignifCurvRegion)
      addSignifFeatureRegion(d,gridsize,SignifCurvRegion.mat,plot.inds,curvCol,
                             dest,lims)
    
    if (addSignifGradData)
      addSignifFeatureData(x.rand,SignifGradData.mat,gradCol)
    
    if (addSignifCurvData)
      addSignifFeatureData(x.rand,SignifCurvData.mat,curvCol)
    
    if (addData)
    {
      if (jitterRug)
        x.rug <- jitter(x.rand)
      else
        x.rug <- x.rand
      rug(x.rug)
    }  
  }
  else if (d==2)
  {
     par(bg=bgCol)
     x.grid.1 <- dest$x.grid[[1]] ; x.grid.2 <- dest$x.grid[[2]]
    
    if (addKDE)
    {  
      image(x.grid.1[plot.inds[[1]]],x.grid.2[plot.inds[[2]]],
             dest$est[plot.inds[[1]],plot.inds[[2]]],col=densCol,
            xlim=xlim, ylim=ylim, xlab=xlab,ylab=ylab,bty="n")
    }
    else
    {
      image(x.grid.1[plot.inds[[1]]],x.grid.2[plot.inds[[2]]],
             dest$est[plot.inds[[1]],plot.inds[[2]]],col="transparent",
            xlim=xlim, ylim=ylim, xlab=xlab,ylab=ylab,bty="n")
    }
    ## Add a border around the image
    
    x1.bor.low <- min(x.grid.1[plot.inds[[1]]]) 
    x1.bor.upp <- max(x.grid.1[plot.inds[[1]]]) 
    x2.bor.low <- min(x.grid.2[plot.inds[[2]]]) 
    x2.bor.upp <- max(x.grid.2[plot.inds[[2]]]) 

    lines(c(x1.bor.low,x1.bor.upp),rep(x2.bor.low,2),lwd=2,col="black")
    lines(c(x1.bor.low,x1.bor.upp),rep(x2.bor.upp,2),lwd=2,col="black")
    lines(rep(x1.bor.low,2),c(x2.bor.low,x2.bor.upp),lwd=2,col="black")
    lines(rep(x1.bor.upp,2),c(x2.bor.low,x2.bor.upp),lwd=2,col="black")
    
    SignifGradRegion.mat <- fs$grad
    SignifCurvRegion.mat <- fs$curv
    
    if (addSignifGradData)  
      SignifGradData.mat <- SignifFeatureData(x.rand, d, dest,SignifGradRegion.mat)
    
    if (addSignifCurvData)
      SignifCurvData.mat <- SignifFeatureData(x.rand, d, dest,SignifCurvRegion.mat)
    
    if (addSignifGradRegion)
      addSignifFeatureRegion(d,gridsize,SignifGradRegion.mat,plot.inds,gradCol,
                             dest,lims)
      
    if (addSignifCurvRegion)
      addSignifFeatureRegion(d,gridsize,SignifCurvRegion.mat,plot.inds,curvCol,
                             dest,lims)
    
    if (addSignifGradData)
      addSignifFeatureData(x.rand,SignifGradData.mat,gradCol)
    
    if (addSignifCurvData)
      addSignifFeatureData(x.rand,SignifCurvData.mat,curvCol)
    
    if (addData)
      points(x.rand, col=dataCol)  
  }
  else if (d==3)
  {
    clear3d()
    rgl.viewpoint(theta=0, phi=-90)
    rgl.bg(col=bgCol)
    pop3d(type="lights")
    light3d(theta=0, phi=30)

    material3d(alpha=1)
    material3d(back="fill")

    num.levs <- length(densCol)
    x.gd.1 <- dest$x.grid[[1]] ; x.gd.2 <- dest$x.grid[[2]]
    x.gd.3 <- dest$x.grid[[3]]
    
    if (addKDE)
    { 
      alph <- seq(0.1,0.5,length=num.levs)
      lev.vals <- seq(0, max(dest$est), length=num.levs+2)[-c(1, num.levs+2)]

      for (il in 1:num.levs)
        contour3d(dest$est,level=lev.vals[il],
                  x=x.gd.1,y=x.gd.2,z=x.gd.3,color=densCol[il],alpha=alph[il],
                  add=(il!=1))   
    }
    SignifGradRegion.mat <- fs$grad
    SignifCurvRegion.mat <- fs$curv
    if (addSignifGradData)
      SignifGradData.mat <- SignifFeatureData(x.rand, d, dest,SignifGradRegion.mat)
       
    if (addSignifCurvData)
      SignifCurvData.mat <- SignifFeatureData(x.rand, d, dest,SignifCurvRegion.mat)
    
    if (addSignifGradRegion)
      addSignifFeatureRegion(d,gridsize,SignifGradRegion.mat,plot.inds,gradCol,
                             dest,lims,trans.alpha=gradRegionAlpha)
     
    if (addSignifCurvRegion)
      addSignifFeatureRegion(d,gridsize,SignifCurvRegion.mat,plot.inds,curvCol,
                             dest,lims,trans.alpha=curvRegionAlpha)
   
    if (addSignifGradData)
      addSignifFeatureData(x.rand,SignifGradData.mat,gradCol, trans.alpha=gradDataAlpha)
    if (addSignifCurvData)
      addSignifFeatureData(x.rand,SignifCurvData.mat,curvCol, trans.alpha=curvDataAlpha)

    if (addData)
      points3d(x.rand[,1],x.rand[,2],x.rand[,3],size=3,col=dataCol, alpha=dataAlpha)

    if (addAxes3d)
    {
      lines3d(xlim, rep(ylim[1],2), rep(zlim[1],2), size=3, color=axisCol, alpha=1)
      lines3d(rep(xlim[1],2), ylim, rep(zlim[1],2), size=3, color=axisCol, alpha=1)
      lines3d(rep(xlim[1],2), rep(ylim[1],2), zlim, size=3, color=axisCol, alpha=1)
      
      texts3d(xlim[2],ylim[1],zlim[1],xlab,size=3, color=axisCol, adj=0, alpha=1)
      texts3d(xlim[1],ylim[2],zlim[1],ylab,size=3, color=axisCol, adj=1, alpha=1)
      texts3d(xlim[1],ylim[1],zlim[2],zlab,size=3, color=axisCol, adj=1, alpha=1)
    }
    
  }
}
  
