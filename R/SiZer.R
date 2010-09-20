SiZer <- function(x, bw, gridsize, scaleData=FALSE, signifLevel=0.05,  plotSiZer=TRUE, logbw=TRUE,  xlim, xlab, addLegend=TRUE, posLegend="topright") 
{ 
  if (!is.vector(x))
    stop("SiZer is currently only available for 1-dimensional data")
 
  tau <- 5
  
  d <- 1
  n <- length(x)

  if (scaleData)
    for (i in 1:d)
      x[,i] <- (x[,i]-min(x[,i]))/(max(x[,i]) - min(x[,i]))
  
  if (missing(gridsize)) gridsize <- rep(101,2)
  if (length(gridsize)==1) gridsize <- rep(gridsize, 2)
  
  gs <- gridsize[1]
  
  ## Set some defaults

  if (missing(bw))
  {
     bw.range  <- dfltBWrange(x,tau)
     bw <- matrix(unlist(bw.range), nrow=2, byrow=FALSE)
  }
  else
     bw <- matrix(bw, ncol=1, nrow=2)
    
  if (missing(xlim))
  {
    h.low <- bw[1,]
    h.upp <- bw[2,]
    hmix.prop <- 1/4
    h.init <- h.low^(hmix.prop)*h.upp^(1-hmix.prop) ##sqrt(h.low*h.upp)
    
    xlim <- c(min(x)- h.init,max(x)+ h.init)
  }  
    
  if (missing(xlab)) {xlab <- names(x); xlab[is.null(xlab)] <- "" }
    
  dfltCounts.out  <- dfltCounts(x,gridsize=gs, apply(bw, 2, max))
  range.x <-dfltCounts.out$range.x
  gcounts <- dfltCounts.out$counts
  x.SiZer <- seq(range.x[[1]][1], range.x[[1]][2], length=gs) 
  bw <- seq(log(bw[1,1]), log(bw[2,1]), length=gridsize[2])
  SiZer.map <- matrix(0, ncol=length(bw), nrow=length(x.SiZer))
    
  i <- 0
  for (logh in bw) 
  {
    h <- exp(logh)
    i <- i + 1

    est.dens <- drvkde(gcounts,drv=0,bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE)
    est.dens$est[est.dens$est<0] <- 0
    ESS <- n*est.dens$est*prod(h)*(sqrt(2*pi)^d)
    sig.ESS <- ESS >= 5
     
    sig.grad <- SignifFeatureRegion(n,d,gcounts,gridsize=gs, est.dens, h,signifLevel, range.x, grad=TRUE, curv=FALSE)$grad

    est.grad <- drvkde(gcounts, drv=1, bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE)$est
          
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
    
  if (logbw)
     image(x.SiZer, bw, SiZer.map, breaks=c(-1,0,1,2,3), col=c("grey", "purple", "blue", "red"), ylab="log(bandwidth)", xlab=xlab, xlim=xlim)
    else
      image(x.SiZer, exp(bw), SiZer.map, breaks=c(-1,0,1,2,3), col=c("grey", "purple", "blue", "red"), ylab="bandwidth", xlab=xlab, xlim=xlim)

  if(addLegend)
    legend(posLegend, legend=c("sparse data", "zero grad", "+ve grad", "-ve grad"), fill=c("grey", "purple", "blue", "red"), bty="n")
  
  feat <- list(x=x, x.grid=x.SiZer, bw=exp(bw), SiZer=SiZer.map)
  
  invisible(feat)
}
