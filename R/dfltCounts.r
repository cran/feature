######### R-function:dfltCounts  #########
 
# Obtain default set of grid counts from a 
# multivariate point cloud 'x'.

# Last changed: 18 JUL 2005

dfltCounts <- function(x,gridsize=rep(64,NCOL(x)),h=rep(0,NCOL(x)))
{
   x <- as.matrix(x)
   d <- ncol(x)

   range.x <- list()
   for (id in 1:d)
      range.x[[id]] <- c(min(x[,id])-1.5*h[id],max(x[,id])+1.5*h[id])  

   a <- unlist(lapply(range.x,min))
   b <- unlist(lapply(range.x,max))

   gpoints <- list()
   for (id in 1:d)
      gpoints[[id]] <- seq(a[id],b[id],length=gridsize[id])  
 
   if ((d!=1)&(d!=2)&(d!=3)&(d!=4)) stop("currently only for d=1,2,3,4")

   if (d==1)
      gcounts <- linbin(x,gpoints[[1]])

   if (d==2)
      gcounts <- linbin2D(x,gpoints[[1]],gpoints[[2]])

   if (d==3)
      gcounts <- linbin3D(x,gpoints[[1]],gpoints[[2]],gpoints[[3]])
   
   if (d==4)
      gcounts <- linbin4D(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]])
   
   return(list(counts=gcounts,range.x=range.x))
}

######## End of dfltCounts ########
