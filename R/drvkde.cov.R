
########## R-function: drvkde.cov ##########

# Computes the mth derivative and covariances
# of a binned
# d-variate kernel density estimate based
# on grid counts.

# Last changed: 28 OCT 2005

drvkde.cov <- function(x,drv,bandwidth,gridsize,range.x,
                   binned=FALSE,drv2,se=TRUE)
{  
   d <- length(drv)

   if (d==1) x <- as.matrix(x)

   if (binned)
   {
      if (any(dim(x)==1))
         gridsize <- max(dim(x))  

      if (!any(dim(x)==1))
         gridsize <- dim(x)
      gcounts <- x
   }

   if ((!binned)&(missing(gridsize))) gridsize <- rep(64,d)

   # Rename common variables

   h <- bandwidth
   tau <- 4 + max(drv)    

   if (length(h)==1) h <- rep(h,d)

   if (missing(range.x)&(binned==FALSE)) 
   {
      range.x <- list()
      for (id in 1:d)
         range.x[[id]] <- c(min(x[,id])-1.5*h[id],max(x[,id])+1.5*h[id])  
   }

   a <- unlist(lapply(range.x,min))
   b <- unlist(lapply(range.x,max))

   M <- gridsize
   gpoints <- list()

   for (id in 1:d)
      gpoints[[id]] <- seq(a[id],b[id],length=M[id])

   # Bin the data if not already binned

   if (binned==FALSE)
   {
      if (d==1) 
        gcounts <- linbin(x,gpoints[[1]])
      if (d==2) 
        gcounts <- linbin2D(x,gpoints[[1]],gpoints[[2]])
      if (d==3) 
        gcounts <- linbin3D(x,gpoints[[1]],gpoints[[2]],gpoints[[3]])
      if (d==4)
        gcounts <- linbin4D(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]])
   }
   else
      gcounts <- x      

   n <- sum(gcounts)

   if (missing(drv2))  # need to specify drv2 for covariance estimation
   {  
     kapmid <- list()
     for (id in (1:d))
     {
       Lid <- min(floor(tau*h[id]*(M[id]-1)/(b[id]-a[id])),M[id]) 
       lvecid <- (0:Lid)
       facid  <- (b[id]-a[id])/(h[id]*(M[id]-1))
       argid <- lvecid*facid
       kapmid[[id]] <- dnorm(argid)/(h[id]^(drv[id]+1))
       hmold0 <- 1
       hmold1 <- argid
       if (drv[id]==0) hmnew <- 1
       if (drv[id]==1) hmnew <- argid
       if (drv[id] >= 2) 
         for (ihm in (2:drv[id])) 
         {
           hmnew <- argid*hmold1 - (ihm-1)*hmold0
           hmold0 <- hmold1   # Compute drv[id] degree Hermite polynomial
           hmold1 <- hmnew    # by recurrence.
         }
       kapmid[[id]] <- hmnew*kapmid[[id]]*(-1)^drv[id]
     }
   }
   else
   {
     kapmid1 <- list()
     kapmid2 <- list()
     for (id in (1:d))
     {
       Lid <- min(floor(tau*h[id]*(M[id]-1)/(b[id]-a[id])),M[id]) 
       lvecid <- (0:Lid)
       facid  <- (b[id]-a[id])/(h[id]*(M[id]-1))
       argid <- lvecid*facid
       kapmid1[[id]] <- dnorm(argid)/(h[id]^(drv[id]+1))
       kapmid2[[id]] <- dnorm(argid)/(h[id]^(drv2[id]+1))

       hmold0 <- 1
       hmold1 <- argid
       if (drv[id]==0) hmnew <- 1
       if (drv[id]==1) hmnew <- argid
       if (drv[id] >= 2) 
         for (ihm in (2:drv[id])) 
         {
           hmnew <- argid*hmold1 - (ihm-1)*hmold0
           hmold0 <- hmold1   # Compute drv[id] degree Hermite polynomial
           hmold1 <- hmnew    # by recurrence.
         }
       kapmid1[[id]] <- hmnew*kapmid1[[id]]*(-1)^drv[id]

       hmold0 <- 1
       hmold1 <- argid
       if (drv2[id]==0) hmnew <- 1
       if (drv2[id]==1) hmnew <- argid
       if (drv2[id] >= 2) 
         for (ihm in (2:drv2[id])) 
         {
           hmnew <- argid*hmold1 - (ihm-1)*hmold0
           hmold0 <- hmold1   # Compute drv[id] degree Hermite polynomial
           hmold1 <- hmnew    # by recurrence.
         }
       kapmid2[[id]] <- hmnew*kapmid2[[id]]*(-1)^drv2[id] 
     }     
   }

   if (d==1)
      kappam <- kapmid[[1]]/n

   if (d==2)
   {
     if (missing(drv2)) # need to specify drv2 for covariance estimation
       kappam <- outer(kapmid[[1]],kapmid[[2]])/n
     else
     {
       kappam1 <- outer(kapmid1[[1]],kapmid1[[2]])/n
       kappam2 <- outer(kapmid2[[1]],kapmid2[[2]])/n
     }
   }
   
   if (d==3)
   {
     if (missing(drv2)) # need to specify drv2 for covariance estimation
       kappam <- outer(kapmid[[1]],outer(kapmid[[2]],kapmid[[3]]))/n
     else
     {
       kappam1 <- outer(kapmid1[[1]],outer(kapmid1[[2]],kapmid1[[3]]))/n
       kappam2 <- outer(kapmid2[[1]],outer(kapmid2[[2]],kapmid2[[3]]))/n
     }
   } 

   if (d==4)
   {
     if(missing(drv2)) # need to specify drv2 for covariance estimation
       kappam <- outer(kapmid[[1]],
                       outer(kapmid[[2]],outer(kapmid[[3]],kapmid[[4]])))/n
     else
     {
       kappam1 <- outer(kapmid1[[1]],
                        outer(kapmid1[[2]],outer(kapmid1[[3]],kapmid1[[4]])))/n
       kappam2 <- outer(kapmid2[[1]],
                        outer(kapmid2[[2]],outer(kapmid2[[3]],kapmid2[[4]])))/n
     }
   }

   if (!any(c(d==1,d==2,d==3,d==4))) stop("only for d=1,2,3,4")

   if (d==1) 
   { 
      kappam <- as.vector(kappam)
      est <- symconv(kappam,gcounts,skewflag=(-1)^drv)
      if (!se)
        return(list(x.grid=gpoints,est=est))
      
      est.var <- ((symconv((n*kappam)^2,gcounts)/n) - est^2)/(n-1)
      est.var[est.var<0] <- 0
      return(list(x.grid=gpoints,est=est,se=sqrt(est.var)))
   }

   if (d==2) 
   {     
     if (missing(drv2))
     {  
       est <- symconv2D(kappam,gcounts,skewflag=(-1)^drv)
       if (!se)
        return(list(x.grid=gpoints,est=est))
       
       est.var <- ((symconv2D((n*kappam)^2,gcounts)/n) - est^2)/(n-1)
       est.var[est.var<0] <- 0
       return(list(x.grid=gpoints,est=est,se=sqrt(est.var)))
     }
     else
     {
       est1 <- symconv2D(kappam1,gcounts,skewflag=(-1)^drv)
       if (!se)
        return(list(x.grid=gpoints,est=est1))
       
       est.var1 <- ((symconv2D((n*kappam1)^2,gcounts)/n) - est1^2)/(n-1)
       est.var1[est.var1<0] <- 0
       est2 <- symconv2D(kappam2,gcounts,skewflag=(-1)^drv2)
      
       est.cov <- ((symconv2D(n*kappam1*n*kappam2, gcounts,
                              skewflag=(-1)^(drv+drv2))/n) - est1*est2)/(n-1)

       return(list(x.grid=gpoints,est=est1,se=sqrt(est.var1), cov=est.cov))
     }
       
   }

   if (d==3)
   {
     if (missing(drv2))
     {
       est <- symconv3D(kappam,gcounts,skewflag=(-1)^drv) 
       if (!se)
        return(list(x.grid=gpoints,est=est))

       est.var <- ((symconv3D((n*kappam)^2,gcounts)/n) - est^2)/(n-1)
       est.var[est.var<0] <- 0
       return(list(x.grid=gpoints,est=est,se=sqrt(est.var))) 
     }
     else
     {
       est1 <- symconv3D(kappam1,gcounts,skewflag=(-1)^drv)
       if (!se)
        return(list(x.grid=gpoints,est=est1))
       
       est.var1 <- ((symconv3D((n*kappam1)^2,gcounts)/n) - est1^2)/(n-1)
       est.var1[est.var1<0] <- 0
       est2 <- symconv3D(kappam2,gcounts,skewflag=(-1)^drv2)
      
       est.cov <- ((symconv3D(n*kappam1*n*kappam2, gcounts,
                              skewflag=(-1)^(drv+drv2))/n) - est1*est2)/(n-1)

       return(list(x.grid=gpoints,est=est1,se=sqrt(est.var1), cov=est.cov))
     }
     
   }
   if (d==4)
   {
     if (missing(drv2))
     {
       est <- symconv4D(kappam,gcounts,skewflag=(-1)^drv) 
       if (!se)
         return(list(x.grid=gpoints,est=est))

       est.var <- ((symconv4D((n*kappam)^2,gcounts)/n) - est^2)/(n-1)
       est.var[est.var<0] <- 0
       return(list(x.grid=gpoints,est=est,se=sqrt(est.var))) 
     }
     else
     {
       est1 <- symconv4D(kappam1,gcounts,skewflag=(-1)^drv)

       if (!se)
        return(list(x.grid=gpoints,est=est1))
       
       est.var1 <- ((symconv4D((n*kappam1)^2,gcounts)/n) - est1^2)/(n-1)
       est.var1[est.var1<0] <- 0
       est2 <- symconv4D(kappam2,gcounts,skewflag=(-1)^drv2)
      
       est.cov <- ((symconv4D(n*kappam1*n*kappam2, gcounts,
                              skewflag=(-1)^(drv+drv2))/n) - est1*est2)/(n-1)

       return(list(x.grid=gpoints,est=est1,se=sqrt(est.var1), cov=est.cov))
     }
   }
}

########## End of drvkde.cov #########
