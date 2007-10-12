######### R-function dfltBWrange  #########
 
# Obtain default set of grid counts from a 
# multivariate point cloud 'x'.

# Last changed: 22 JUL 2005

dfltBWrange <- function(x,gridsize,tau,scale.fac=1.2)
{
   d <- NCOL(x)
   if (d==1) x <- as.matrix(x)

   dmn.fac <- (4/(d+2))^(1/(d+4))
   samp.size.fac <- nrow(x)^(-1/(d+4))
   cmb.fac <- dmn.fac*samp.size.fac

   # Compute the scale in each direction

   st.devs <- sqrt(apply(x,2,var))
   Q1.vals <- apply(x,2,quantile,1/4)
   Q3.vals <- apply(x,2,quantile,3/4)
   corr.fac <- qnorm(3/4) - qnorm(1/4)
   IQR.vals <- (Q3.vals - Q1.vals)/corr.fac
   sig.hats <- apply(cbind(st.devs,IQR.vals),1,min)
   range.vals <- apply(x,2,max) - apply(x,2,min)

   range.h <- list()
   for (id in 1:d)
   {
      h.upp <- cmb.fac*scale.fac*sig.hats[id]
      h.low <- 3*(range.vals[id] + 2*tau*h.upp)/((gridsize[id]-1)*tau)
      range.h[[id]] <- c(h.low,h.upp)
  }
   return(range.h)
}

######## End of dfltBWrange ########
