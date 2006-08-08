########## R-function: symconv2D ##########

# Computes the discrete two-dimensional convolution of
# a symmetric response matrix rr and data matrix ss.

# Last changed: 20 MAY 2005
 
symconv2D <- function(rr,ss,skewflag=rep(1,2))

{  
   L <- dim(rr)-1
   M <- dim(ss) 
   L1 <- L[1]
   L2 <- L[2]               # find dimensions of r,s
   M1 <- M[1]
   M2 <- M[2]
   P1 <- 2^(ceiling(log(M1+L1)/log(2))) # smallest power of 2 >= M1+L1         
   P2 <- 2^(ceiling(log(M2+L2)/log(2))) # smallest power of 2 >= M2+L2         

   rp <- matrix(0,P1,P2)
   rp[1:(L1+1),1:(L2+1)] <- rr
   if (L1>0)
      rp[(P1-L1+1):P1,1:(L2+1)] <- skewflag[1]*rr[(L1+1):2,]
   if (L2>0)
   {
      rp[1:(L1+1),(P2-L2+1):P2] <- skewflag[2]*rr[,(L2+1):2]
      rp[(P1-L1+1):P1,(P2-L2+1):P2] <-  prod(skewflag)*rr[(L1+1):2,(L2+1):2]   
   }
                                       # wrap around version of rr
   sp <- matrix(0,P1,P2)
   sp[1:M1,1:M2] <- ss                 # zero-padded version of ss
                                       
   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2))[1:M1,1:M2]) # return normalized truncated tt
}

######## End of symconv2D ########
