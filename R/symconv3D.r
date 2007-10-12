########### R-function: symconv3D ##########

# Computes the discrete three-dimensional convolution of a symmetric  
# response array rr and a data matrix ss.

# Last changed: 01 JUN 2005
 
symconv3D <- function(rr,ss,skewflag=rep(1,3))

{  
   L <- dim(rr) - 1

   M <- dim(ss) 
   P <- 2^(ceiling(log(M+L)/log(2))) # smallest powers of 2 >= M+L
   L1 <- L[1] ; L2 <- L[2] ; L3 <- L[3]
   M1 <- M[1] ; M2 <- M[2] ; M3 <- M[3]               
   P1 <- P[1] ; P2 <- P[2] ; P3 <- P[3]
   sf <- skewflag

   rp <- array(0,P) 
   rp[1:(L1+1),1:(L2+1),1:(L3+1)] <- rr
   rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1)] <- sf[1]*rr[(L1+1):2,1:(L2+1),1:(L3+1)]
   rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1)] <- sf[2]*rr[1:(L1+1),(L2+1):2,1:(L3+1)]
   rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3] <- sf[3]*rr[1:(L1+1),1:(L2+1),(L3+1):2]
   rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1)] <-  
                                   sf[1]*sf[2]*rr[(L1+1):2,(L2+1):2,1:(L3+1)]
   rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3] <-  
                                   sf[2]*sf[3]*rr[1:(L1+1),(L2+1):2,(L3+1):2]
   rp[(P1-L1+1):P1,1:(L2+1),(P3-L3+1):P3] <-  
                                   sf[1]*sf[3]*rr[(L1+1):2,1:(L2+1),(L3+1):2]
   rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3] <-  
                             sf[1]*sf[2]*sf[3]*rr[(L1+1):2,(L2+1):2,(L3+1):2]

   sp <- array(0,P)
   sp[1:M1,1:M2,1:M3] <- ss            # zero-padded version of ss

   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)   
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2*P3))[1:M1,1:M2,1:M3]) # return normalized truncated tt
}

########## End of symconv3D ###########
