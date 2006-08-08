########## R-function: symconv ##########

# Computes the discrete convolution of
# a symmetric or skew-symmetric response 
# vector r and a data vector s.
# If r is symmetric then "skewflag"=1.
# If r is skew-symmetric then "skewflag"=-1.

# Last changed: 03 AUG 2005
 
symconv <- function(r,s,skewflag=1)

{ 
   L <- length(r)-1
   M <- length(s) 
   P <- 2^(ceiling(log(M+L)/log(2))) # smallest power of 2>=M+L         
   r <- c(r,rep(0,P-2*L-1),skewflag*r[(L+1):2])
                                     # wrap-around version of r
   s <- c(s,rep(0,P-M))              # zero-padded version of s
   R <- fft(r)                       # Obtain FFT's of r and s  
   S <- fft(s)
   t <- fft(R*S,TRUE)               # invert element-wise product of FFT's 
   return((Re(t)/P)[1:M])            # return normalized truncated t
}

########## End of symconv ##########

