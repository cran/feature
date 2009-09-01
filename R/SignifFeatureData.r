########## R function: SignifFeatureData ##########

# For determining data points that lie in the
# region of significant
# feature for a particular bandwidth and
# significance level.

# Last changed: 08 DEC 2005

SignifFeatureData <- function(x, d, dest, SignifFeature)
{
  n <- nrow(x)
  x.ind <- matrix(0, ncol=d, nrow=n)
  signif <- rep(FALSE, length=n)

  for (i in 1:n)
  {  
    for (j in 1:d)
      x.ind[i,j] <- sum(x[i,j] >= dest$x.grid[[j]])

    ##if (!is.null(SignifFeature))
    if (d==1)
      sigf <- SignifFeature[x.ind[i,1]]
    else if (d==2)
      sigf <- SignifFeature[x.ind[i,1], x.ind[i,2]]
    else if (d==3)
      sigf <- SignifFeature[x.ind[i,1], x.ind[i,2], x.ind[i,3]]
    else if (d==4)
      sigf <- SignifFeature[x.ind[i,1], x.ind[i,2], x.ind[i,3], x.ind[i,4]]
    if (!is.null(sigf)) signif[i] <- sigf  
  }  
  
  return(signif)
}

########## End of SignifFeatureData ##########
