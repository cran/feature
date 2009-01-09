barToCross <- function(lower,upper)
{
   ng <- length(lower)

    prod.l <- lower[1:(ng - 1)] * lower[2:ng]

    z.c.l <- numeric(0)

    if(any(prod.l < 0))

      z.c.l <- (1:ng)[prod.l < 0]

    prod.u <- upper[1:(ng - 1)] * upper[2:ng]

    z.c.u <- numeric(0)

    if(any(prod.u < 0)) 
       
      z.c.u <- (1:ng)[prod.u < 0]

   # Determine crossing sequence

    cross.vec <- "o"

    if(lower[1] > 0)

	cross.vec <- "a"

    if(upper[1] < 0)

	cross.vec <- "b"

    z.c <- cbind(z.c.l, rep(0, length(z.c.l)))	

    # Then crossings of lower CI

    z.c <- rbind(z.c, cbind(z.c.u, rep(1, length(z.c.u))))

    if(length(z.c) > 0) 
       {
	#z.c <- sortcol(z.c, 1) 
	z.c <- z.c[order(z.c[,1]),]
	
	if(length(z.c) == 2) z.c <- t(as.matrix(z.c))
	
	# Check for ties

	num.ties <- length(z.c[, 1]) - length(unique(z.c[, 1]))

	if(num.ties > 0) 
        {
         # Ties exist in crossing that needs to be
         # sort out

	diffs <- z.c[2:(nrow(z.c)), 1] - z.c[1:(nrow(z.c) - 1), 1]

	tie.pos <- (1:length(diffs))[diffs == 0]

	for(itp in 1:length(tie.pos)) 
        {
	 if(tie.pos[itp] == 1) 
          {
	    if(lower[z.c[1, 1]] > 0) 
             {
	      z.c[tie.pos[itp], 2] <- 0
              z.c[tie.pos[itp] + 1, 2] <- 1
             }

        if(upper[z.c[1, 1]] < 0) 
         {
          z.c[tie.pos[itp], 2] <- 1
          z.c[tie.pos[itp] + 1, 2] <- 0
	 }
	}

	if(tie.pos[itp] > 1) 
         {
	  if(z.c[tie.pos[itp] - 1, 2] == 0) 
            {
            z.c[tie.pos[itp], 2] <- 0
	    z.c[tie.pos[itp] + 1, 2] <- 1
	    }
	  if(z.c[tie.pos[itp] - 1, 2] == 1) 
            {
	     z.c[tie.pos[itp], 2] <- 1
	     z.c[tie.pos[itp] + 1, 2] <- 0
	    }
	  }
	}
       }

   c.v.add <- z.c[, 2]
   c.v.add[c.v.add == 0] <- "l"
   c.v.add[c.v.add == 1] <- "u"

   for(icv in 1:length(c.v.add))

     cross.vec <- cbind(cross.vec, c.v.add[icv])

   }
	
   if(lower[ng] > 0)

	cross.vec <- cbind(cross.vec,"a") # Significance at last gridpoint

   if(upper[ng] < 0)

	cross.vec <- cbind(cross.vec, "b")


   if((lower[ng] <= 0) & (upper[ng] >= 0))

	cross.vec <- cbind(cross.vec, "o")


    bdy.inds <- c(1, sort(c(z.c.l, z.c.u)), ng)
 
    return(list(bdy.inds=bdy.inds,cross.vec=cross.vec))

  }
