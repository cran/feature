   crossToCol <- function(cross.vec)
   {
    col.vec <- rep("purple", length = length(cross.vec))	
	
  # vector of colour indicators
  # elements initilialised to "purple"
	
   if(cross.vec[1] == "a") col.vec[1] <- "blue"	# before first crossing
	
   if(cross.vec[1] == "b")
	
	col.vec[1] <- "red"
	
   if(cross.vec[length(cross.vec)] == "a") 

        col.vec[length(cross.vec)] <-"blue" # after last crossing

   if(cross.vec[length(cross.vec)] == "b")

	col.vec[length(cross.vec)] <- "red"

   if(cross.vec[1] == "a") 
   {

    # Case: l_1 >0

     for (i in 2:(length(cross.vec) - 1))
     {
	if(cross.vec[i] == "l") 

        {

	 if (sum(cross.vec[1:i] == "l") %% 2 == 0)
				  
             col.vec[i] <- "blue"
	}


	if(cross.vec[i] == "u") 
        
         {

	 if(sum(cross.vec[1:i] == "u") %% 2 == 1)

            col.vec[i] <- "red"
	 }
      }
    }
	
   if(cross.vec[1] == "b") 
    {

    # Case: u_1<0

    for(i in 2:(length(cross.vec) - 1))  

    {
      if(cross.vec[i] == "l") 
      {
	if(sum(cross.vec[1:i] == "l") %% 2 == 1)
         
          col.vec[i] <- "blue"
       }

     if(cross.vec[i] == "u") 
     {
	if(sum(cross.vec[1:i] == "u") %% 2 == 0)

           col.vec[i] <- "red"
     }
		
     }
   }


  if(cross.vec[1] == "o") 

  {
  # Case: l_1< 0 and u_1 > 0

   for(i in 2:(length(cross.vec) - 1)) 

     {
	if(cross.vec[i] == "l") 
      
         {
	  if(sum(cross.vec[1:i] == "l") %% 2 == 1)

             col.vec[i] <- "blue"
	 }

	if(cross.vec[i] == "u") 
        {
	 if(sum(cross.vec[1:i] == "u") %% 2 == 1)

	   col.vec[i] <- "red"
	}
      }
  }
	
  return(col.vec)
  }
