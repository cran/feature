

fsplot3d <- function(panel)
{
x <- rnorm(panel$n)
y <- rnorm(panel$n)
z <- rnorm(panel$n)
plot3d(x,y,z)
return(panel)
}

ncpanel.width <- 850
ncpanel.height <- 650

ncpanel <- rp.control(panelname="nuclocR.panel", title="nuclocR GUI: Gene maps for cell nuclei", size=c(ncpanel.width, ncpanel.height))
rp.textentry(ncpanel, n, initval="100", pos=c(10, 100, 100, 200))

ncpanel.button.width <- (ncpanel.width - 50)/4

ncpanel.button.start <- (ncpanel.width - 4*ncpanel.button.width - 3*10)/2  ##x.start1 ##+ (ncpanel.width - x.start2 - (x.width3 + x.width4))/2

rp.button(ncpanel, action=fsplot3d, title="plot", pos=c(ncpanel.button.start, ncpanel.height-60, ncpanel.button.width, 50))
 
