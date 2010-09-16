# TODO: Add comment
# 
# Author: twutz
###############################################################################

x <- seq(0,1,length.out=41)[-c(1,41)]
dx <- diff(x[1:2])	#differential approximation
lx <- logit(x)			

theta <- twCoefLogitnorm(c(0.7,0.9), perc=c(0.5,0.975))
pcx <- plogitnorm(x,mean=theta[1],sd=theta[2])	#percentiles function
px <- dlogitnorm(x,mean=theta[1],sd=theta[2])	#density function

z <- rlogitnorm(1e6, mean=theta[1],sd=theta[2])
dz <- density(z)
fz <- splinefun(dz$x, dz$y)
(xmax <- optimize(fz,range(dz$x), maximum=TRUE)$maximum)

theta[1]-theta[2]^2*(1-2*xmax) - logit(xmax)

#mtrace(.ofLogitnormMLE)
.ofLogitnormMLE(theta,mle=xmax,quant=0.9)


pmax <- invlogit( +theta[2]^2*(1+theta[1]) )
invlogit( theta[1]-theta[2]^2 )
plot(px~x); abline=c(0.7,)

plot(px~x); abline(v=c(0.7,xmax),col="gray")
lines(dz$y ~ dz$x)

quantile(z,probs=c(0.5,0.975))
mean(z)

sum(px*dx)
plot( cumsum(px*dx) ~ I(x+x[1]/2))	# average across the previous interval
lines( pcx ~ x)
lines( cumsum(dz$y*c(diff(dz$x),1)) ~ dz$x, col="blue")



