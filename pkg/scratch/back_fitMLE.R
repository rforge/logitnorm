# TODO: Add comment
# 
# Author: twutz
###############################################################################


.ofLogitnormMLE <- function(
	### Objective function used by \code{\link{coefLogitnormMLE}}. 
	theta				##<< theta[1] is the mean, theta[2] is the sd
	,mle				##<< the mode of the density distribution
	,quant				##<< q are the quantiles for perc
	,perc=c(0.975) 
){
	# there is an analytical expression for the gradient, but here just use numerical gradient
	
	# calculating percentiles should be faster than calculating quantiles of the normal distr.
	tmp.predp = plogitnorm(quant, mean=theta[1], sd=theta[2] )
	tmp.diff = tmp.predp - perc
	tmp.diff.mle <- theta[1]-theta[2]^2*(1-2*mle)   - logit(mle)
	sum(tmp.diff^2) + tmp.diff.mle^2
}

twCoefLogitnormMLE <- function(
	### Estimating coefficients of logitnormal distribution from mode and upper quantile	
	mle						##<< numeric vector: the mode of the density function
	,quant					##<< numeric vector: the upper quantile value
	,perc=0.975				##<< numeric vector: the probability for which the quantile was specified
	#,method="BFGS"			##<< method of optimization (see \code{\link{optim}})
	,theta0=c(mean=0,sd=1)	##<< starting parameters
	,returnDetails=FALSE	##<< if TRUE, the full output of optim is attached as attributes resOptim
	, ... 					##<< Further parameters to \code{\link{optim}}, like method or \code{control=list(maxit=1000)}
){
	##seealso<< \code{\link{logitnorm}}
	# twCoefLogitnormMLE
	names(theta0) <- c("mean","sd")
	nc <- c(length(mle),length(quant),length(perc)) 
	n <- max(nc)
	res <- matrix( as.numeric(NA), n,2, dimnames=list(NULL,c("mean","sd")))
	resOptim <- list()
	for( i in 1:n ){
		i0 <- i-1
		tmp <- optim( theta0, .ofLogitnormMLE, mle=(mleI<-mle[1+i0%%nc[1]]), quant=(quantI<-quant[1+i0%%nc[2]]), perc=(percI<-perc[1+i0%%nc[3]]), ...)
		if( tmp$convergence == 0)
			res[i,] <- tmp$par
		else
			warning(paste("coefLogitnorm: optim did not converge. theta0=",paste(theta0,collapse=","),"; mle=",mleI,"; quant=",quantI,"; perc=",percI,sep=""))
		resOptim[[i]] <- tmp
	}
	if( returnDetails ) attr(res,"resOptim") <- resOptim
	res
	### numeric matrix with columns \code{c("mean","sd")}
	### rows correspond to rows in mle, quant, and perc
}
attr(twCoefLogitnormMLE,"ex") <- function(){
	x <- seq(0,1,length.out=41)[-c(1,41)]	# plotting grid
	
	# estimate the parameters, with mode 0.7 and upper quantile 0.9
	(theta <- twCoefLogitnormMLE(0.7,0.9))
	
	px <- plogitnorm(x,mean=theta[1],sd=theta[2])	#percentiles function
	plot(px~x); abline(v=c(0.7,0.9),col="gray"); abline(h=c(0.975),col="gray")
	
	dx <- dlogitnorm(x,mean=theta[1],sd=theta[2])	#density function
	plot(dx~x); abline(v=c(0.7,0.9),col="gray")
	
	# vectorized, including an impossible combination of mean and quantile
	(theta <- twCoefLogitnormMLE(seq(0.4,0.8,by=0.1),0.9))
}


