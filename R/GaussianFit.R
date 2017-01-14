
GaussFunBiSig <- function(amp,pos,sig1,sig2,vector)
{
	z.par1 <- (vector-pos)/sig1
	z.par2 <- (vector-pos)/sig2
	
	gaussian.function <- rep(0,length(vector))
	gaussian.function[1:pos] <- amp*(exp(-(z.par1^2)/2))[1:pos]
	gaussian.function[(pos+1):length(vector)] <- amp*(exp(-(z.par2^2)/2))[(pos+1):length(vector)]
	gaussian.function
}

fitGaussian <- function(x)
{
	
	#x <- mzQ.chrom[,1]
	
	x.fitted <- x*0
	Ap <- max(x)
	pos <- which.max(x)
	vect <- 1:length(x)
	nls.fit <- try(nls(x~GaussFunBiSig(Ap,pos,sig1,sig2,vect), start=list(sig1=1,sig2=1)),silent=T)
	if(class(nls.fit)!="try-error")
	{
		sigma1 <- summary(nls.fit)$coefficients[1,1]
		sigma2 <- summary(nls.fit)$coefficients[2,1]
		mask <- GaussFunBiSig(1,pos,sigma1,sigma2,vect)
		mask[normalize(mask)<0.05] <- 0
		mask[mask!=0] <- 1
		x.fitted <-  x*mask	
	}else{
		x.fitted <- x
	} 
	x.fitted	
}
