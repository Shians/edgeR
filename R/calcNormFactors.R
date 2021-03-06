calcNormFactors <- function(object, ...)
UseMethod("calcNormFactors")

calcNormFactors.DGEList <- function(object, method=c("TMM","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Scale normalization of RNA-Seq data, for DGEList objects
#	Created 2 October 2014.  Last modified 27 August 2015.
{
	object$samples$norm.factors <- calcNormFactors(object=object$counts, lib.size=object$samples$lib.size, method=method, refColumn=refColumn, logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff, p=p)
	object
}

calcNormFactors.default <- function(object, lib.size=NULL, method=c("TMM","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Scale normalization of RNA-Seq data, for count matrices
#	Mark Robinson.  Edits by Gordon Smyth.
#	Created October 22 October 2009 by Mark Robinson.
#	Last modified 31 July 2015.
{
#	Check object
	x <- as.matrix(object)
	if(any(is.na(x))) stop("NA counts not permitted")

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(x)
	if(any(is.na(lib.size))) stop("NA lib.sizes not permitted")

#	Check method
	method <- match.arg(method)

#	Remove all zero rows
	allzero <- .rowSums(x>0, nrow(x), ncol(x)) == 0
	if(any(allzero)) x <- x[!allzero,,drop=FALSE]

#	Degenerate cases
	if(nrow(x)==0 || ncol(x)==1) method="none"

#	Calculate factors
	f <- switch(method,
		TMM = {
			f75 <- .calcFactorQuantile(data=x, lib.size=lib.size, p=0.75)
			if( is.null(refColumn) ) refColumn <- which.min(abs(f75-mean(f75)))
			if(length(refColumn)==0 | refColumn < 1 | refColumn > ncol(x)) refColumn <- 1
			f <- rep(NA,ncol(x))
			for(i in 1:ncol(x))
				f[i] <- .calcFactorWeighted(obs=x[,i],ref=x[,refColumn], libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
			f
		},
		RLE = .calcFactorRLE(x)/lib.size,
		upperquartile = .calcFactorQuantile(x,lib.size,p=p),
		none = rep(1,ncol(x))
	)

#	Factors should multiple to one
	f <- f/exp(mean(log(f)))

#	Output
	f
}


.calcFactorRLE <- function (data)
{
	gm <- exp(rowMeans(log(data)))
	apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
{
#	i <- apply(data<=0,1,all)
#	if(any(i)) data <- data[!i,,drop=FALSE]
	y <- t(t(data)/lib.size)
	f <- apply(y,2,function(x) quantile(x,p=p))
#	f/exp(mean(log(f)))
}

.calcFactorWeighted <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
#	TMM between two libraries
{
	obs <- as.numeric(obs)
	ref <- as.numeric(ref)

	if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
	if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

	logR <- log2((obs/nO)/(ref/nR))			# log ratio of expression, accounting for library size
	absE <- (log2(obs/nO) + log2(ref/nR))/2	# absolute expression
	v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref	 # estimated asymptotic variance

#	remove infinite values, cutoff based on A
	fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

	logR <- logR[fin]
	absE <- absE[fin]
	v <- v[fin]

	if(max(abs(logR)) < 1e-6) return(1)

#	taken from the original mean() function
	n <- length(logR)
	loL <- floor(n * logratioTrim) + 1
	hiL <- n + 1 - loL
	loS <- floor(n * sumTrim) + 1
	hiS <- n + 1 - loS

#	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
#	a fix from leonardo ivan almonacid cardenas, since rank() can return
#	non-integer values when there are a lot of ties
	keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

	if(doWeighting)
		f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
	else
		f <- mean(logR[keep], na.rm=TRUE)

#	Results will be missing if the two libraries share no features with positive counts
#	In this case, return unity
	if(is.na(f)) f <- 0
	2^f
}

