#  Estimating dispersion using weighted likelihood empirical Bayes.

estimateDisp <- function(y, ...)
UseMethod("estimateDisp")

estimateDisp.DGEList <- function(y, design=NULL, prior.df=NULL, trend.method="locfit", mixed.df=FALSE, tagwise=TRUE, span=NULL, min.row.sum=5, grid.length=21, grid.range=c(-10,10), robust=FALSE, winsor.tail.p=c(0.05,0.1), tol=1e-06, ...)
#  Yunshun Chen. Created 16 March 2016.
{
	y <- validDGEList(y)
	group <- y$samples$group
	lib.size <- y$samples$lib.size * y$samples$norm.factors

	d <- estimateDisp(y=y$counts, design=design, group=group, lib.size=lib.size, offset=getOffset(y), prior.df=prior.df, trend.method=trend.method, mixed.df=mixed.df, tagwise=tagwise, span=span, min.row.sum=min.row.sum, grid.length=grid.length, grid.range=grid.range, robust=robust, winsor.tail.p=winsor.tail.p, tol=tol, weights=y$weights, ...)

	if(!is.null(design)) y$design <- design
	y$common.dispersion <- d$common.dispersion
	y$trended.dispersion <- d$trended.dispersion
	if(tagwise) y$tagwise.dispersion <- d$tagwise.dispersion
	y$AveLogCPM <- aveLogCPM(y)
	y$trend.method <- trend.method
	y$prior.df <- d$prior.df
	y$prior.n <- d$prior.n
	y$span <- d$span
	y
}

estimateDisp.default <- function(y, design=NULL, group=NULL, lib.size=NULL, offset=NULL, prior.df=NULL, trend.method="locfit", mixed.df=FALSE, tagwise=TRUE, span=NULL, min.row.sum=5, grid.length=21, grid.range=c(-10,10), robust=FALSE, winsor.tail.p=c(0.05,0.1), tol=1e-06, weights=NULL, ...)
#  Use GLM approach if a design matrix is given, and classic approach otherwise.
#  It calculates a matrix of likelihoods for each gene at a set of dispersion grid points, and then calls WLEB() to do the shrinkage.
#  Yunshun Chen, Aaron Lun, Gordon Smyth. Created July 2012. Last modified 03 Oct 2016.
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	if(ntags==0) return(numeric(0))
	nlibs <- ncol(y)
	
#	Check trend.method
	trend.method <- match.arg(trend.method, c("none", "loess", "locfit", "movingave"))
	
#	Check group
	if(is.null(group)) group <- rep(1, nlibs)
	if(length(group)!=nlibs) stop("Incorrect length of group.")
	group <- dropEmptyLevels(group)

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(y)
	if(length(lib.size)!=nlibs) stop("Incorrect length of lib.size.")
	
#	Check offset
	offset <- .compressOffsets(y, lib.size=lib.size, offset=offset)

#	Check weights
	weights <- .compressWeights(weights)

#	Check for genes with small counts
	sel <- rowSums(y) >= min.row.sum
	sely <- .subsetMatrixWithoutCopying(y, i=sel)
	seloffset <- .subsetMatrixWithoutCopying(offset, i=sel)
	selweights <- .subsetMatrixWithoutCopying(weights, i=sel)
	
#	Spline points
	spline.pts <- seq(from=grid.range[1], to=grid.range[2], length=grid.length)
	spline.disp <- 0.1 * 2^spline.pts
	grid.vals <- spline.disp/(1+spline.disp)
	l0 <- matrix(0, sum(sel), grid.length)

#	Classic edgeR
	if(is.null(design)){
		# One way
		cat("Design matrix not provided. Switch to the classic mode.\n")
		if(length(levels(group))==1)
			design <- matrix(1, nlibs, 1)
		else
			design <- model.matrix(~group)
		if( all(tabulate(group)<=1) ) {
			warning("There is no replication, setting dispersion to NA.")
			return(list(common.dispersion=NA, trended.dispersion=NA, tagwise.dispersion=NA))
		}
		
		eq <- equalizeLibSizes(y, group=group, dispersion=0.01, lib.size=lib.size)
		y.pseudo <- eq$pseudo.counts[sel, , drop=FALSE]
		y.split <- splitIntoGroups(y.pseudo, group=group)
		delta <- optimize(commonCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y.split, der=0)
		delta <- delta$maximum
		disp <- delta/(1-delta)

		eq <- equalizeLibSizes(y, group=group, dispersion=disp, lib.size=lib.size)
		y.pseudo <- eq$pseudo.counts[sel, , drop=FALSE]
		y.split <- splitIntoGroups(y.pseudo, group=group)
	
		for(j in 1:grid.length) for(i in 1:length(y.split)) 
			l0[,j] <- condLogLikDerDelta(y.split[[i]], grid.vals[j], der=0) + l0[,j]
	}
	# GLM edgeR 
	else {
		design <- as.matrix(design)
		if(ncol(design) >= nlibs) {
			warning("No residual df: setting dispersion to NA")
			return(list(common.dispersion=NA, trended.dispersion=NA, tagwise.dispersion=NA))
		}

		# Identify which observations have means of zero (weights aren't needed here).
		glmfit <- glmFit(sely, design, offset=seloffset, dispersion=0.05, prior.count=0)
		zerofit <- .areFittedValuesZero(y=glmfit$counts, mu=glmfit$fitted.values)
		by.group <- .comboGroups(zerofit)

		for (subg in by.group) { 
			cur.nzero <- !zerofit[subg[1],]
			if (!any(cur.nzero)) { next } 

			# Removing samples with zero means from design, to avoid attempts to converge to -Inf.
			if (all(cur.nzero)) { 
				redesign <- design
			} else {
				redesign <- design[cur.nzero,,drop=FALSE]
				QR <- qr(redesign)
				redesign <- redesign[,QR$pivot[1:QR$rank],drop=FALSE]
				if (nrow(redesign) == ncol(redesign)) { next }
			}

			cury <- .subsetMatrixWithoutCopying(sely, i=subg, j=cur.nzero)
			curo <- .subsetMatrixWithoutCopying(seloffset, i=subg, j=cur.nzero)
			curw <- .subsetMatrixWithoutCopying(selweights, i=subg, j=cur.nzero)
		   
			# Using the last fit to hot-start the next fit
			last.beta <- NULL
			for(i in seq_len(grid.length)) {
				out <- adjustedProfileLik(spline.disp[i], y=cury, design=redesign, 
					offset=curo, weights=curw, start=last.beta, get.coef=TRUE)
				l0[subg,i] <- out$apl
				last.beta <- out$beta
			}
		}
	}

	overall <- maximizeInterpolant(spline.pts, matrix(colSums(l0), nrow=1))
	common.dispersion <- 0.1 * 2^overall
	AveLogCPM <- aveLogCPM(y, lib.size=lib.size, dispersion=common.dispersion, weights=weights)

	out.1 <- WLEB(theta=spline.pts, loglik=l0, covariate=AveLogCPM[sel], trend.method=trend.method, 
		mixed.df=mixed.df, span=span, overall=FALSE, individual=FALSE, m0.out=TRUE)
	disp.trend <- 0.1 * 2^out.1$trend
	trended.dispersion <- rep( disp.trend[which.min(AveLogCPM[sel])], ntags )
	trended.dispersion[sel] <- disp.trend

	if(!tagwise) return(list(common.dispersion=common.dispersion, trended.dispersion=trended.dispersion))

	# Calculate prior.df
	if(is.null(prior.df)){
		glmfit <- glmFit(sely, offset=seloffset, weight=selweights, design=design, dispersion=disp.trend, prior.count=0)

		# Residual deviances
		df.residual <- glmfit$df.residual

		# Adjust df.residual for fitted values at zero
        zerofit <- .areFittedValuesZero(y=glmfit$counts, mu=glmfit$fitted.values)
		df.residual <- .residDF(zerofit, design)

		# Empirical Bayes squeezing of the quasi-likelihood variance factors
		s2 <- glmfit$deviance / df.residual
		s2[df.residual==0] <- 0
		s2 <- pmax(s2,0)
		s2.fit <- squeezeVar(s2, df=df.residual, covariate=AveLogCPM[sel], robust=robust, winsor.tail.p=winsor.tail.p)
		prior.df <- s2.fit$df.prior
	}
	ncoefs <- ncol(design)
	prior.n <- prior.df/(nlibs-ncoefs)

	if (trend.method!='none') { 
		tagwise.dispersion <- trended.dispersion
	} else {
		tagwise.dispersion <- rep(common.dispersion, ntags)
	}

	# Checking if the shrinkage is near-infinite.
	too.large <- prior.n > 1e6
	if (!all(too.large)) { 
		temp.n <- prior.n
		if (any(too.large)) { 
			temp.n[too.large] <- 1e6 
		}

		# Estimating tagwise dispersions	
		out.2 <- WLEB(theta=spline.pts, loglik=l0, prior.n=temp.n, covariate=AveLogCPM[sel], 
			trend.method=trend.method, mixed.df=mixed.df, span=span, overall=FALSE, trend=FALSE, m0=out.1$shared.loglik)

		if (!robust) { 
			tagwise.dispersion[sel] <- 0.1 * 2^out.2$individual
		} else {
			tagwise.dispersion[sel][!too.large] <- 0.1 * 2^out.2$individual[!too.large]
		}
	}

	if(robust) {
		temp.df <- prior.df
		temp.n <- prior.n
		prior.df <- prior.n <- rep(Inf, ntags)
		prior.df[sel] <- temp.df
		prior.n[sel] <- temp.n
	}

	list(common.dispersion=common.dispersion, trended.dispersion=trended.dispersion, tagwise.dispersion=tagwise.dispersion, span=out.1$span, prior.df=prior.df, prior.n=prior.n)
}


WLEB <- function(theta, loglik, prior.n=5, covariate=NULL, trend.method="locfit", mixed.df=FALSE, span=NULL, 
	overall=TRUE, trend=TRUE, individual=TRUE, m0=NULL, m0.out=FALSE)
#  Weighted likelihood empirical Bayes for estimating a parameter vector theta
#  given log-likelihood values on a grid of theta values
#  Yunshun Chen, Gordon Smyth
#	Created July 2012. Last modified 24 October 2012.
{
#	Check loglik
	loglik <- as.matrix(loglik)
	ntheta <- ncol(loglik)
	ntags <- nrow(loglik)

#	Check covariate and trend
	if(is.null(covariate))
		trend.method <- "none"
	else
		trend.method <- match.arg(trend.method, c("none", "loess", "locfit", "movingave"))

	if(trend.method=="locfit" & mixed.df) trend.method <- "locfit_mix"

#	Set span
	if(is.null(span)) if(ntags<=50) span <- 1 else span <- 0.25+0.75*(50/ntags)^0.5

#	Output	
	out <- list()
	out$span <- span

#	overall prior
	if(overall)
		out$overall <- maximizeInterpolant(theta, matrix(colSums(loglik), nrow=1))

#	trended prior
	if(is.null(m0))
	m0 <- switch(trend.method,
		"movingave" = {
			o <- order(covariate)
			oo <- order(o)
			movingAverageByCol(loglik[o,], width=floor(span*ntags))[oo,]
		},
		"loess" = loessByCol(loglik, covariate, span=span)$fitted.values,
		"locfit" = locfitByCol(loglik, covariate, span=span, degree=0),
		"locfit_mix" = {
			deg0 <- locfitByCol(loglik, covariate, span=span, degree=0)
			deg1 <- locfitByCol(loglik, covariate, span=span, degree=1)
			r <- range(covariate)
			w <- pbeta((covariate-r[1])/(r[2]-r[1]),shape1=2,shape2=2)
			w*deg0 + (1-w)*deg1
		},
		"none" = matrix(colMeans(loglik), ntags, length(theta), byrow=TRUE)
	)

#	make sure each row of m0 is unimodal
	if(trend.method=="locfit_mix"){
		for(i in ncol(m0):3){
			diff1_m0 <- m0[,i] - m0[,i-1]
			diff2_m0 <- m0[,i-1] - m0[,i-2]
			k <- which(diff1_m0>0 & diff2_m0<0)
			m0[k,1:(i-2)] <- m0[k,(i-1)]
		}
	}

	if(trend)
		out$trend <- maximizeInterpolant(theta, m0)

#	weighted empirical Bayes posterior estimates
	if(individual){
		stopifnot(all(is.finite(prior.n)))
		l0a <- loglik + prior.n*m0
		out$individual <- maximizeInterpolant(theta, l0a)
	}

	if(m0.out) out$shared.loglik <- m0

	out
}

.subsetMatrixWithoutCopying <- function(x, i, j) 
# This will attempt to subset the matrix without any copying if
# it detects that 'i' and 'j' don't modify the ordering of the matrix.
# This reduces the memory footprint for large matrices.
#
# written by Aaron Lun
# created 29 September 2016
{
	if (is(x, "compressedMatrix")) {
		# Subsetting is ignored if rows/columns are repeated.
		if (attributes(x)$repeat.row) i <- TRUE
		if (attributes(x)$repeat.col) j <- TRUE
	} else if (!is.matrix(x)) {
		stop("'x' must be a matrix for subsetting")
	}

	isokay <- TRUE
	if (!missing(i)) {
		# Most flexible way of handling different types of subset vectors;
		# try it out and see if it gives the same results.
		example <- cbind(seq_len(nrow(x)))
		rownames(example) <- rownames(x)
		if (!identical(example, example[i,,drop=FALSE])) isokay <- FALSE
	}
	if (!missing(j)) {
		example <- rbind(seq_len(ncol(x)))
		colnames(example) <- colnames(x)
		if (!identical(example, example[,j,drop=FALSE])) isokay <- FALSE
	}

	if (isokay) {
		# Avoids copying if no modification incurred.
		return(x)
	} else {
		return(x[i,j,drop=FALSE])
	}
}

.areFittedValuesZero <- function(y, mu, tol=1e-4) 
# A slightly more efficient way to check if the fitted values are zero,
# based on both the count and fitted value being zero for an observation.
{
    if (!is.double(mu)) storage.mode(mu) <- "double"
    out <- .Call(.cR_check_zero_fitted, y, mu, as.double(tol))
    if (is.character(out)) stop(out)
    return(out)
}
