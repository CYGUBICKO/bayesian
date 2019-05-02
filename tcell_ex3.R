#### ---- T-cell receptor simulation ----
#### ---- Task: Analysis using DP ----
#### ---- Date: 2019 May 01 (Wed) ----

# The dataset contains counts of Tcell receptors 

## Sampling model
# p(yi|F) = G(yi), yi>0

## Prior
# F = DP(M G0)

#######################################################################
## Without censoring -- DP for y[i], y>0 only
#######################################################################

## Data

yf <- c(37,11,5,2) # frequencies
xf <- c(1,2,3,4)   # counts
y <- rep(xf,yf)    # raw data
n <- length(y)
k <- n             # initialize

## hyperparameters
## a0 <- 1; b0 <- 1   # hyperprior b ~ Ga(a0,b0)
a <- 1; b <- .05     # G0(mu) = Ga(a,b)
lambda <- 300        # k ~ Poi(lambda)

M <- 1
H <- 10
N <- 25; p=8

# Genrate x ~ Dir(a,...a)(n-dim)
rdiric <- function(n, a){
	p <- length(a)
	m <- matrix(nrow = n, ncol = p)
	for (i in 1:p){
		m[,i] <- rgamma(n, a[[i]])
	}
	sumvec <- m %*% rep(1,p)
   m / as.vector(sumvec)
}


cdfplt <- function(x,y,lo=NULL,hi=NULL,
                    lw=1,lt=1,cl=1)
  {
    p <- length(x)
    if (!is.null(hi)) # final line segment
      if(hi > x[p])
        lines(c(x[p],hi), c(1,1),col=cl,lwd=lw,lty=lt)
    if (!is.null(lo)){ # initial line segment
      if (lo<x[1]){
        lines(c(lo,x[1]), c(0,0),col=cl,lwd=lw,lty=lt)
        if (y[1]>0)   # prob mass at x[1]
          points(x[1],0,pch=1)
      }
    }# initial seg
    ylag <- c(0,y)    # prev prob
    for(i in 1:p){
      if (i<p)
        lines(x[c(i,i+1)],y[c(i,i)],col=cl,lwd=lw,lty=lt)
      if (y[i]>ylag[i]){ # prob mass @ x[i]
        points(x[i],y[i],pch=19)
        points(x[i],ylag[i],pch=1)
      }
    }# for i
  }

# Generate posterior p(G|y) for yi ~ G
## G0: prior mean
## Ghat: empirical mean
## Gbar: posterior mean
## G: posterior sample

sampleDP <- function(N = 10
	, p = 8
	, plt.G0 = FALSE
	, plt.spl = FALSE
	, plt.Ghat = FALSE
	, cdf = TRUE
	, ltype = "n"){
	xgrid <- 1:p
	r <- 1/(1 - dpois(0, lambda = 2))
	G0 <- dpois(xgrid, lambda = 2) * r	# Prior base measure
	G1 <- M * G0								# Posterior base measure
	G1[xf] <- G1[xf] + yf
	G <- rdiric(N, G1)
	Gcdf <- apply(G, 1, cumsum)
	n <- sum(yf)
	Gbar <- G1/(n + M)
	Gbarcdf <- cumsum(Gbar)
	
	if (cdf){
		matplot(xgrid
			, Gcdf
			, type = ltype
			, bty = "l"
			, xlim = c(1, 10)
			, ylim = c(0, 1)
			, xlab = "X"
			, ylab = "G"
		)
	} else{
		matplot(xgrid
			, t(G)
			, type = ltype
			, bty = "l"
			, xlim = c(1, 8)
			, ylim = c(0, 1)
			, xlab = "COUNT"
			, ylab = "G"
		)		
	}

	if (plt.spl){
		for (i in 1:N){
			if (cdf)
				cdfplt(xgrid, Gcdf[, i], lw = 1, hi = 10)
			else
				lines(xgrid, G[i,], lw = 1, lt = i)
		}
	}

	G0cdf <- cumsum(G0)
	if (plt.G0){
		if (cdf)
			cdfplt(xgrid, G0cdf, lt = 3, lw = 3, cl = 1)
		else
			lines(xgrid, G0, lt = 3, lw = 3, col = 1)
	}

	Ghat <- rep(0, p)
	Ghat[xf] <- yf/n
	Ghatcdf <- cumsum(Ghat)
	if (plt.Ghat){
		if (cdf){
			cdfplt(xgrid, Ghatcdf, lt = 1, lw = 3, cl = 1)
			cdfplt(xgrid, Gbarcdf, lt = 2, lw = 3, cl = 1)
		} else{
			xg <- as.numeric(names(table(y)))
			lines(table(y)/n, lwd = 1)
			points(xg, table(y)/n, pch = 19)
			lines(xgrid, Gbar, lty = 1, lwd = 3)
		}
	}# plt.Ghat

	list(r = r
		, G0 = G0
		, G1 = G1
		, G = G
		, Gcdf = Gcdf
		, Gbar = Gbar
		, G0cdf = G0cdf
	)
}
yf
xf
sampleDP(ltype = "l")
sampleDP(cdf = FALSE, ltype = "l")
sampleDP(plt.spl = TRUE)
sampleDP(plt.G0 = TRUE)
sampleDP(plt.Ghat = TRUE)

