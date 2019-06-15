library(MASS)

# Tring to implement MCMC algorithm for two categorical variables

y1_beta0 <- 0
y1_beta1 <- 5
y1_trueSD <- 5
y2_beta0 <- 2
y2_beta1 <- 3
y2_trueSD <- 5
sampleSize <- 100

## Generate data
#x <- rnorm(sampleSize)
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y1 <- rbinom(sampleSize, 1, plogis(y1_beta0 + y1_beta1 * x + rnorm(n = sampleSize, sd = y1_trueSD)))
y2 <- rbinom(sampleSize, 1, plogis(y2_beta0 + y2_beta1 * x + rnorm(n = sampleSize, sd = y2_trueSD)))

Y <- cbind(y1, y2)
X <- cbind(rep(1, sampleSize), x)

## MCMC implementation
s <- ncol(Y)
n <- nrow(Y)
p <- ncol(X)

## Starting values

### Overall sd
### Intercepts
y1_beta0 <- 3
y1_beta1 <- 2
y2_beta0 <- y1_beta0
y2_beta1 <- y1_beta1

betas <- matrix(
	c(y1_beta0
		, y1_beta1
 		, y2_beta0
 		, y2_beta1
	), 2, 2
)

### Constants and latent vars for logistic approximation
nu <- 7.3				# df
sigma <- pi * sqrt((nu - 2)/(3 * nu))
phis <- rep(1, n)

### Constract correlation matrix
cor_b1b2 <- 0.20
R <- matrix(
	c(1, cor_b1b2
		, cor_b1b2, 1
	), 2, 2
)


### Compute conditional distibution of zi
i <- 1
XB <- X[i,] %*% betas
XB
zi <- rnorm(XB, sigma^2/phis[i] * R)
zi
### Sample from the full conditional distribution of phi
phis <- rgamma((nu + p)/2, (nu + sigma^(-0.5) * as.vector(t(zi - XB)) * ((R)^(-1)) * as.vector(zi - XB))/2)
phis
quit()
## Compute likelihoods by sampling from the probability densities
likelihood <- function(params){
	beta0 <- params[1]
	beta1 <- params[2]
	sd <- params[3]

	# Predictions
	preds <- beta0 + beta1 * x
	logll <- dnorm(y, mean = preds, sd = sd, log = TRUE)
	sumlogll <- sum(logll)
	return(sumlogll)
}

likelihood(c(0, 5, 5))

## Plot loglikehood profiles example
beta1_vals <- seq(3, 7, by = 0.05)
slopellhoods <- lapply(beta1_vals, function(b){
		return(likelihood(c(beta0, b, trueSD)))
	}
)

plot(data.frame(beta = beta1_vals, ll = unlist(slopellhoods)), type = "l", xlab = "Beta1", y = "Log likelihood")

## Define priopr distribution
prior <- function(params){
	beta0 <- params[1]
	beta1 <- params[2]
	sd <- params[3]
	beta0_prior <- dnorm(beta0, sd = 5, log = TRUE)
	beta1_prior <- dunif(beta1, min = 0, max = 10, log = TRUE)
	sd_prior <- dunif(sd, min = 0, max = 30, log = TRUE)
	return(sum(beta0_prior, beta1_prior, sd_prior))
}
prior(c(0, 5, 10))

## Define posterior
posterior <- function(params){
	return(likelihood(params) + prior(params))
}
posterior(c(0, 5, 10))

## Metropolis Algorithm
target <- function(params){
	return(rnorm(3, mean = params, sd = c(0.5, 0.1, 0.3)))
}

metropolisMCMC <- function(startvalue, iters){
	chain <- array(dim = c(iters + 1, 3))
	chain[1, ] = startvalue
	for (i in 1:iters){
		proposal <- target(chain[i, ])
		prob <- exp(posterior(proposal) - posterior(chain[i, ]))
		if (runif(1) < prob){
			chain[i + 1, ] <- proposal
		}else{
			chain[i+1, ] <- chain[i, ]
		}
	}
	return(chain)
}

