# This script describes step by step implementation of Metropolis-Hastings MCMC

beta0 <- 0
beta1 <- 5
trueSD <- 5
sampleSize <- 100

## Generate data
#x <- rnorm(sampleSize)
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y <- beta0 + beta1 * x + rnorm(n = sampleSize, sd = trueSD)
plot(x, y, main = "Test data")

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

startvalue <- c(0, 4, 10)
chain <- metropolisMCMC(startvalue, 100000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
acceptance

## Plot
par(mfrow = c(2,3))
hist(chain[-(1:burnIn),2],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = beta1, col="red" )
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = beta0, col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]) )
abline(v = trueSD, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = beta1, col="red" )
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = beta0, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSD, col="red" )

# for comparison:
summary(lm(y~x))


