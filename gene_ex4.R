#### ---- Gene expression ----
#### ---- Task: Analysis using DPM ----
#### ---- Date: 2019 May 04 (Sat) ----

# MODEL:
## y_i | th_i ~ N(th_i, sig^2)
## th_i ~ G 
## G ~ DP(M, G0) with G0 = N(0, 4)

# Hyperprior:
## 1/sig ~ Ga(a, b), a=1, b=1

# Inference with DPpackage

library(DPpackage)

a <- 1; b <- 1 # 1/sig ~ Ga(a, b)
m0 <- -3; B0 <- 4 # G0 = N(m0, B0)
M <- 1
H <- 10

## data: EIG 121 data
X <- read.table("EIG.txt", header = 1)
head(X, n=20)

y <- X$recorded
y <- log(y[!is.na(y)]) 
n <- length(y); n

# DPdensity
state <- NULL
nburn <- 10; nsave <- 1000; nskip <- 10; ndisplay <- 100 ## MCMC parameters
mcmc <- list(nburn = nburn, nsave = nsave, nskip = nskip, ndisplay = ndisplay)

## Fixing alpha, m1, and Psi1
prior1 <- list(alpha = 1
	, m1 = rep(0, 1)
	, psiinv1 = diag(0.5, 1)
	, nu1 = 4
	, tau1 = 1
	, tau2 = 100
)

prior4 <- list(a0 = 2
	, b0 = 1
	, m2 = rep(0, 1)
	, s2 = diag(100000, 1)
	, psiinv1 = solve(diag(0.5, 1))
	, nu1 = 4
	, tau1 = 1
	, tau2 = 100
)

## Fit models

fit1 <- DPdensity(y = y
	, prior = prior1
	, mcmc = mcmc
	, state = state
	, status = TRUE
)

summary(fit1)

fit4 <- DPdensity(y = y
	, prior = prior4
	, mcmc = mcmc
	, state = state
	, status = TRUE
)
summary(fit4)

## Plot estimated densities
plot(fit1$x1, fit1$dens, type = "l")
plot(fit4$x1, fit4$dens, type = "l")

## Plot the parameters
plot(fit1, output = "param", nfigr = 2, nfigc = 2)
plot(fit4, output = "param", nfigr = 2, nfigc = 2)

## Plot the posterior for specific parameters
plot(fit4, output = "param" ,param = "ncluster", nfigr = 1, nfigc = 2)

## Extracting the posterior mean of the specific means and covariance matrices
th <- DPrandom(fit1)[[1]]
th

## Ploting predictive information about the specific means and
## covariance matrices with HPD and Credibility intervals

plot(DPrandom(fit1,predictive=TRUE), ask = FALSE, hpd = TRUE)


