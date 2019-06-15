#
# Read dat from local file Epildata.txt and define variables
#
library(DPpackage)
Epildata <- dget("Epildata.txt")
n <- Epildata$N
y <- t(matrix(c(Epildata$y),ncol=59))
x1 <- Epildata$Trt
x2 <- log(Epildata$Base/4)
x3 <- log(Epildata$Age)
x4 <- Epildata$V4
#
# Reformat data to use with DPpackage
#
yy <- c(t(y))
xx1 <- NULL
xx2 <- NULL
xx3 <- NULL
xx4 <- NULL
id <- NULL
for (i in 1:length(x1)) {
  xx1 <- c(xx1,rep(x1[i],4))
  xx2 <- c(xx2,rep(x2[i],4))
  xx3 <- c(xx3,rep(x3[i],4))
  xx4 <- c(xx4,c(0,0,0,1))
  id <- c(id,rep(i,4))
}
xx5 <- xx1*xx2
#
# Define hyperparameters
#
beta0 <- rep(0,5)
Sbeta0 <- diag(1000,5)
tinv <- diag(1)
prior <- list(a0=2,b0=0.1,nu0=3,tinv=tinv,mub=rep(0,1),Sb=diag(1000,1),
              beta0=beta0,Sbeta0=Sbeta0)
state <- NULL
#
# Now specify MCMC stuff
#
nburn <- 200000
nsave <- 3000
nskip <- 200
ndisplay <- 10
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
#
# Now call DPglmm
#
fit.DPglmm <- DPglmm(fixed=yy~xx1+xx2+xx3+xx4+xx5,random=~1|id,
                     family=poisson(log),prior=prior,mcmc=mcmc,
                     state=state,status=TRUE)
#
# Please see DPpackage manual to learn how to extract components for inference.
#

summart(fit.DPglmm)
