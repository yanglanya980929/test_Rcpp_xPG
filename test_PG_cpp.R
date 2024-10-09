install.packages("Rcpp")
library(Rcpp)
library(devtools)
# remove.packages("dlm")
devtools::install_local("/home/yangl35/STOR-601-env/PGcpp/dlm", force = TRUE)
library(dlm)
library(mvtnorm)
library(ggplot2)
library(MASS)
library(coda)

source("/home/yangl35/STOR-601-env/PGcpp/rest_of_dlm.R")

setwd("/home/yangl35/STOR-601-env/ParticleGibbs")
source("sampling_functions.R")
source("utils.R")
source("x_particle_Gibbs.R")
source("x_conditional_particle_filter.R")
source("particle_Gibbs.R")
source("conditional_particle_filter.R")

################ dynamic linear model #######################
d <- 1
M <- 250
N <- 100

# parameters
v0 <- 1 # starting from a point far away from zero
vx <- 1
vy.low <- 0.25
vy.high <- 4
d <- 1
a <- 0.8
b <- 0

# arguments 
simX0_params_2 <- c(v0)
simXt_params_2 <- c(a, b, vx)
loglike_params_2_low <- c(vy.low)
loglike_params_2_high <- c(vy.high)

ts <- seq(from = 1, to = N, by = 1)
xs.2 <- rep(NA, N)
set.seed(123456)
xs.2[1] <- simX0_2(simX0_params_2, 1, d)
print(xs.2[1])
for(i in 1:(N-1)){
  xs.2[i+1] <- simXt_2(xs.2[i], simXt_params_2, 1, d)
}

# simulate two observed processes with low and high variance
set.seed(123456)
ys.2.low <- rnorm(N, mean = xs.2, sd = sqrt(vy.low))
set.seed(123456)
ys.2.high <- rnorm(N, mean = xs.2, sd = sqrt(vy.high))

plot(ts, ys.2.high, type="b",pch=1,lwd=2,xlab="t",ylab="x",col="orange") # high
lines(ts, xs.2, type="b",pch=4,lwd=2,xlab="t",ylab="x",col="blue") # hidden process
lines(ts, ys.2.low, type="b",pch=1,lwd=2,xlab="t",ylab="x",col="red") # low
abline(h=0,col="black",lty=2)
legend("bottomleft",c("Yt.low","Yt.high","Xt"),col=c("red","orange","blue"), lty = c(1, 1), lwd = c(2, 2))

#################################################################
path0 <- matrix(xs.2, nrow = 1, ncol = N)
num <- 10000
a0 <- 0.2
eps <- 0.6
simXt_params_2 <- c(a0, b, vx) 

set.seed(123456)
# timing 
start_time <- Sys.time() 
xPG.1 <- xParticleGibbs(simPara.a, Xproposal_0, Xproposal_t, loglike_2, simX0_params_2, simXt_params_2, loglike_params_2_low, ys.2.low, path0, M, num, g_inverse_0, g_inverse_1, eps)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)
# Time difference of 2.689357 mins

xTheta.1 <- xPG.1$Theta
xPaths.1 <- xPG.1$Paths

plot(xTheta.1, type = "l", xlab = "Iteration", ylab = "Theta", main = "Trace Plot of Theta") # theta
plot(xPaths.1[,N,], type = "l", xlab = "Iteration", ylab = "X_T", main = "Trace Plot of X_T") # X_T
plot(xPaths.1[,N/2,], type = "l", xlab = "Iteration", ylab = "X_1", main = "Trace Plot of X_{N/2}")
plot(xPaths.1[,1,], type = "l", xlab = "Iteration", ylab = "X_1", main = "Trace Plot of X_1") # X_1
# plot(xPaths.1[,50,], type = "l", xlab = "Iteration", ylab = "X_{T/2}", main = "Trace Plot of X_{T/2}") # X_{T/2}
B = 10
a.hat <- mean(xTheta.1[B:length(xTheta.1)])
hist(xTheta.1[B:length(xTheta.1)])
quantile(xTheta.1[B:length(xTheta.1)],probs = c(0.05,0.5,0.95))
# 5%       50%       95% 
# 0.6696280 0.7796335 0.8828776  
