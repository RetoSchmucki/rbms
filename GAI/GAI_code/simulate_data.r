
# This file contains code to simulate example data to demonstrate the GAI model approach

# To fit to real data, the model requires fitting to data for each separately. Code for bootstrapping, regressing parameters or the hierarchical approach is available on request.

# Set working directory, containing model_code.r 
 setwd(...)
# Load file containing functions to fit the model
source(file="GAI/GAI_code/model_code.r")

###########################################
# We simulate data for the P/N2 case (for a single year)
###########################################
# Number of sites
nS <- 500
# Number of visits
nT <- 25

# Set parameter values within {a}
mu1 <- 8
mud <- 10
w <- .6
sigma <- 2.5
# Relative abundance values
N <- rpois(nS,150)
# Simulate data
y <- matrix(rpois(nS*nT,rep(N,each=nT)*(w*dnorm(1:nT,mu1,sigma) + (1-w)*dnorm(1:nT,mu1+mud,sigma))),nrow=nS,ncol=nT,byrow=TRUE)
# Create missing values
y[sample(1:(nS*nT),.3*nS*nT)] <- NA

plot(1:25, y[1,], ylim = c(0, (2*max(y[1,], na.rm = TRUE))))
for (i in 2:nS){
    points(1:25, y[i, ])
}

##########################################
# Model fitting
##########################################
# dist.type can be "P", "NB" or "ZIP"b  (e.g. Poisson, Negative Binomial, Zero-Inflated Poisson)
dist.type <- "NB"
# a.type can be "N", "SO" or "S" (e.g Normal, Stop-Over, Spline)
a.type <- "S"
# If a.type is "N" or "SO"
B <- 1
mu.type <- "common"
mu.diff.type <- "common"
sigma.type <- "hom"
w.type <- "common"

if( a.type == "S"){
    degf <- 12
    deg <- 3
}

# Specify number of random starts 
nstart <- 3
# Output is a list of length 3: the best model of multiple starts, output from the multiple starts, and the log-likelihood values for the multiple starts
output <- fit_it_model()

#debugonce(fit_it_model)
# To estimate the index value G (for a single year)
mean(output[[1]]$N.est)

output[[1]]$afunc.out[1,]
points(1:25, output[[1]]$Fitted[1,], col = 'blue', type = 'l')

cy <- data.frame(count =as.vector(t(y)), week =rep(1:nT, nS), site = rep(1:nS, each = nT))

library(mgcv)

m1 <- gam(count ~ s(week, bs = "cr") + factor(site) -1, data = cy, family = "nb")
cy$fitted <- predict(m1, newdata = cy, type = "response")

points(1:nT, cy$fitted[cy$site == 1], col = 'magenta', type = 'l')

m2 <- gam(count ~ s(week, bs = "cr") - 1, data = cy, family = "nb")
cy$fitted2 <- predict(m2, newdata = cy, type = "response")

points(1:nT, cy$fitted2[cy$site == 1], col = "brown", type = "l")


pt <- proc.time()
m1 <- gam(count ~ s(week, bs = "cr") + factor(site) -1, data = cy, family = "nb")
cy$fitted <- predict(m1, newdata = cy, type = "response")
proc.time()-pt

pt <- proc.time()
output <- fit_it_model()
proc.time() - pt


sum(cy$fitted)/nS
sum(cy$fitted2)/nS