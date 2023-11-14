library(nimble)
library(BayesNSGP)
nimbleOptions(verbose = FALSE)

########### MDR example: comparing computational times
Nvec <- c(50,100,200,500,1000,2000,5000,10000)
tau_knot_coords  <- as.matrix(expand.grid(seq(0,1,length = 10),seq(0,1,length = 10)))
fullGPtime <- NNGPtime <- SGVtime <- rep(NA, length(Nvec))
for(n in 1:length(Nvec)){
  N <- Nvec[n]
  print(N)
  coords <- matrix(5*runif(2*N), ncol = 2)
  Xmat1 <- cbind(rep(1,N),coords[,1])
  Xmat2 <- cbind(rep(1,N),coords[,2])
  constants <- list( X_sigma = Xmat1, X_Sigma = Xmat2, X_mu = Xmat1, tau_knot_coords = tau_knot_coords, k = 10 )
  
  # fullGP
  if(N <= 1000){
    Rmodel_fullGP <- nsgpModel( likelihood = "fullGP", 
                                sigma_model = "logLinReg", Sigma_model = "compRegIso",
                                mu_model = "linReg", tau_model = "approxGP", 
                                constants = constants, coords = coords, data = rnorm(N) ) #, returnModelComponents = TRUE )
    fullGPtime[n] <- system.time(Rmodel_fullGP$calculate())[3]
    rm(Rmodel_fullGP)
  }
  
  # NNGP
  Rmodel_NNGP <- nsgpModel( likelihood = "NNGP", 
                            sigma_model = "logLinReg", Sigma_model = "compRegIso",
                            mu_model = "linReg", tau_model = "approxGP", 
                            constants = constants, coords = coords, data = rnorm(N) )
  NNGPtime[n] <- system.time(Rmodel_NNGP$calculate())[3]
  rm(Rmodel_NNGP)
  
  # SGV
  Rmodel_SGV <- nsgpModel( likelihood = "SGV", 
                           sigma_model = "logLinReg", Sigma_model = "compRegIso",
                           mu_model = "linReg", tau_model = "approxGP", 
                           constants = constants, coords = coords, data = rnorm(N) )
  SGVtime[n] <- system.time(Rmodel_SGV$calculate())[3]
  rm(Rmodel_SGV)
  
}

timeDF <- data.frame(
  N = Nvec,
  exactGP = fullGPtime,
  SGV = SGVtime,
  NNGP = NNGPtime
)
#       N exactGP   SGV  NNGP
# 1    50   0.034 0.070 0.060
# 2   100   0.045 0.081 0.061
# 3   200   0.085 0.123 0.078
# 4   500   0.918 0.242 0.130
# 5  1000   3.965 0.421 0.209
# 6  2000      NA 0.969 0.370
# 7  5000      NA 3.892 0.965
# 8 10000      NA 7.780 1.724

plot(Nvec, log(fullGPtime), type = "b", pch = "+")
lines(Nvec, log(NNGPtime), type = "b", pch = "+", col = 2)
lines(Nvec, log(SGVtime), type = "b", pch = "+", col = 4)

plotDF <- data.frame(
  N = rep(Nvec, 3),
  Time = c(fullGPtime, SGVtime, NNGPtime),
  Likelihood = rep(c("Exact GP", "SGV", "NNGP"), each = length(Nvec))
)
library(ggplot2)
ggplot(plotDF, aes(x = N, y = Time, color = Likelihood)) + geom_point() + geom_line() +
  scale_y_continuous(trans="log", breaks=10^(-2:3), limits = c(0.01,100), name = "Time (s)") +
  scale_x_continuous(breaks = c(50,1000,2000,5000,10000))

ggsave(
  "Timings.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)




