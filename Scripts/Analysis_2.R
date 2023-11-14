#-----Libraries----

rm(list = ls())
library(BayesNSGP)
library(ggplot2)
library(readr)
library(dplyr)
library(maps)
nimbleOptions(verbose = FALSE) # Silences warnings

# Paper uses RW block samplers but package throws warnings

#---------Setup----

set.seed(333) # For replication purposes

# Load raw data
CONUS_WY2018 <- read_csv("Desktop/Spatial/Project/CONUS_WY2018.csv")
CONUS_predGrid <- read_csv("Desktop/Spatial/Project/CONUS_predGrid.csv")
CONUSprecip <- as.data.frame(CONUS_WY2018)

# Various constants and distance matrices
coords <- as.matrix(CONUSprecip[,c("longitude", "latitude")])
data <- CONUSprecip$logPR
CONUSprecip$Zelevation <- scale(CONUSprecip$Xelevation)
CONUSprecip$Zlongitude <- scale(CONUSprecip$longitude)
Xmat <- unname(lm(logPR ~ Zelevation*Zlongitude, x = TRUE,
                   data = CONUSprecip)$x) # Design Matrix

x_min <- min(CONUSprecip[,"longitude"])
x_max <- max(CONUSprecip[,"longitude"])
y_min <- min(CONUSprecip[,"latitude"])
y_max <- max(CONUSprecip[,"latitude"])
knot_coords <- expand.grid(
  lon = seq(from = x_min + 1, to = x_max - 2, length = 10),
  lat = seq(from = y_min + 2, to = y_max, length = 7)
  )
knot_coords <- as.matrix(knot_coords)
knot_coords <- knot_coords[-c(1:4,6:7,9:12,19:21,30,40,58,67:70),]
knot_coords <- as.data.frame(knot_coords)
knot_groups <- matrix(c(1,3:5,10:12,19:21,
                        2,6:8,13:16,22:25,
                        9,17:18,26:27,36:37,45:46,
                        28:30,38:40,47:49,
                        31:35,41:44,50),
                      ncol = 5)

N <- 1500 # Number of iterations (~ 2 hours for full code on 16GB M1 Macbook)

#---Exploratory----

US <- map_data("usa")

Knot_A <- knot_coords[knot_groups[,1],]
Knot_B <- knot_coords[knot_groups[,2],]
Knot_C <- knot_coords[knot_groups[,3],]
Knot_D <- knot_coords[knot_groups[,4],]
Knot_E <- knot_coords[knot_groups[,5],]

ggplot() +
  geom_polygon(data = US,
               aes(x=long, y=lat, group=group),
               fill = "white", colour = "black") + 
  geom_point(data = CONUSprecip,
             aes(x=longitude, y=latitude,color=logPR)) + xlab('') + ylab('') +
  geom_point(data = Knot_A, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 0, size = 3, alpha = 0.75) +
  geom_point(data = Knot_B, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 1, size = 3, alpha = 0.75) +
  geom_point(data = Knot_C, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 2, size = 3, alpha = 0.75) +
  geom_point(data = Knot_D, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 5, size = 3, alpha = 0.75) +
  geom_point(data = Knot_E, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 6, size = 3, alpha = 0.75) +
  scale_colour_viridis_c(begin = 0.3, direction = -1) + coord_map()

ggsave(
  "Precipitation_2018.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

ggplot() +
  geom_polygon(data = US,
               aes(x=long, y=lat, group=group),
               fill = "white", colour = "black") +
  geom_point(data = CONUSprecip,
             aes(x=longitude, y=latitude,color=elevation)) +
  scale_colour_gradientn(colours = terrain.colors(20)) +
  xlab('') + ylab('') + coord_map()

ggsave(
  "Elevation_2018.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

#--------Models----

# Constants

constants_A <- list(nu = 0.5, k = 15, tau_HP1 = 10,
                  sigma_knot_coords = knot_coords,
                  sigma_HP1 = 10, sigma_HP2 = 5,
                  sigma_HP3 = 10, sigma_HP4 = 10,
                  X_Sigma = Xmat, Sigma_HP1 = 5,
                  maxAnisoDist = max(dist(coords)),
                  X_mu = Xmat, mu_HP1 = 10 )

constants_B <- list(nu = 0.5, k = 15, tau_HP1 = 10, sigma_HP1 = 10,
                    Sigma_HP1 = max(dist(coords)), X_mu = Xmat, mu_HP1 = 10)

# Models

Rmodel_A <- nsgpModel(likelihood = "SGV", constants = constants_A,
                    coords = coords, data = data, 
                    tau_model = "constant", sigma_model = "approxGP",
                    mu_model = "linReg", Sigma_model = "compReg")

Rmodel_B <- nsgpModel(likelihood = "SGV", constants = constants_B,
                      coords = coords, data = data, 
                      tau_model = "constant", sigma_model = "constant",
                      mu_model = "linReg", Sigma_model = "constant")

# Configurations

conf_A <- configureMCMC(Rmodel_A)

conf_A$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
conf_A$removeSamplers(c("Sigma_coef1[1]","Sigma_coef1[2]",
                      "Sigma_coef2[1]","Sigma_coef2[2]",
                      "Sigma_coef3[1]","Sigma_coef3[2]"))
conf_A$removeSamplers(c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"))
conf_A$addSampler(target = c("beta[1]","beta[2]","beta[3]","beta[4]"),
                  type = "RW_block")
conf_A$addSampler(target = c("Sigma_coef1[1]","Sigma_coef1[2]",
                           "Sigma_coef2[1]", "Sigma_coef2[2]",
                           "Sigma_coef3[1]","Sigma_coef3[2]"),
                  type = "RW_block")
conf_A$addSampler(target = c("sigmaGP_mu","sigmaGP_phi","sigmaGP_sigma"),
                  type = "RW_block")
conf_A$removeSamplers("w_sigma[1:50]")
for(h in 1:ncol(knot_groups)){
  conf_A$addSampler(target = c(paste0("w_sigma[",knot_groups[,h],"]")),
                    type = "RW_block" )
}

conf_B <- configureMCMC(Rmodel_B)

conf_B$removeSamplers(c("beta[1]","beta[2]","beta[3]","beta[4]"))
conf_B$removeSamplers(c("Sigma_coef1","Sigma_coef2","Sigma_coef3"))
conf_B$addSampler(target = c("beta[1]","beta[2]","beta[3]","beta[4]"),
                  type = "RW_block")
conf_B$addSampler(target = c("Sigma_coef1",
                             "Sigma_coef2",
                             "Sigma_coef3"),
                  type = "RW_block")

Rmcmc_A <- buildMCMC(conf_A) # Build the MCMC
Cmodel_A <- compileNimble(Rmodel_A) # Compile the model
Cmcmc_A <- compileNimble(Rmcmc_A, project = Rmodel_A) # Compile the MCMC

Rmcmc_B <- buildMCMC(conf_B) # Build the MCMC
Cmodel_B <- compileNimble(Rmodel_B) # Compile the model
Cmcmc_B <- compileNimble(Rmcmc_B, project = Rmodel_B) # Compile the MCMC

samples_A <- runMCMC(Cmcmc_A, niter = N) # Run

write_csv(as.data.frame(samples_A),
          "/Users/sanketdutta/Desktop/Spatial/Project/Data/MCMC_A_2.csv")

samples_B <- runMCMC(Cmcmc_B, niter = N) # Run

write_csv(as.data.frame(samples_B),
          "/Users/sanketdutta/Desktop/Spatial/Project/Data/MCMC_B_2.csv")

#----Predictive----

predCoords <- as.data.frame(CONUS_predGrid)

# Define the grid size based on RAM usage
grid_size <- 0.80

predCoords$lat_grid <- round(predCoords$latitude / grid_size) * grid_size
predCoords$long_grid <- round(predCoords$longitude / grid_size) * grid_size

# Smoothing the covariates to estimate values for subsamples grid
grid_pred <- predCoords %>%
  group_by(lat_grid, long_grid) %>%
  summarise(avg_elevation = mean(Xelevation))

predZ <- scale(grid_pred)[,2:3]
colnames(predZ) <- c("Zlongitude","Zelevation")
pred_scaled <- as.data.frame(cbind(grid_pred[,1:2], predZ))

# Slightly different design matrix needed for nsgpPredict
PXmat <- unname(lm(rnorm(nrow(pred_scaled)) ~ Zelevation*Zlongitude, x = TRUE, 
                    data = pred_scaled)$x)

pred_A <- nsgpPredict(model = Rmodel_A,
                    samples = samples_A[seq(1,N),],
                    coords.predict = grid_pred[,1:2],
                    PX_Sigma = PXmat, PX_mu = PXmat)

write_csv(as.data.frame(pred_A$pred),
          "/Users/sanketdutta/Desktop/Spatial/Project/Data/Pred_A_2.csv")

pred_B <- nsgpPredict(model = Rmodel_B,
                      samples = samples_B[seq(1,N),],
                      coords.predict = grid_pred[,1:2],
                      PX_mu = PXmat)

write_csv(as.data.frame(pred_B$pred),
          "/Users/sanketdutta/Desktop/Spatial/Project/Data/Pred_B_2.csv")

#-----Posterior----

Means_A <- apply(pred_A$pred, 2, mean)
Means_B <- apply(pred_B$pred, 2, mean)

# Posterior Means in a single dataframe for facet based plot
Mean_Matrix <- rbind(cbind(as.matrix(grid_pred)[,1:2],
                           cbind(rep("NGP_SGV",length(Means_A)), Means_A)),
                     cbind(as.matrix(grid_pred)[,1:2],
                           cbind(rep("AGP_SGV", length(Means_B)), Means_B)))

colnames(Mean_Matrix) <- c("Lat","Lon","Experiment", "Mean")

SD_A <- apply(pred_A$pred, 2, sd)
SD_B <- apply(pred_B$pred, 2, sd)

# Posterior SDs in a single dataframe for facet based plot
SD_Matrix <- rbind(cbind(as.matrix(grid_pred)[,1:2],
                         cbind(rep("NGP_SGV",length(SD_A)), SD_A)),
                   cbind(as.matrix(grid_pred)[,1:2],
                         cbind(rep("AGP_SGV", length(SD_B)), SD_B)))

colnames(SD_Matrix) <- c("Lat","Lon","Experiment", "SD")

ggplot(as.data.frame(Mean_Matrix),
       aes(x = as.numeric(Lon), 
           y = as.numeric(Lat), 
           fill = as.numeric(Mean))) + 
  geom_tile() + facet_wrap(~Experiment) + 
  labs(x = NULL, y = NULL, fill = "") + 
  scale_fill_viridis_c(direction = -1) + coord_map()


ggsave(
  "Means_6.2.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

ggplot(as.data.frame(SD_Matrix),
       aes(x = as.numeric(Lon), 
           y = as.numeric(Lat), 
           fill = as.numeric(SD))) + 
  geom_tile() + facet_wrap(~Experiment) + 
  labs(x = NULL, y = NULL, fill = "") + 
  scale_fill_viridis_c(option = "A", direction = -1) + coord_map()

ggsave(
  "SDs_6.2.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

#----Eigen_plot----

Posterior_Dist <- as.data.frame(samples_A)

credible_interval <- function(column_data) {
  ci <- quantile(column_data, probs = c(0.025, 0.975))
  mean_value <- mean(column_data)
  c(ci, mean = mean_value)
}

ci_data <- sapply(Posterior_Dist, credible_interval)

# Selecting the design matrix coressponding parameters for eigenvalues
ci_Eigen <- as.data.frame(rbind(t(ci_data[,2:4]), t(ci_data[,6:8])))

ci_Eigen$variable <- c("Elevation", "Longitude", "Interaction",
                       "Elevation", "Longitude", "Interaction")
ci_Eigen$value <- c("Eigenvalue 1", "Eigenvalue 1", "Eigenvalue 1",
                    "Eigenvalue 2", "Eigenvalue 2", "Eigenvalue 2")

colnames(ci_Eigen) <- c("lower", "upper", "mean", "variable", "value")

ggplot(ci_Eigen, aes(x = variable)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0, color = "black") +
  geom_point(aes(y = mean), color = "black", size = 3) +
  labs(y = "Change per standardized unit") + 
  geom_hline(yintercept=0, color = "red") + facet_wrap(~value)

ggsave(
  "Eigen_CI.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

