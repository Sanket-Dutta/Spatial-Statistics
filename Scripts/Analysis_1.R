#---------Libraries----

rm(list = ls())
library(BayesNSGP)
library(ggplot2)
library(readr)
library(dplyr)
library(maps)
library(ggforce)
nimbleOptions(verbose = FALSE) # Silences warnings

# Paper uses RW block samplers but package throws warnings

#-------------Setup----

set.seed(333) # For replication purposes

# Load raw data
COprecip1981 <- read.csv("COprecip1981.csv")
COpredDF <- read_csv("COpredDF.csv")
COprecip <- as.data.frame(COprecip1981)

# Various constants and distance matrices
coords <- as.matrix(COprecip[,c("Longitude", "Latitude")])
data <- COprecip$logPrecip
N <- nrow(coords)
Xmat <- unname(lm(logPrecip ~ Zelevation*Zslope10, x = TRUE, 
                  data = COprecip)$x) # Design Matrix

x_min <- min(coords[,1]); x_max <- max(coords[,1])
y_min <- min(coords[,2]); y_max <- max(coords[,2])
knot_coords <- expand.grid(
  lon = seq(from = x_min + 0.5*(x_max - x_min)/8,
            to = x_max - 0.5*(x_max - x_min)/8, length = 8),
  lat = seq(from = y_min + 0.5*(y_max - y_min)/8,
            to = y_max - 0.5*(y_max - y_min)/8, length = 8)
)

knot_groups <- matrix(c(1:4,9:12,17:20,25:28,
                        5:8,13:16,21:24,29:32,
                        32+c(1:4,9:12,17:20,25:28),
                        32+c(5:8,13:16,21:24,29:32)), ncol = 4)

N <- 2500 # Number of iterations (~ 2.5 hours for full code on 16GB M1 Macbook)

#-------Exploratory----

Colorado <- map_data("state")

Knot_A <- knot_coords[knot_groups[,1],]
Knot_B <- knot_coords[knot_groups[,2],]
Knot_C <- knot_coords[knot_groups[,3],]
Knot_D <- knot_coords[knot_groups[,4],]

ggplot() +
  geom_polygon(data = Colorado,
               aes(x=long, y=lat, group=group),
               fill = "white", colour = "black") +
  geom_polygon(aes(x=c(-109.5,-101.5,-101.5,-109.5), 
                   y=c(36.5,36.5,41.5,41.5)),
               fill = NA, colour = "red", stroke = 3) +
  coord_map() + xlab('') + ylab('')

ggsave(
  "Study_Area.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

ggplot() +
  geom_polygon(data = Colorado,
               aes(x=long, y=lat, group=group),
               fill = "white", colour = "black") +
  geom_point(data = Knot_A, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 0, size = 3, alpha = 0.75) +
  geom_point(data = Knot_B, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 1, size = 3, alpha = 0.75) +
  geom_point(data = Knot_C, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 2, size = 3, alpha = 0.75) +
  geom_point(data = Knot_D, color = "red", stroke = 1.5,
             aes(x=lon, y=lat), shape = 5, size = 3, alpha = 0.75) +
  geom_point(data = COprecip,
             aes(x=Longitude, y=Latitude,color=logPrecip), size = 3) +
  scale_colour_viridis_c() + xlab('') + ylab('') +
  coord_map(ylim = c(36.5,41.5), xlim = c(-109.5,-101.5))

ggsave(
  "CO_precip_1981.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

ggplot() +
  geom_polygon(data = Colorado,
               aes(x=long, y=lat, group=group),
               fill = "white", colour = "black") +
  geom_point(data = COprecip,
             aes(x=Longitude, y=Latitude, color=Elevation), size = 3) +
  scale_colour_gradientn(colours = terrain.colors(20))+ xlab('') + 
  ylab('') + coord_map(ylim = c(36.5,41.5), xlim = c(-109.5,-101.5))

ggsave(
  "Elevation_km.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

ggplot() +
  geom_polygon(data = Colorado,
               aes(x=long, y=lat, group=group),
               fill = "white", colour = "black") +
  geom_point(data = COprecip,
             aes(x=Longitude, y=Latitude,color=Slope10), size = 3) +
  scale_colour_gradient2(low = "red",
                         mid = "grey",
                         high = "cornflowerblue") + 
  xlab('') + ylab('') + coord_map(ylim = c(36.5,41.5), xlim = c(-109.5,-101.5))

ggsave(
  "Slope.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

#---------PacScherv----

constants_A <- list(nu = 2, mu_HP1 = 10,
                    Sigma_knot_coords = knot_coords,
                    Sigma_HP1 = c(10,10), Sigma_HP2 = rep(5,2),
                    Sigma_HP3 = rep(3.85,2), Sigma_HP4 = c(10,20),
                    maxAnisoDist = 16)

Rmodel_A <- nsgpModel(likelihood = "fullGP", coords = coords, 
                      data = data, constants = constants_A, 
                      Sigma_model = "npApproxGP")

conf_A <- configureMCMC(Rmodel_A)
conf_A$removeSamplers("w1_Sigma[1:64]")
conf_A$removeSamplers("w2_Sigma[1:64]")
conf_A$removeSamplers("w3_Sigma[1:64]")
for(h in 1:ncol(knot_groups)){
  conf_A$addSampler(target = c(paste0("w1_Sigma[",knot_groups[,h],"]")),
                    type = "RW_block")
  conf_A$addSampler(target = c(paste0("w2_Sigma[",knot_groups[,h],"]")), 
                    type = "RW_block")
  conf_A$addSampler(target = c(paste0("w3_Sigma[",knot_groups[,h],"]")),
                    type = "RW_block")
}

Rmcmc_A <- buildMCMC(conf_A) # Build the MCMC
Cmodel_A <- compileNimble(Rmodel_A) # Compile the model
Cmcmc_A <- compileNimble(Rmcmc_A, project = Rmodel_A) # Compile the MCMC

samples_A <- runMCMC(Cmcmc_A, niter = N)

write_csv(as.data.frame(samples_A),
          "/Users/sanketdutta/Desktop/Spatial/Project/Data/MCMC_A_1.csv")

#----------RissCald----

constants_B <- list(nu = 0.5, X_mu = Xmat, mu_HP1 = 10,
                    X_sigma = Xmat, log_sigma_HP1 = 10,
                    X_Sigma = Xmat, Sigma_HP1 = c(10,10),
                    Sigma_HP2 = c(2,2), maxAnisoDist = 16)

Rmodel_B <- nsgpModel(likelihood = "fullGP",
                      constants = constants_B,
                      coords = coords, data = data,
                      mu_model = "linReg", sigma_model = "logLinReg", 
                      Sigma_model = "covReg" )

conf_B <- configureMCMC(Rmodel_B)
conf_B$removeSamplers( c("psi11", "psi22", "rho") )
conf_B$addSampler( target = c("psi11", "psi22", "rho"), type = "RW_block" )

Rmcmc_B <- buildMCMC(conf_B) # Build the MCMC
Cmodel_B <- compileNimble(Rmodel_B) # Compile the model
Cmcmc_B <- compileNimble(Rmcmc_B, project = Rmodel_B) # Compile the MCMC

samples_B <- runMCMC(Cmcmc_B, niter = N)

write_csv(as.data.frame(samples_B),
          "/Users/sanketdutta/Desktop/Spatial/Project/Data/MCMC_B_1.csv")

#--------Predictive----

predCoords <- as.data.frame(COpredDF)

# Define the grid size based on RAM usage
grid_size <- 0.20

predCoords$lat_grid <- round(predCoords$latitude / grid_size) * grid_size
predCoords$long_grid <- round(predCoords$longitude / grid_size) * grid_size

# Smoothing the covariates to estimate values for subsamples grid
grid_pred <- predCoords %>%
  group_by(lat_grid, long_grid) %>%
  summarise(avg_elevation = mean(elevation),
            avg_slope = mean(slope))

predZ <- scale(grid_pred)[,3:4]
colnames(predZ) <- c("Zelevation", "Zslope10")
pred_scaled <- as.data.frame(cbind(grid_pred[,1:2], predZ))

# Slightly different design matrix needed for nsgpPredict
Xmat_pred <- unname(lm(rnorm(nrow(pred_scaled)) ~ Zelevation*Zslope10, x = TRUE, 
                       data = pred_scaled)$x)

pred_constants <- list( PX_sigma = Xmat_pred,
                        PX_Sigma = Xmat_pred,
                        PX_mu = Xmat_pred)

pred_A <- nsgpPredict(model = Rmodel_A, 
                      coords.predict = grid_pred[,1:2],
                      samples = samples_A[seq(1,N),])

write_csv(as.data.frame(pred_A$pred),
          "/Users/sanketdutta/Desktop/Spatial/Project/Data/Pred_A_1.csv")

pred_B <- nsgpPredict(model = Rmodel_B, 
                      coords.predict = grid_pred[,1:2],
                      samples = samples_B[seq(1,N),],
                      constants = pred_constants)

write_csv(as.data.frame(pred_B$pred),
          "/Users/sanketdutta/Desktop/Spatial/Project/Data/Pred_B_1.csv")

#---------Posterior----

Means_A <- apply(pred_A$pred, 2, mean)
Means_B <- apply(pred_B$pred, 2, mean)

# Posterior Means in a single dataframe for facet based plot
Mean_Matrix <- rbind(cbind(as.matrix(grid_pred)[,1:2],
                           cbind(rep("PacScherv",length(Means_A)), Means_A)),
                     cbind(as.matrix(grid_pred)[,1:2],
                           cbind(rep("RissCald", length(Means_B)), Means_B)))

colnames(Mean_Matrix) <- c("Lat","Lon","Experiment", "Mean")

SD_A <- apply(pred_A$pred, 2, sd)
SD_B <- apply(pred_B$pred, 2, sd)

# Posterior SDs in a single dataframe for facet based plot
SD_Matrix <- rbind(cbind(as.matrix(grid_pred)[,1:2],
                         cbind(rep("PacScherv",length(SD_A)), SD_A)),
                   cbind(as.matrix(grid_pred)[,1:2],
                         cbind(rep("RissCald", length(SD_B)), SD_B)))

colnames(SD_Matrix) <- c("Lat","Lon","Experiment", "SD")

ggplot(as.data.frame(Mean_Matrix),
       aes(x = as.numeric(Lon), 
           y = as.numeric(Lat), 
           fill = as.numeric(Mean))) + 
  geom_tile() + facet_wrap(~Experiment) + 
  labs(x = NULL, y = NULL, fill = "") + 
  scale_fill_viridis_c() + coord_map()

ggsave(
  "Means_6.1.pdf",
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
  scale_fill_viridis_c(option = "A") + coord_map()

ggsave(
  "SDs_6.1.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)
  
#-------RissCald_CI----

Posterior_Dist <- as.data.frame(samples_B)

# Calculating Mean + 95% CI

credible_interval <- function(column_data) {
  ci <- quantile(column_data, probs = c(0.025, 0.975))
  mean_value <- mean(column_data)
  c(ci, mean = mean_value)
}

ci_data <- sapply(Posterior_Dist, credible_interval)
ci_df <- data.frame(t(ci_data))
ci_df$variable <- rownames(ci_df)

ggplot(ci_df, aes(x = variable)) +
  geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                width = 0, color = "darkgreen") +
  geom_point(aes(y = mean), color = "darkgreen", size = 3) +
  labs(y = "Posterior Mean + 95% CI", x = "Parameter")


ggsave(
  "RissCald_CI.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)

#----PacSherv_Sigma----

Posterior_Sigma <- as.matrix(colMeans(samples_A[,10:201]))

Posterior_Sigma <- cbind(cbind(Posterior_Sigma[1:64,],
                               Posterior_Sigma[65:128,]),
                         Posterior_Sigma[129:192,])

rownames(Posterior_Sigma) <- seq(1,64)

# Outputs of the model
colnames(Posterior_Sigma) <- c("log_Eigen_1", "log_Eigen_2",
                               "Logit_Scaled_Rotation")

Sigma_Matrix <- matrix(0, nrow = 64, ncol = 8)

# Parameters needed to recreate Sigma matrix
colnames(Sigma_Matrix) <- c("Eigen_1", "Eigen_2", "Eigen_diff", "Cos", "Sin", 
                            "Cos^2", "Sin^2", "Sin*Cos")

Sigma_Matrix[1:64,1] <- exp(Posterior_Sigma[1:64,1])
Sigma_Matrix[1:64,2] <- exp(Posterior_Sigma[1:64,2])
Sigma_Matrix[1:64,3] <- Sigma_Matrix[1:64,1] - Sigma_Matrix[1:64,2]
Sigma_Matrix[1:64,4] <- cos(pi*expit(Posterior_Sigma[1:64,3])/2)
Sigma_Matrix[1:64,5] <- sin(pi*expit(Posterior_Sigma[1:64,3])/2)
Sigma_Matrix[1:64,6] <- Sigma_Matrix[1:64,4]**2
Sigma_Matrix[1:64,7] <- Sigma_Matrix[1:64,5]**2
Sigma_Matrix[1:64,8] <- Sigma_Matrix[1:64,4]*Sigma_Matrix[1:64,5]

c <- sqrt(qchisq(0.05, df = 2))

Ellipse <- matrix(0, nrow = 64, ncol = 5)

# Parameters needed to plot geom_ellipse
colnames(Ellipse) <- c("x0", "y0", "a", "b", "angle")

Ellipse[1:64,1] <- knot_coords$lon
Ellipse[1:64,2] <- knot_coords$lat
Ellipse[1:64,3] <- c*sqrt((Sigma_Matrix[1:64,1]*Sigma_Matrix[1:64,6]) +
                          (Sigma_Matrix[1:64,2]*Sigma_Matrix[1:64,7]) +
                          (Sigma_Matrix[1:64,3]*Sigma_Matrix[1:64,8]))
Ellipse[1:64,4] <- c*sqrt((Sigma_Matrix[1:64,1]*Sigma_Matrix[1:64,7]) +
                          (Sigma_Matrix[1:64,2]*Sigma_Matrix[1:64,6]) +
                          (Sigma_Matrix[1:64,3]*Sigma_Matrix[1:64,8]))
Ellipse[1:64,5] <- pi*expit(Posterior_Sigma[1:64,3])/2

Ellipse[is.nan(Ellipse)] <- 0

Ellipse <- as.data.frame(Ellipse)

ggplot() +
  geom_polygon(data = Colorado,
               aes(x=long, y=lat, group=group),
               fill = "white", colour = "black") +
  geom_ellipse(data = Ellipse, color = "red",
               aes(x0 = x0, y0 = y0, a = a, b = b, angle = angle)) +
  geom_point(data = knot_coords, color = "black", aes(x=lon, y=lat), 
             shape = 3, size = 2, stroke = 0.5) + xlab('') + ylab('') +
  coord_map(ylim = c(36.5,41.5), xlim = c(-109.5,-101.5))


ggsave(
  "Ellipse.pdf",
  device = "pdf",
  path = "/Users/sanketdutta/Desktop/Spatial/Project/Plots/",
  dpi = "retina"
)
