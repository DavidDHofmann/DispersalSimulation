################################################################################
#### Dispersal Simulation
################################################################################
# Description: Based on derived selection coefficients, we now simulate
# dispersers moving across the study area.

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(viridis)        # For nice colors
library(rgdal)          # To load shapefiles
library(rgeos)          # For geometry manipulation
library(pbmcapply)      # For parallel computing (on linux/mac)
library(foreach)        # For parallel computing (on windows)
library(doSNOW)         # For parallel computing (on windows)
library(sf)             # To plot spatial things with ggplot

# Set working directory
setwd("/home/david/ownCloud/DispersalSimulation")

# Load custom functions
source("00_Functions.R")

# Load covariate layers and estimated coefficients
cov <- stack("CovariateLayers.grd")
nps <- readOGR("NationalParks.shp")
beta <- read_csv("Estimates.csv")

# Make sure covariate layers are loaded into memory for maximal efficiency
cov <- readAll(cov)
inMemory(cov)

# Create a polygon for the extent of the study area
ext <- extent(cov)
ext <- as(ext, "SpatialPolygons")

# Define area on which simulated individuals are allowed to move. Note that we
# allow virtual dispersers to move through the buffer zone.
ext_move <- extent(-20, 120, -20, 120)
ext_move <- as(ext_move, "SpatialPolygons")

# Visualize covariates
cov %>%
  as.data.frame(xy = T) %>%
  gather(key = covariate, value = value, 3:ncol(.)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = value)) +
    geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
    geom_sf(data = st_as_sf(ext_move), col = "red", fill = NA, lty = 2) +
    facet_wrap("covariate") +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    coord_sf()

################################################################################
#### Simulate Dispersers
################################################################################
# Based on the derived selection coefficients, we can now simulate dispersers
formula <- ~ elev + dist + cos_ta
prefs     <- beta$Estimate
sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 0.5))
ta_dist   <- list(name = "vonmises", params = list(kappa = 0, mu = 0))
n_rsteps  <- 25
n_steps   <- 200
stop      <- F

# Let's also sample a source point. This time, we will initiate 500 individuals
# within each national park
pts <- lapply(1:nrow(nps), function(x) {
  spsample(nps[x, ], type = "random", n = 500)
})
pts <- do.call(rbind, pts)
pts <- coordinates(pts)

# Visualize
plot(cov[["elev"]])
points(pts, pch = 20, cex = 0.2)

################################################################################
#### Multiple Trajectories (in parallel on mac/linux)
################################################################################
# Now we expand the code from above to simulate multiple individuals. To make
# bookkeeping easier, we prepare a tibble.
sims <- tibble(ID = 1:nrow(pts))

# Simulate movement for each individual (this will take pretty long)
sims$simulations <- pbmclapply(1:nrow(sims)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x){
    move(
        xy       = rbind(pts[x, ])
      , covars   = cov
      , formula  = formula
      , prefs    = prefs
      , sl_dist  = sl_dist
      , ta_dist  = ta_dist
      , ext      = ext_move
      , n_steps  = n_steps
      , n_rsteps = n_rsteps
      , stop     = stop
    )
})

# Let's take a look at the final object
print(sims)

################################################################################
#### Multiple Trajectories (in parallel on windows)
################################################################################
# Now we expand the code from above to simulate multiple individuals (let's say
# 10). To make bookkeeping easier, we prepare a tibble.
sims <- tibble(ID = 1:nrow(pts)

# Prepare cluster
cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

# Prepare progressbar
pb <- txtProgressBar(max = nrow(sims), style = 3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress = progress)

# Run the simulation
sims$simulations <- foreach(x = 1:nrow(sims), .packages = "raster", .options.snow = opts) %dopar% {
  move(
      xy       = rbind(pts[x, ])
    , covars   = cov
    , formula  = formula
    , prefs    = prefs
    , sl_dist  = sl_dist
    , ta_dist  = ta_dist
    , ext      = ext_move
    , n_steps  = n_steps
    , n_rsteps = n_rsteps
    , stop     = stop
  )
}

# Close cluster
registerDoSEQ()
stopCluster(cl)

################################################################################
#### Cleaning and Visualizing Simulations
################################################################################
# Unnest the simulations
sims <- sims %>% unnest(simulations)

# Visualize them
ids <- unique(sims$ID)
plot(cov[["elev"]])
for (i in 1:length(unique(sims$ID))){
  sub <- subset(sims, ID == ids[i])
  points(sub$y ~ sub$x, col = i, pch = 16, cex = 0.4)
  lines(sub$y ~ sub$x, col = i, lwd = 0.3)
}

# Store simulations
write_csv(sims, "SimulatedMovements.csv")
