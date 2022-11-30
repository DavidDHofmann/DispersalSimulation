################################################################################
#### Simulating Movement with Known Preferences
################################################################################
# Description: In this script, we simulate movement data with known preferences.
# However, this is not yet the data we will use to conduct connectivity
# analyses. Rather, it is supposed to generate the kind of data you would
# collect in the field and then use to fit a step selection model. To simulate
# the data, we will use the exact same framework (i.e. step selection functions)
# that we will later use to simulate virtual dispersers to assess connectivity.

# Clear R's brain
rm(list = ls())

# Load required packages
library(RandomFields)   # To simulate covariates
library(raster)         # To handle spatial data
library(foreach)        # For parallel computing
library(doSNOW)         # For parallel computing
library(parallel)       # For parallel computing
library(tidyverse)      # For data wrangling
library(sf)             # For nice plots
library(lubridate)      # To handle times
library(rgdal)          # To save shapefile

# Set working directory
# setwd("/home/david/ownCloud/DispersalSimulation")

# Load custom functions
source("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/DispersalSimulation/00_Functions.R")

################################################################################
#### Simulate Covariates
################################################################################
# Set seed for reproducability
set.seed(123)

# Create a polygon for the extent of the study area
ext <- extent(c(0, 100, 0, 100))
ext <- as(ext, "SpatialPolygons")

# Simulate an elevation layer
elev <- RMexp(var = 5, scale = 30) +
  RMnugget(var = 0.5) +
  RMtrend(mean = 0)
elev <- RFsimulate(elev, x = 1:500, y = 1:500)
elev <- raster(elev)
elev <- setExtent(elev, ext)

# Let's also create a gradient that we can add to the elevation layer
x <- seq(-1, 1, len = ncol(elev))
y <- seq(-1, 1, len = nrow(elev))
gradient <- raster(elev)
values(gradient) <- outer(x, y, function(x, y){0 * x + 1 * y})

# Now add the gradient to our elevation raster
elev <- elev + 5 * gradient

# Assume there are three national parks inside the study area
np1 <- Polygon(data.frame(
    x <- c(27.16, 27.13, 25.74, 27.04, 30.63, 33.42, 34.61, 31.88, 29.36, 28.64, 27.16)
  , y <- c(67.75, 69.35, 71.19, 74.07, 75.37, 73.92, 71.28, 70.96, 69.95, 67.66, 67.75)
))
np2 <- Polygon(data.frame(
    x <- c(33.30, 35.53, 34.81, 30.75, 27.19, 28.94, 32.23, 33.92, 32.26, 33.30)
  , y <- c(27.56, 27.59, 22.79, 23.11, 25.90, 30.32, 33.49, 31.77, 29.43, 27.56)
))
np3 <- Polygon(data.frame(
    x <- c(75.21, 78.33, 79.71, 75.97, 73.57, 75.03, 68.72, 66.45, 66.67, 69.74, 71.96, 75.21)
  , y <- c(49.60, 49.87, 47.15, 43.73, 43.06, 38.48, 38.88, 43.37, 45.68, 46.75, 49.07, 49.60)
))
np1 <- Polygons(list(np1), 1)
np2 <- Polygons(list(np2), 2)
np3 <- Polygons(list(np3), 3)
nps <- SpatialPolygons(list(np1, np2, np3))
nps$ParkName <- c("Northern NP", "Southern NP", "Eastern NP")

# Also put a spatial point inside each national park. We'll later use them as
# points of attraction.
pts <- SpatialPoints(rbind(
    c(30, 73)
  , c(30, 28)
  , c(72, 45)
))

# We will later assume that animals prefer to stend close to these national
# parks. Hence, let's generate a layer rendering the (square-rooted) distance to
# each national park.
dist <- distanceFromPoints(elev, pts)
dist <- sqrt(dist)

# Standardize (center and scale) each of the covariates
elev <- (elev - cellStats(elev, mean)) / cellStats(elev, sd)
dist <- (dist - cellStats(dist, mean)) / cellStats(dist, sd)

# Put covariate layers into a stack
cov <- stack(elev, dist)
names(cov) <- c("elev", "dist")

# We slightly expand each covariate layer and randomize covariate values within
# the buffer zone.
cov <- extendRaster(cov, extent(-20, 120, -20, 120))

# We also want to specify the extent on which our animals are allowed to move.
# For now we will limit movement to the main study area, but we won't allow any
# movement through the buffer zone. Later, when we simulate movement to assess
# connectivity, we will use the buffer zone to mitigate edge effects. Here,
# however, this is not required.
ext_move <- extent(c(0, 100, 0, 100))
ext_move <- as(ext_move, "SpatialPolygons")

# Visualize all covariates + extent
cov %>%
  as.data.frame(xy = T) %>%
  gather(key = covariate, value = value, 3:ncol(.)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = value)) +
    geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
    geom_sf_text(data = st_as_sf(nps), aes(label = ParkName), col = "white"
      , nudge_y = 10, size = 3) +
    geom_sf(data = st_as_sf(ext_move), col = "red", fill = NA, lty = 2) +
    facet_wrap("covariate") +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    coord_sf()

################################################################################
#### Simulating a Single Trajectory
################################################################################
# Let's start simple and simulate a single trajectory. For this, we first need
# to specify some simulation parameters.
# formula   : model formula that is used to span the design matrix
# prefs     : preferences for each term in the resulting model matrix
# sl_dist   : step length distribution used to draw random steps (same format as in the amt package)
# ta_dist   : turning angle distribution used to draw random steps (same format as in the amt package)
# n_rsteps  : number of random steps proposed at each iteration
# n_steps   : number of steps simulated in total (realized)
# stop      : true/false whether the simulation should stop in case an individual hits a map boundary

# Simulation Parameters
formula <- ~ elev + dist + cos_ta  # They need to match the covariates
prefs     <- c(0.5, -0.2, 1)       # They need to match the formula!!!
sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 0.5))
ta_dist   <- list(name = "vonmises", params = list(kappa = 0, mu = 0))
n_rsteps  <- 25
n_steps   <- 200
stop      <- T

# Sample a random starting location within any national park
pts <- coordinates(spsample(nps, type = "random", n = 1))

# Apply the "move" function an let a virtual disperser move across the
# landscape. The move function simulates a trajectory using the step selection
# framework.
sim <- move(
    xy       = pts
  , covars   = cov
  , formula  = formula
  , prefs    = prefs
  , sl_dist  = sl_dist
  , ta_dist  = ta_dist
  , ext      = ext
  , n_steps  = n_steps
  , n_rsteps = n_rsteps
  , stop     = stop
)

# And plot the trajectory
cov[[1]] %>%
  as.data.frame(xy = T) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = elev)) +
    geom_point(data = sim, inherit.aes = F, aes(x = x, y = y), size = 0.8) +
    geom_path(data = sim, inherit.aes = F, aes(x = x, y = y), size = 0.3) +
    geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
    geom_sf_text(data = st_as_sf(nps), aes(label = ParkName), col = "white"
      , nudge_y = 10, size = 2) +
    geom_sf(data = st_as_sf(ext_move), col = "red", fill = NA, lty = 2) +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    coord_sf()

# We can even plot a nice animation of the trajectory (we need to plot it to an
# external windo though)
range_x <- range(sim$x)
range_y <- range(sim$y)
x11()
for (i in 1:nrow(sim)){
  plot(sim$y[1:i] ~ sim$x[1:i]
    , type = "o"
    , pch  = 20
    , xlim = range_x
    , ylim = range_y
    , xlab = "x"
    , ylab = "y"
    , main = paste0("A First Simulation - Step ", i)
  )
  Sys.sleep(0.05)
}

################################################################################
#### Multiple Trajectories
################################################################################
# Now we expand the code from above to simulate multiple individuals (let's say
# 10). To make bookkeeping easier, we prepare a tibble.
sims <- tibble(ID = 1:25)

# Also sample sourcepoint for each individual within national park. For
# simplicity, we will place source points randomly. However, you could also
# initiate the same number of individuals within each national park.
pts <- coordinates(spsample(nps, type = "random", n = nrow(sims)))

# Also sample random start points for each individual
plot(cov[["elev"]], asp = 1.03)
points(pts, pch = 20, cex = 0.2)

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

# Let's take a look at the final object
print(sims)

# You can see that not all trajectories consist of 200 steps. This is because
# some individuals hit a map boundary and we specified to stop the simulation in
# such cases.
summary(sapply(sims$simulations, nrow))

################################################################################
#### Cleaning and Visualizing Simulations
################################################################################
# Unnest the simulations
sims <- sims %>% unnest(simulations)

# Make sure ID's are treated as factors
sims$ID <- as.factor(sims$ID)

# Assign a unique ID to each step
sims$step_id <- 1:nrow(sims)

# Also assign a timestamp to each step (with som randomness)
sims$timestamp <- ymd_hms("2021-01-01 11:00:00") +
  sims$step_number * hours(1) +
  seconds(rnorm(n = nrow(sims), mean = 0, sd = 300))

# Take a look at the simulated data
head(sims, 20)

# Visualize simulated tracks
cov[[1]] %>%
  as.data.frame(xy = T) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = elev)) +
    geom_point(data = sims, inherit.aes = F, aes(x = x, y = y, col = ID), size = 0.8) +
    geom_path(data = sims, inherit.aes = F, aes(x = x, y = y, col = ID), size = 0.3) +
    geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
    geom_sf_text(data = st_as_sf(nps), aes(label = ParkName), col = "white"
      , nudge_y = 10, size = 3) +
    geom_sf(data = st_as_sf(ext_move), col = "red", fill = NA, lty = 2) +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    coord_sf() +
    theme(legend.position = "none")

# In reality, we don't observe step lengths, turning angles etc. but have to
# derive them from xy coordinates. Hence, let's assume that we only observed xy
# data + ID and remove the rest.
obs <- dplyr::select(sims, x, y, ID, step_number, step_id, timestamp)

# Also, in reality we might have some missing fixes. Let's randomly remove some
# fixes (here, 5%).
remove <- sample(1:nrow(obs), size = nrow(obs) * 0.05, replace = F)
obs <- obs[-remove, ]

# Let's look at the final dataframe
head(obs, 20)

# This is the type of data that we would observe in reality. Let's store it to
# file so that we can later analyse it. Also store the simulated spatial layers
# to file, as well as the national parks.
write_rds(obs, "99_ObservedMovements.rds")
write_rds(cov, "99_CovariateLayers.rds")
write_rds(nps, "99_NationalParks.rds")

# Let's also store the true preferences used for the simulation. Our goal later
# will be to estimate these using a step selection model.
true_prefs <- data.frame(
    Coefficient      = prefs
  , Covariate        = c("elev", "dist", "cos_ta")
  , stringsAsFactors = F
)
true_prefs
write_rds(true_prefs, "99_TruePreferences.rds")
