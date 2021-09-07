################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Use ISSFs to simulate movement trajectories with known
# preferences

# Clear R's brain
rm(list = ls())

# Load required packages
library(RandomFields)   # To simulate covariates
library(raster)         # To handle spatial data
library(pbmcapply)      # For parallel computing (on linux/mac)
library(foreach)        # For parallel computing (on windows)
library(doSNOW)         # For parallel computing (on windows)
library(tidyverse)      # For data wrangling
library(sf)             # For nice plots
library(lubridate)      # To handle times
library(rgdal)          # To save shapefile

# Set working directory
setwd("/home/david/ownCloud/DispersalSimulation")

# Load custom functions
source("00_Functions.R")

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
elev <- as.matrix(elev)
elev <- raster(elev)
elev <- setExtent(elev, ext)

# Let's also create a gradient that we can add to the elevation layer
x <- seq(-1, 1, len = ncol(elev))
y <- seq(-1, 1, len = nrow(elev))
gradient <- raster(elev)
values(gradient) <- outer(x, y, function(x, y){0 * x + 1 * y})

# Now add the gradient to our elevation raster
elev <- elev + 3 * gradient

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

# Also put three points of attraction inside each national park
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

# Normalize each of the covariates
elev <- (elev - cellStats(elev, mean)) / cellStats(elev, sd)
dist <- (dist - cellStats(dist, mean)) / cellStats(dist, sd)

# Put covariate layers into a stack
cov <- stack(elev, dist)
names(cov) <- c("elev", "dist")

# We expand each covariate layer to this extent and randomize covariate values
# within the buffer zone
cov <- extendRaster(cov, extent(-20, 120, -20, 120))

# In any case, for now we will limit movement to the main study area
ext_move <- extent(c(0, 100, 0, 100))
ext_move <- as(ext_move, "SpatialPolygons")

# Visualize all covariates + extent
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
#### Simulating a Single Trajectory
################################################################################
# Simulation Parameters
formula <- ~ elev + dist + cos_ta  # They need to match the covariates
prefs     <- c(0.5, -0.2, 1)       # They need to match the formula!!!
sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 0.5))
ta_dist   <- list(name = "vonmises", params = list(kappa = 0, mu = 0))
n_rsteps  <- 25
n_steps   <- 200
stop      <- T

# To try out the function, sample a random starting location within a national
# park
pts <- coordinates(spsample(nps, type = "random", n = 1))

# Apply the simulation function an let a virtual disperser move across the
# landscape
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

# And plot it
cov[[1]] %>%
  as.data.frame(xy = T) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = elev)) +
    geom_point(data = sim, inherit.aes = F, aes(x = x, y = y), size = 0.8) +
    geom_path(data = sim, inherit.aes = F, aes(x = x, y = y), size = 0.3) +
    geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
    geom_sf(data = st_as_sf(ext_move), col = "red", fill = NA, lty = 2) +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    coord_sf()

# We can even plot a nice animation of the trajectory
range_x <- range(sim$x)
range_y <- range(sim$y)
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
}

################################################################################
#### Multiple Trajectories (in parallel on mac/linux)
################################################################################
# Now we expand the code from above to simulate multiple individuals. To make
# bookkeeping easier, we prepare a tibble.
sims <- tibble(ID = 1:100)

# Also sample sourcepoint for each individual within national park. For
# simplicity, we will place source points randomly. However, you could also
# initiate the same number of individuals within each national park.
pts <- coordinates(spsample(nps, type = "random", n = nrow(sims)))

# Visualize
plot(cov[["elev"]])
points(pts, pch = 20, cex = 0.2)

# Simulate movement for each individual
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
sims <- tibble(ID = 1:100)

# Also sample sourcepoint for each individual within national park. For
# simplicity, we will place source points randomly. However, you could also
# initiate the same number of individuals within each national park.
pts <- coordinates(spsample(nps, type = "random", n = nrow(sims)))

# Also sample random start points for each individual
plot(cov[["elev"]])
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
sims$timestamp <- as.POSIXct("2021-01-01 11:00:00") +
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
    geom_sf(data = st_as_sf(ext_move), col = "red", fill = NA, lty = 2) +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    coord_sf() +
    theme(legend.position = "none")

# In reality, we don't observe step lengths, turning angles etc. and we have to
# derive them from xy coordinates. Hence, let's assume that we only observed xy
# data + ID. Remove the rest.
obs <- dplyr::select(sims, x, y, ID, step_number, step_id, timestamp)

# Also, in reality we might have some missing fixes. So let's randomly remove
# some fixes (here, 5%).
remove <- sample(1:nrow(obs), size = nrow(obs) * 0.05, replace = F)
obs <- obs[-remove, ]

# This is the type of data that we would observe in reality. Let's store it to
# file. Also store the simulated spatial layers to file
write_csv(obs, "ObservedMovements.csv")
writeRaster(cov, "CovariateLayers.grd", overwrite = T)
writeOGR(nps, ".", "NationalParks", driver = "ESRI Shapefile", overwrite = T)

# Let's also store the true preferences
true_prefs <- data.frame(
    Coefficient      = prefs
  , Covariate        = c("elev", "dist", "cos_ta")
  , stringsAsFactors = F
)
write_csv(true_prefs, "TruePreferences.csv")
