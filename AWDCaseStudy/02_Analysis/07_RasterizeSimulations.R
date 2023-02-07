################################################################################
#### Rasterization of Simulated Dispersal Trajectories
################################################################################
# Description: In this script, we rasterize the simulated dispersal
# trajectories and create "heatmaps".

# Clear R's brain
rm(list = ls())

# Load required packages
library(terra)        # For quick raster manipulation
library(raster)       # For general raster manipulation
library(rgdal)        # To read spatial data
library(tidyverse)    # For data wrangling
library(rgeos)        # For manipulating vector data
library(lubridate)    # For working with dates
library(viridis)      # For nicer colors
library(pbmcapply)    # To show progress bar in mclapply calls
library(tmap)         # For nice spatial plots
library(davidoff)     # Custom functions
library(spatstat)     # For quick rasterization
library(maptools)     # For quick rasterization
library(reproj)       # For quick reprojection of coordinates

################################################################################
#### Load and Prepare Data
################################################################################
# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Load the simulated dispersal trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Ungroup them
sims <- ungroup(sims)

# Remove undesired columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "Area")]

# Reproject coordinates
sims[, c("x", "y")] <- reproj_xy(
    x      = sims[, c("x", "y")]
  , source = 4326
  , target = 32734
)

# Prepare extent that encompassess all coordinates + some buffer
ext <- extent(min(sims$x), max(sims$x), min(sims$y), max(sims$y)) +
  c(-1000, +1000, -1000, +1000)

# Span a raster with desired resolution
r <- raster(ext, res = 1000)
values(r) <- runif(ncell(r))
crs(r) <- CRS("+init=epsg:32734")
plot(r)

# Collect garbage
gc()

# Check out the number of rows
nrow(sims) / 1e6

################################################################################
#### Function to Rasterize Tracks
################################################################################
# Function to rasterize trajectories after desired number of steps and from
# desired source area
rasterizeSims <- function(
      simulations = NULL      # Simulated trajectories
    , raster      = NULL      # Raster onto which we rasterize
    , steps       = 500       # How many steps should be considered
    , area        = "Main"    # Simulations from which areas?
    , messages    = T         # Print update messages?
    , mc.cores    = detectCores() - 1
  ) {

  # Subset to corresponding data
  sub <- simulations[which(
      simulations$StepNumber <= steps
    & simulations$Area %in% area
  ), ]

  # Make sure raster values are all 0
  values(raster) <- 0

  # Create spatial lines
  sub_traj <- sims2tracks(
      simulations = sub
    , id          = "TrackID"
    , messages    = messages
    , mc.cores    = mc.cores
  )

  # Remove data by converting into spatial lines
  sub_traj <- as(sub_traj, "SpatialLines")
  crs(sub_traj) <- CRS("+init=epsg:32734")

  # Rasterize lines onto the cropped raster
  if (messages){
    cat("Rasterizing spatial lines...\n")
  }
  heatmap <- rasterizeSpatstat(
      l        = sub_traj
    , r        = raster
    , mc.cores = 1
  )

  # Store heatmap to temporary file
  heatmap <- writeRaster(heatmap, tempfile())
  crs(heatmap) <- CRS("+init=epsg:32734")

  # Return the resulting heatmap
  return(heatmap)
}

################################################################################
#### Rasterize Trajectories
################################################################################
# Create a dataframe with all source points and points in time at which we want
# to rasterize trajectories
rasterized <- as_tibble(
  expand.grid(
      steps            = c(125, 500, 2000)
    , area             = unique(sims$Area)
    , stringsAsFactors = F
  )
)

# Add a column for temporary but unique filename. Make sure the tempdir has
# plenty of storage.
rasterized$filename <- tempfile(
    pattern = paste0(
        "steps_", rasterized$steps
      , "_area_", rasterized$area
      , "_"
    )
  , fileext = ".tif"
)

# Loop through the study design and reasterize trajectories
heatmaps <- list()
for (i in 1:nrow(rasterized)) {

  # Create heatmap
  heatmaps[[i]] <- rasterizeSims(
      simulations = sims
    , raster      = r
    , steps       = rasterized$steps[i]
    , area        = rasterized$area[i]
    , messages    = T
    , mc.cores    = detectCores() - 1
  )

  # Clean garbage
  gc()

  # Print update
  cat(i, "/", nrow(rasterized), "done...\n")

}

# Combine maps
combined <- stack(heatmaps)

# Reproject them
combined <- rast(combined)
combined <- terra::project(combined, CRS("+init=epsg:4326"), method = "bilinear")
combined <- stack(combined)

# Crop them to the buffer
buffer <- readOGR("03_Data/03_Results/99_BufferArea.shp")
combined <- crop(combined, buffer)

# Store to file
writeRaster(combined, "03_Data/03_Results/99_Heatmaps.grd", overwrite = T)

# Add maps to the tibble
rasterized <- mutate(rasterized, heatmap = lapply(1:nlayers(combined), function(x){
  combined[[x]]
}))

# Store to file
write_rds(rasterized, "03_Data/03_Results/99_Heatmaps.rds")
rasterized <- read_rds("03_Data/03_Results/99_Heatmaps.rds")
