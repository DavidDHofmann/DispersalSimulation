################################################################################
#### Betweenness Analysis
################################################################################
# Description: In this script we will analyse the simulated trajectories using a
# network approach. Ultimately, we will generate a map illustrating the
# betweenness of each pixel in our study area.

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(parallel)       # To run on multiple cores
library(pbmcapply)      # To run on multiple cores with progress bar
library(igraph)         # For network analysis
library(viridis)        # For nice colors
library(RColorBrewer)   # For more nice colors
library(davidoff)       # Custom functions
library(spatstat)       # For the interpolation

# Set a seed
set.seed(12345)

# Suppress scientific notation
options(scipen = 999)

################################################################################
#### Functions to Calculate Network Metrics
################################################################################
# We will create and use networks in order to derive network metrics that we can
# plot spatially. Let's prepare a function that we can use to derive all desired
# network metrics and depict them on a raster.
netMet <- function(
    network   = NULL    # The network based on which the metrics are calculated
  , raster    = NULL    # The raster onto which the metrics are calculated
  , tempfile  = F       # Should the resulting raster go to a temporary file?
  , metrics = c("betweenness", "closeness", "degree")
  ) {

    # Calculate the desired network metrics and return them as rasters
    result <- vector(mode = "list", length = 3)
    if ("betweenness" %in% metrics){
      betweenness <- raster
      values(betweenness) <- betweenness(network)
      names(betweenness) <- "betweenness"
      result[[1]] <- betweenness
    }
    if ("closeness" %in% metrics){
      closeness <- raster
      values(closeness) <- closeness(network)
      names(closeness) <- "closeness"
      result[[2]] <- closeness
    }
    if ("degree" %in% metrics){
      degree <- raster
      values(degree) <- degree(network)
      names(degree) <- "degree"
      result[[3]] <- degree
    }

    # Remove NULLs from the list
    result <- plyr::compact(result)

    # Put all into a stack
    result <- stack(result)

    # In case the user desires to store a temporary file, do so
    if (tempfile){
      result <- writeRaster(result, tempfile())
    }

    # Return the final stack
    return(result)
}

# Function to retrieve the visitation history from a sequence of values
visitHist <- function(x, singlecount = F) {
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections = 1
  }
  return(transitions)
}

# Function to interpolate a path
interpolatePath <- function(x, y, eps = 0.1) {
  inter <- lapply(1:(length(x) - 1), function(i) {
    xy_new <- interpolatePointsC(
        x1 = x[i]
      , x2 = x[i + 1]
      , y1 = y[i]
      , y2 = y[i + 1]
      , by = eps
    )
    return(xy_new)
  }) %>% do.call(rbind, .)
  inter <- as.data.frame(inter)
  names(inter) <- c("x", "y")
  return(inter)
}

################################################################################
#### Load and Clean Data
################################################################################
# Load the simulated dispersal trajectories
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")
sims <- ungroup(sims)

# Prepare design through which we want to loop (currently, we will simply use
# the "number of steps" and check how betweenness changes depending on the
# number of simulated steps)
design <- tibble(steps = c(125, 500, 2000))

# Load the reference raster. We'll use it to generate a grid that will
# ultimately serve as "network". The centerpoint of each grid-cell is a network
# node and is spatially referenced!
r   <- raster("03_Data/02_CleanData/00_General_Raster.tif")
r   <- aggregate(r, fact = 2500 / 250, fun = max)
ver <- 1:ncell(r)
lay <- as.matrix(as.data.frame(r, xy = T)[, c(1, 2)])
r[] <- ver

# Loop through the design and calculate a betweenness map for each treatment
maps <- list()
for (i in 1:nrow(design)) {

  # Subset data to desired steps
  sub <- subset(sims, StepNumber <= design$steps[i], )

  # Nest tracks by their IDs
  sub <- nest(sub, data = -TrackID)

  # Determine the visitation history of each path. Note that we will
  # "interpolate" each of the simulated steps. This will allow us to determine
  # cell-transitions at a much finer scale than if we would simply use the start
  # and endpoint of each step.
  cat("Getting visitation history...\n")
  history <- pbmclapply(
      X                  = sub$data
    , ignore.interactive = T
    , mc.cores           = 1
    , FUN                = function(path) {
      path   <- interpolatePath(path$x, path$y, eps = 500 / 111000)
      visits <- extract(r, path)
      visits <- visitHist(visits, singlecount = T)
      return(visits)
    }) %>%
    do.call(rbind, .) %>%
    group_by(from, to) %>%
    summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
    ungroup() %>%
    mutate(weight = mean(TotalConnections) / TotalConnections)

  # Create network
  cat("Creating graph...\n")
  net <- graph_from_data_frame(history, vertices = ver)

  # Calculate Betweenness
  cat("Calculating betweenness...\n")
  met <- netMet(
      network  = net
    , raster   = r
    , metrics  = "betweenness"
    , tempfile = T
  )

  # Put into list
  maps[[i]] <- met

  # Update
  cat(i, "out of", nrow(design), "done\n")

}

# Store them
maps <- stack(maps)
plot(sqrt(maps), col = magma(100))

# Store maps to file
writeRaster(maps, "03_Data/03_Results/99_Betweenness.grd")
