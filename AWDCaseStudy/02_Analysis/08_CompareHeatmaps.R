################################################################################
#### Comparison of Heatmap of Simulated Trajectories and Least Cost Corridors
################################################################################
# Description: In this script, we are going to compare the rasterized simulate
# trajectories (heatmaps) to the least-cost corridors that we predicted in an
# earlier paper. Basically, we want to know how similar the heatmaps are to the
# permeability and least-cost maps.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)     # For general raster manipulation
library(tidyverse)  # For data wrangling
library(davidoff)   # Custom functions

################################################################################
#### Prepare Data
################################################################################
# Reload rasterized maps. We have to load both, the tibble (containing some
# useful information), as well as the rasters themselves. The reason the two
# files are seperate is that in the tibble the heatmaps are not linked correctly
# anymore.
rasterized <- read_rds("03_Data/03_Results/99_RasterizedSimulationsBootstrap.rds")
heatmaps <- stack("03_Data/03_Results/99_RasterizedSimulationsBootstrap.tif")

# Remove unneccessary columns
rasterized$filename <- NULL
rasterized$heatmap2 <- NULL

# Load least-cost corridors and the permeability map that we want to compare our
# heatmaps to.
perm <- raster("03_Data/03_Results/99_PermeabilityMap.tif")
corr <- raster("03_Data/03_Results/99_LeastCostCorridors.tif")

# Aggregate corridor- and permeability-map to same resolution as heatmaps
perm <- aggregate(perm, fact = 10, fun = "mean")
corr <- aggregate(corr, fact = 10, fun = "mean")

# In the heatmaps there is one row missing. Let's extend them to match the
# permeability and corridor maps. Note that we'll replace the missing row simply
# with the last row containing values. Not the most beautiful solution but it
# will do the job.
heatmaps <- extend(heatmaps, corr, value = 0)
heatmaps[nrow(heatmaps), ] <- heatmaps[nrow(heatmaps)- 1, ]

# Some metrics require the maps to be on the same scale. Let's therefore
# normalize the maps.
perm <- normalizeMap(perm)
corr <- normalizeMap(corr)
for (i in 1:nlayers(heatmaps)){
  heatmaps[[i]] <- normalizeMap(heatmaps[[i]])
  heatmaps[[i]] <- writeRaster(heatmaps[[i]], tempfile())
}

################################################################################
#### Compare Heatmaps to Permeability and Corridor Maps
################################################################################
# Function to calculate Bhattacharyyas affinity index
bhattacharyya_affinity <- function(x, y){

  # Make sure map values add up to one
  map1 <- x / sum(values(x))
  map2 <- y / sum(values(y))

  # Calculate affinity
  affinity <- sum(values(sqrt(map1) * sqrt(map2)))

  # Return the result
  return(affinity)
}

# Let's loop through all maps and identify the correlation and bhattacharyyas
# affinity
rasterized <- lapply(1:nrow(rasterized), function(x){
  data.frame(
      CorrelationCorr   = cor(values(corr), values(heatmaps[[x]]))
    , CorrelationPerm   = cor(values(perm), values(heatmaps[[x]]))
    , BhattacharyyaAffinityCorr = bhattacharyya_affinity(corr, heatmaps[[x]])
    , BhattacharyyaAffinityPerm = bhattacharyya_affinity(perm, heatmaps[[x]])
  )
}) %>% do.call(rbind, .) %>% cbind(rasterized, .) %>% as_tibble()

# Store results to file
write_rds(rasterized, "03_Data/03_Results/99_HeatmapMetrics.rds")

################################################################################
#### Compare Heatmaps among themselves
################################################################################
# Load all heatmaps (not bootstrapped ones)
heatmaps <- stack("03_Data/03_Results/99_RasterizedSimulations.tif")
rasterized <- read_rds("03_Data/03_Results/99_RasterizedSimulations.rds")

# Maps to compare
comp <- as.data.frame(matrix(1:nrow(rasterized), ncol = 2))
comp$steps <- unique(rasterized$steps)

# Calculate correlation for all pairs
comp$Correlation <- lapply(1:nrow(comp), function(x){
  cor(values(heatmaps[[comp[x, 1]]]), values(heatmaps[[comp[[x, 2]]]]))
}) %>% do.call(rbind, .) %>% as.numeric()

# Calculate bhattacharyya's affinity for all pairs
comp$Affinity <- lapply(1:nrow(comp), function(x){
  bhattacharyya_affinity(heatmaps[[comp[x, 1]]], heatmaps[[comp[[x, 2]]]])
}) %>% do.call(rbind, .) %>% as.numeric()

# Remove unnecessary columns
comp <- select(comp, -c(V1, V2))

# Check results
comp

# Store results
write_rds(comp, "03_Data/03_Results/99_HeatmapMetrics2.rds")
