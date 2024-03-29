################################################################################
#### Connectivity Maps
################################################################################
# Description: We now make use of the simulated dispersal trajectories to derive
# a set of connectivity maps. This includes a heatmap, a betweenness map, and a
# map of inter-patch connectivity.

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(viridis)        # For nice colors
library(rgdal)          # To load shapefiles
library(rgeos)          # For geometry manipulation
library(sf)             # To plot spatial things with ggplot
library(spatstat)       # To quickly rasterize (and count) spatial lines
library(maptools)       # To quickly rasterize (and count) spatial lines
library(igraph)         # For network analysis
library(ggpubr)         # To arrange multiple ggplots
library(ggnetwork)      # To plot a network using ggplot

# Set working directory
setwd("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/DispersalSimulation/SimulationStudy")

# Load custom functions
source("00_Functions.R")

# Load covariate layers and simulated movements
cov <- read_rds("99_CovariateLayers.rds")
nps <- read_rds("99_NationalParks.rds")
sims <- read_rds("99_SimulatedMovements.rds")

# Create a polygon for the extent of the study area
ext <- extent(extent(0, 100, 0, 100))
ext <- as(ext, "SpatialPolygons")

# Assign IDs to NPs
nps$ID <- 1:nrow(nps)

################################################################################
#### Heatmap
################################################################################
# In order to generate a heatmap, we first need to convert all trajectories into
# proper spatial lines. So let's go through each simulated individual and
# convert its coordinates into a spatial line
tracks <- lapply(unique(sims$ID), function(x) {
  sub <- sims[sims$ID == x, ]
  coordinates(sub) <- c("x", "y")
  lines <- spLines(sub)
  return(lines)
})
tracks <- do.call(rbind, tracks)
tracks$ID <- 1:length(tracks)

# Visualize them
plot(cov[["elev"]])
plot(tracks, add = T, col = tracks$ID)

# Create an empty raster onto which we can rasterize the lines. I'll use the
# covariate raster as reference
heatmap <- raster(cov)

# Rasterize and count the tracks (we use a custom function that is much quicker
# than raster::rasterize(... fun = "count"))
heatmap <- rasterizeSpatstat(tracks, heatmap)

# Crop the layer to the main study area
heatmap <- crop(heatmap, ext)

# Plot the resulting heatmap
plot_heatmap <- ggplot(as.data.frame(heatmap, xy = T)) +
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
  scale_fill_gradientn(
      colors = hcl.colors(n = 100, palette = "Spectral", rev = T)
    , name   = "Traversal Frequency"
    # , trans  = "sqrt"
    , guide   = guide_colorbar(
        show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap") +
  coord_sf()
plot_heatmap

################################################################################
#### Betweenness
################################################################################
# Overlay the study area with a regular grid (we could again use the covariate
# layer as reference, but a coarser resolution will suffice for now)
grid <- raster(ncol = 100, nrow = 100, crs = NA)
grid <- setExtent(grid, ext)
grid[] <- 1:ncell(grid)
plot(grid)

# Prepare network (the center of each grid cell is a node in the final network)
vertices <- values(grid)
lay <- as.matrix(as.data.frame(grid, xy = T)[, c(1, 2)])

# Interpolate simulated data to detect cell-transitions at finer scale
sims_inter <- sims %>%
  nest(data = -c(SourceArea, ID)) %>%
  mutate(data = map(data, function(x) {
    interpolatePath(x$x, x$y, eps = 0.1)
  })) %>%
  unnest(data)

# At each coordinate of the simulated trajectories we now extract the cell IDs
# from the grid (this allows us to see from which to which cell an individual
# moved)
visits <- data.frame(
    ID          = sims_inter$ID
  , x           = sims_inter$x
  , y           = sims_inter$y
  , StepID      = sims_inter$segment_id
  , CellID      = raster::extract(grid, sims_inter[, c("x", "y")])
)

# Now we write a function that we can use to retrieve the visitation history
# along a single trajectory. That is, the function tells us all transitions from
# one cell to another. The option "singlecount = T" ensures that re-visits are
# not counted twice.
visitHist <- function(x, singlecount = T){
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections = 1
  }
  return(transitions)
}

# Try what it does!
visitHist(c(1, 2, 2, 2, 3, 1), singlecount = T)
visitHist(c(1, 2, 2, 2, 3, 1), singlecount = F)

# We now want to retrieve the visitation history for each individual. For this,
# we first nest the data on visits
visits <- visits %>% nest(data = -ID)

# Then apply the function
visits <- mutate(visits, history = map(data, function(x) {
  visitHist(x$CellID, singlecount = T)
}))

# Let's check the visitation history of one of the individuals
visits$history[[1]]

# Summarize across individuals
history <- visits %>%
  dplyr::select(ID, history) %>%
  unnest(history) %>%
  group_by(from, to) %>%
  summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
  ungroup() %>%
  mutate(weight = mean(TotalConnections) / TotalConnections)

# Use this to compute (weighted) betweenness
net <- graph_from_data_frame(history, vertices = vertices)
is.weighted(net)

# Calculate betweenness and put the resulting values into a raster
betweenness <- grid
values(betweenness) <- betweenness(net)

# Crop the layer to the main study area
betweenness <- crop(betweenness, ext)

# Plot the betweenness (we'll make low scores better visible by applying a
# square root transformation to the color scale)
plot_betweenness <- ggplot(as.data.frame(betweenness, xy = T)) +
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
  scale_fill_viridis(
      option = "magma"
    , name   = "Betweenness"
    # , trans  = "sqrt"
    , guide   = guide_colorbar(
        show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ggtitle("Betweenness") +
  coord_sf()
plot_betweenness

################################################################################
#### Inter-Patch Connectivity
################################################################################
# Rasterize the naional parks
npsr <- rasterize(nps, disaggregate(grid, fact = 10), field = "ID")

# Use the raster to identify through which areas each trajectory passes
visits <- data.frame(
    ID          = sims_inter$ID
  , x           = sims_inter$x
  , y           = sims_inter$y
  , StepID      = sims_inter$segment_id
  , SourceArea  = sims_inter$SourceArea
  , CurrentArea = raster::extract(npsr, sims_inter[, c("x", "y")])
)

# Count the number of connections from one area to another and the duration it
# takes
conns <- visits %>%
  subset(SourceArea != CurrentArea) %>%
  group_by(ID, SourceArea, CurrentArea) %>%
  summarize(Steps = min(StepID), .groups = "drop") %>%
  group_by(SourceArea, CurrentArea) %>%
  summarize(MeanSteps = mean(Steps), TotalConnections = n(), .groups = "drop")

# Generate network
net <- graph_from_data_frame(conns, vertices = nps$ID)
lay <- coordinates(gCentroid(nps, byid = T))

# Prepare networks for ggplotting with ggplot
net_p <- ggnetwork(net, layout = lay, arrow.gap = 1, scale = F)

# Plot the network
plot_interpatch <- ggplot(as.data.frame(cov[[1]], xy = T)) +
  geom_raster(aes(x = x, y = y), fill = "black") +
  geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
  geom_edges(
      data      = net_p
    , mapping   = aes(x = x, y = y, xend = xend, yend = yend , size = TotalConnections, col = MeanSteps)
    , curvature = 0.2
    , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 30)
  ) +
  geom_nodes(data = net_p, aes(x = x, y = y), col = "white") +
  scale_color_viridis_c(
      option = "plasma"
    , name   = "Duration (No. Steps)"
    , begin  = 0.3
    , breaks = c(110, 130, 150, 170)
    , guide  = guide_colorbar(
        show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(3.0, "cm")
    )
  ) +
  scale_size_area(
      name     = "Total Connections"
    , max_size = 0.5
    , breaks   = 1:5
    # , trans    = "exp"
    , guide = guide_legend(
        title.position = "top"
      , title.hjust    = 0.5
    )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ggtitle("Inter-Patch Connectivity") +
  coord_sf()
plot_interpatch

################################################################################
#### Combining all Plots
################################################################################
# Remove the legends for now and rotate y-axes
p1 <- plot_heatmap + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
p2 <- plot_betweenness + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
p3 <- plot_interpatch + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

# Arrange the plots
p <- ggarrange(p1, p2, p3)
ggsave("ConnectivityMetrics.png", width = 8, height = 8, bg = "white")
p
