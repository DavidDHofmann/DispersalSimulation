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

# # Set working directory
# setwd("/home/david/ownCloud/DispersalSimulation")

# Load custom functions
source("00_Functions.R")

# Load covariate layers and simulated movements
cov <- read_rds("99_CovariateLayers.rds")
nps <- read_rds("99_NationalParks.rds")
sims <- read_rds("99_SimulatedMovements.rds")

# Create a polygon for the extent of the study area
ext <- extent(extent(0, 100, 0, 100))
ext <- as(ext, "SpatialPolygons")

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
  scale_fill_viridis(option = "magma", name = "Traversal Frequency") +
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

# At each coordinate of the simulated trajectories we now extract the cell IDs
# from the grid (this allows us to see from which to which cell an individual
# moved)
visits <- data.frame(
    ID          = sims$ID
  , step_number = sims$step_number
  , x           = sims$x
  , y           = sims$y
  , CellID      = raster::extract(grid, sims[, c("x", "y")])
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
visits <- mutate(visits, history = map(data, function(x){
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
  scale_fill_viridis(option = "magma", name = "Betweenness", trans = "sqrt") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ggtitle("Betweenness") +
  coord_sf()
plot_betweenness

################################################################################
#### Inter-Patch Connectivity
################################################################################
# Finally, we want to visualize the relative frequency at which simulated
# individuals are moving between national parks. For this, let's give each
# national park an ID.
nps$ID <- 1:nrow(nps)

# In order to determine interpatch connectivity between national parks, we need
# to know from which national park trajectories leave and into which other
# national parks they go to. We can use the first coordinate of each simulated
# trajectory to determine from which national park the trajectory leaves
first <- subset(sims, step_number == 1)
coordinates(first) <- c("x", "y")
plot(first, pch = 20, cex = 0.3)

# Determine with which national park each startpoint intersects (can only be
# one)
from <- gIntersects(nps, first, byid = T)
from <- apply(from, 1, which)
from <- as.vector(from)

# Let's also check whith which national parks each trajectory intersects
to <- gIntersects(nps, tracks, byid = T)
to <- as.data.frame(to)
names(to) <- 1:3

# Put all into a single dataframe
inter <- cbind(first$ID, from, to)
names(inter)[1:2] <- c("ID", "from")

# Let's count how many individuals were released at each national park. Here,
# this is not too interesting as we put the same number of dispersers in each
# national park.
inter %>% count(from)

# Let's count the number of connections from one park to another
conns <- inter %>%
  gather(key = to, value = reached, 3:ncol(.)) %>%
  subset(reached & from != to) %>%
  group_by(from, to) %>%
  summarize(total_connections = n(), .groups = "drop")

# Generate network
net <- graph_from_data_frame(conns, vertices = nps$ID)
lay <- coordinates(gCentroid(nps, byid = T))

# Prepare networks for ggplotting with ggplot
net_p <- ggnetwork(net, layout = lay, arrow.gap = 1, scale = F)

# Plot the network
plot_interpatch <- ggplot(as.data.frame(cov[[1]], xy = T)) +
  geom_raster(aes(x = x, y = y), fill = "black") +
  geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
  geom_edges(data = net_p, aes(x = x, y = y, xend = xend, yend = yend
    , size = total_connections), color = magma(20)[15], curvature = 0.2
    , arrow = arrow(length = unit(6, "pt"), type = "closed", angle = 30)) +
  geom_nodes(data = net_p, aes(x = x, y = y), col = "white") +
  scale_size_area(name = "Total Connections", max_size = 0.5, breaks = 1:5
    , trans = "exp") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ggtitle("Inter-Patch Connectivity") +
  coord_sf()
plot_interpatch

################################################################################
#### Combining all Plots
################################################################################
# Remove the legends for now and rotate y-axes
p1 <- plot_heatmap + theme(legend.position = "none"
  , axis.title.y = element_text(angle = 0, vjust = 0.5))
p2 <- plot_betweenness + theme(legend.position = "none"
  , axis.title.y = element_text(angle = 0, vjust = 0.5))
p3 <- plot_interpatch + theme(legend.position = "none"
  , axis.title.y = element_text(angle = 0, vjust = 0.5))

# Arrange the plots
p <- ggarrange(p1, p2, p3)
p

# Let's store the arranged plots to file
ggsave(plot = p, "ConnectivityPlots.png", bg = "white", width = 8, height = 8)
