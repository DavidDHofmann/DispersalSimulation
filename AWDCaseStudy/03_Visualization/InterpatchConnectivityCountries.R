################################################################################
#### Visualize Connections Between National Parks
################################################################################
# Description: In this script, we visualize the interpatch connectivity

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(igraph)         # For network analysis
library(rgdal)          # To load spatial data
library(sf)             # To plot spatial stuff with ggplot
library(ggspatial)      # To add scale bars etc to plots
library(rgeos)          # To calculate areas
library(ggnetwork)      # To plot network using ggplot
library(viridis)        # For nice colors
library(davidoff)       # Custom functions
library(cowplot)        # Grab ggplot legend
library(gtable)         # To combine ggplots
library(grid)           # To combine ggplots
library(gridExtra)      # To combine ggplots

################################################################################
#### Prepare Network
################################################################################
# Reload data on areas reached and country borders
# visits  <- read_rds("03_Data/03_Results/99_DirectInterpatchConnectivityCountries.rds")
visits  <- read_rds("03_Data/03_Results/99_DirectInterpatchConnectivity.rds")
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")

p100 <- ggplot(visits, aes(x = To, y = From)) +
  geom_raster(aes(fill = RelFrequency)) +
  # geom_text(aes(label = round(RelFrequency, 2))) +
  coord_fixed() +
  scale_fill_viridis() +
  scale_x_discrete(position = "top") +
  # scale_y_discrete(limits = rev) +
  theme_minimal()

p101 <- ggplot(visits, aes(x = To, y = From)) +
  geom_raster(aes(fill = MeanStepNumber)) +
  geom_text(aes(label = round(MeanStepNumber))) +
  geom_text(aes(label = paste0("(", round(SDStepNumber), ")")), nudge_y = -0.2, size = 3, fontface = 3) +
  coord_fixed() +
  scale_fill_viridis() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_minimal()

# Replace country names with integers
visits$From <- as.numeric(factor(visits$From
  , levels = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")))
visits$To <- as.numeric(factor(visits$To
  , levels = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")))

# Create a network
net <- graph_from_data_frame(
    d        = visits
  , vertices = 1:5
  , directed = T
)

# Prepare layout for countries
lay <- cbind(
    x = c(20.39, 23.94, 20.07, 25.99, 28.22)
  , y = c(-15.28, -21.80, -19.39, -14.52, -18.9)
)

################################################################################
#### Prepare Additional Plotting Data
################################################################################
# Load additional shapefiles and raster for the background
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
r       <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Rename
kaza$Name <- "KAZA-TFCA Borders"

# Prepare country labels
labels_countries <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.99, 28.22)
  , y = c(-15.28, -21.80, -19.39, -14.52, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Convert objects to sf
kaza             <- st_as_sf(kaza)
africa           <- st_as_sf(africa)
labels_countries <- st_as_sf(labels_countries)

# Convert heatmap to dataframe
r <- as.data.frame(r, xy = T)

################################################################################
#### Plot for Country Network
################################################################################
# Prepare networks for ggplotting with ggplot
net_p <- ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F)

# Prepare color palette
pal <- colorRampPalette(plasma(100, begin = 0.9, end = 0))

# Main Plot
p1 <- ggplot() +
  geom_sf(
      data        = africa
    , col         = "black"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.5
    , show.legend = F
  ) +
  geom_edges(
      data      = net_p
    , mapping   = aes(
        x    = x
      , y    = y
      , xend = xend
      , yend = yend
      , size = RelFrequency
      , col  = MeanStepNumber
    )
    , curvature = 0.1
    , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
  ) +
  geom_sf_text(
      data     = labels_countries
    , mapping  = aes(label = Label)
    , col      = "black"
    , fontface = 2
    , size     = 5
    , nudge_y  = 0.5
  ) +
  scale_size_area(
      name     = "Relative Frequency"
    , max_size = 2
  ) +
  scale_color_gradientn(
      colors  = pal(100)
    , guide   = guide_colorbar(
        title          = "Duration (Steps)"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(3.0, "cm")
      , order = 1
    )
  ) +
  scale_fill_manual(
    values = c("#70ab70", "#d9f0d3")
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    , title    = "Interpatch Connectivity"
    , subtitle = "In Relation to Dispersal Duration"
  ) +
  guides(
      size  = guide_legend(title.position = "top", order = 2)
  ) +
  theme(
      legend.position      = "bottom"
    , legend.box           = "horizontal"
    , legend.title.align   = 0.5
    , panel.background     = element_blank()
    , panel.border         = element_rect(colour = "black", fill = NA, size = 1)
    , legend.title         = element_text(size = 10),
    , legend.text          = element_text(size = 8)
    , legend.margin        = margin(c(0, 0, 0, 0))
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("black", "white")
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
  )

# Plot for the separate legend of the dots
p2 <- ggplot() +
  geom_point(
      data    = net_p %>% dplyr::select(x, y, name, Simulations) %>% distinct() %>% na.omit()
    , mapping = aes(x = x, y = y, size = Simulations, color = Simulations)
    , col     = "orange"
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    , title    = "Interpatch Connectivity"
    , subtitle = "In Relation to Dispersal Duration"
  ) +
  guides(
    size = guide_legend(title = "Number of Simulations", title.position = "top")
  ) +
  theme(
      legend.position      = "bottom"
    , legend.box           = "horizontal"
    , legend.title.align   = 0.5
    , panel.background     = element_rect(fill = "transparent", colour = NA)
    , panel.grid           = element_blank()
    , legend.margin        = margin(c(0, 0, 13, 0))
    , legend.title         = element_text(size = 10),
    , legend.text          = element_text(size = 8)
  )

# Extract legends from all plots
legend1 <- get_legend(p1)
legend2 <- get_legend(p2)

# Remove the original legend from main plot
p3 <- p1 + theme(legend.position = "none")

# Add circles to main plot
g1 <- ggplotGrob(p3)
g2 <- ggplotGrob(p2)
g2 <- gtable_filter(g2, "panel")
pos <- c(subset(g1$layout, grepl("panel", g1$layout$name), select = t:r))
p4 <- gtable_add_grob(g1, g2, t = pos$t, l = pos$l)
p4 <- ggplotify::as.ggplot(p4)

# Put legends together
legends <- grid.arrange(legend1, legend2, nrow = 1, widths = c(2, 1))
legends <- gtable_add_padding(legends, unit(c(0, 1, 0, 0.3), "cm"))

# Put them below the first plot
p5 <- arrangeGrob(p4, legends, heights = c(10, 1))
p5 <- ggplotify::as.ggplot(p5)

# Visualize
p5

# Store
ggsave("04_Manuscript/99_AreasReachedCountries.png", plot = p6)
