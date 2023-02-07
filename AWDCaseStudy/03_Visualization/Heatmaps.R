################################################################################
#### Plot of Heatmaps
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)         # To handle raster data
library(rgdal)          # To load shapefiles
library(tidyverse)      # To wrangle data
library(davidoff)       # Custom functions
library(RColorBrewer)   # For colors
library(sf)             # To plot spatial objects with ggplot
library(ggspatial)      # For north arrow and scale bar
library(cowplot)        # To grab legends
library(ggpubr)         # To arrange ggplots

################################################################################
#### Load Required Data
################################################################################
# Load heatmaps
rasterized <- read_rds("03_Data/03_Results/99_Heatmaps.rds")

# Load shapefiles
buffer  <- readOGR("03_Data/03_Results/99_BufferArea.shp")
main    <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Subset to national parks
prot <- subset(prot, Desig == "National Park")

# Load reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Prepare country labels
labels_countries <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78, 22.63, 27.92, 28.54)
  , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81, -14.54, -17.76, -20.53)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum", "Liuwa Plains"
    , "Chizarira", "Matobo"
  ), "\nNP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Put heatmaps into a stack
heatmaps <- stack(rasterized$heatmap)

# Crop them to the extent of the main area
heatmaps <- crop(heatmaps, r)

# Prepare a merged map
merged <- heatmaps[[6]] + heatmaps[[12]]

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Convert to sf for plotting with ggplot
kaza                  <- st_as_sf(kaza)
africa                <- st_as_sf(africa)
prot                  <- st_as_sf(prot)
labels_countries      <- st_as_sf(labels_countries)
labels_nationalparks  <- st_as_sf(labels_nationalparks)

# Convert heatmap to dataframe
merged <- as.data.frame(merged, xy = T)

# Also convert the reference raster to a dataframe
r <- as.data.frame(r, xy = T)

################################################################################
#### Plot
################################################################################
# Prepare color palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# Main Plot
p1 <- ggplot() +
  geom_raster(
      data    = merged
    , mapping = aes(x = x, y = y, fill = layer)
  ) +
  geom_sf(
      data        = prot
    , col         = "gray25"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.6
  ) +
  geom_sf_text(
      data     = labels_countries
    , mapping  = aes(label = Label)
    , col      = "black"
    , fontface = 2
    , size     = 5
  ) +
  geom_sf_text(
      data     = labels_nationalparks
    , mapping  = aes(label = Label)
    , col      = "gray25"
    , fontface = 3
    , size     = 2.5
    , alpha    = 0.8
  ) +
  scale_fill_gradientn(
      colours = myPalette(100)
    , labels  = function(x){format(x, big.mark = "'")}
    # , trans   = "sqrt"
    , guide   = guide_colorbar(
      , title          = "Number of Traversing Trajectories"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  geom_sf(
      data        = kaza
    , col         = "black"
    , fill        = NA
    , lty         = 1
    , lwd         = 1
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "black"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.5
    , show.legend = F
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
    # , title    = "Heatmap"
    # , subtitle = "After 2'000 Steps"
  ) +
  theme(
      legend.position  = "bottom"
    , legend.box       = "vertical"
    , panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
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

# Plot for separate legend
p2 <- ggplot() +
  geom_raster(
      data        = merged
    , mapping     = aes(x = x, y = y, fill = layer)
    , show.legend = F
  ) +
  scale_fill_gradientn(
      colours = myPalette(100)
    , labels  = function(x){format(x, big.mark = "'")}
    , guide   = guide_colorbar(
        title          = "Number of Traversing Trajectories"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  geom_sf(
      data        = kaza
    , mapping     = aes(color = "KAZA-TFCA Borders")
    , fill        = NA
    , lty         = 1
    , lwd         = 1
    , show.legend = "line"
  ) +
  geom_sf(
      data        = africa
    , mapping     = aes(color = "Country Borders")
    , fill        = NA
    , lty         = 2
    , lwd         = 0.5
    , show.legend = "line"
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , col      = NULL
    # , title    = "Heatmap"
    # , subtitle = "After 2'000 Steps"
  ) +
  scale_color_manual(
    , breaks = c("Country Borders", "KAZA-TFCA Borders")
    , values = c("black", "black")
    , guide  = guide_legend(
        override.aes = list(
          linetype = c(2, 1)
        , shape    = c(NA, NA)
        , lwd      = c(0.5, 1)
      )
    )
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  theme(
      legend.position       = c(0.20, 0.90)
    , legend.box            = "vertical"
    , legend.background     = element_blank()
    , legend.box.background = element_rect(fill  = "white")
    , legend.margin         = margin(0, 8, 2, 6)
    , legend.text           = element_text(color = "black")
    , legend.key            = element_blank()
    , legend.key.size       = unit(0.8, "lines")
    , legend.key.width      = unit(1.2, "lines")
    , panel.background      = element_blank()
  )

# Extract legend
legend <- get_legend(p2)

# Put into main plot
p3 <- p1 + annotation_custom(
      grob = legend
    , xmin = 18.25
    , xmax = 21
    , ymin = -12.7
    , ymax = -13.5
  )

# Show it
p3

# Store it to file
ggsave("04_Manuscript/99_Heatmap.png", plot = p3, bg = "white")

################################################################################
#### Function to Plot
################################################################################
# Write a function to plot a heatmap
plotHeatmap <- function(x, subtitle = NULL, legend = T, barwidth = 10){

  # Prepare dataframe
  x <- as.data.frame(x, xy = T)
  names(x) <- c("x", "y", "layer")

  # Main Plot
  p1 <- ggplot() +
    geom_raster(
        data    = x
      , mapping = aes(x = x, y = y, fill = layer)
    ) +
    geom_sf(
        data        = prot
      , col         = "gray25"
      , fill        = NA
      , lty         = 1
      , lwd         = 0.1
      , show.legend = F
      , alpha       = 0.6
    ) +
    geom_sf_text(
        data     = labels_countries
      , mapping  = aes(label = Label)
      , col      = "black"
      , fontface = 2
      , size     = 5
    ) +
    geom_sf_text(
        data     = labels_nationalparks
      , mapping  = aes(label = Label)
      , col      = "gray25"
      , fontface = 3
      , size     = 2.5
      , alpha    = 0.8
    ) +
    scale_fill_gradientn(
        colours = myPalette(100)
      , labels  = function(x){format(x, big.mark = "'")}
      , guide   = guide_colorbar(
        , title          = "Number of Traversing Trajectories"
        , show.limits    = T
        , title.position = "top"
        , title.hjust    = 0.5
        , ticks          = T
        , barheight      = unit(0.6, "cm")
        , barwidth       = unit(barwidth, "cm")
      )
    ) +
    geom_sf(
        data        = kaza
      , mapping     = aes(color = "KAZA-TFCA Borders")
      , fill        = NA
      , lty         = 1
      , lwd         = 1
      , show.legend = F
    ) +
    geom_sf(
        data        = africa
      , mapping     = aes(color = "Country Borders")
      , fill        = NA
      , lty         = 2
      , lwd         = 0.5
      , show.legend = F
    ) +
    scale_color_manual(
      , breaks = c("National Parks", "Country Borders", "KAZA-TFCA Borders")
      , values = c("gray25", "black", "black")
      , guide  = guide_legend(
          override.aes = list(
            linetype = c(1, 1, 1)
          , shape    = c(NA, NA, NA)
          , lwd      = c(0.2, 0.5, 1)
        )
      )
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
      , title    = "Heatmap"
      , subtitle = subtitle
    ) +
    theme(
        legend.position      = "top"
      , legend.justification = "right"
      , legend.box           = "vertical"
      , legend.box.margin    = margin(-10, 0, -10, 0)
      , panel.background     = element_blank()
      , panel.border         = element_rect(colour = "black", fill = NA, size = 1)
      , plot.subtitle        = element_text(margin = margin(t = 0, b = -30))
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 0.5
      , height     = unit(0.15, "cm")
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

  if (legend){

    # Plot for separate legend
    p2 <- ggplot() +
      geom_raster(
          data        = x
        , mapping     = aes(x = x, y = y, fill = layer)
        , show.legend = F
      ) +
      scale_fill_gradientn(
          colours = myPalette(100)
        , labels  = function(x){format(x, big.mark = "'")}
        , guide   = guide_colorbar(
            title          = "Number of Traversing Trajectories"
          , show.limits    = T
          , title.position = "top"
          , title.hjust    = 0.5
          , ticks          = T
          , barheight      = unit(0.6, "cm")
          , barwidth       = unit(10.0, "cm")
        )
      ) +
      geom_sf(
          data        = kaza
        , mapping     = aes(color = "KAZA-TFCA Borders")
        , fill        = NA
        , lty         = 1
        , lwd         = 0.5
        , show.legend = "line"
      ) +
      geom_sf(
          data        = africa
        , mapping     = aes(color = "Country Borders")
        , fill        = NA
        , lty         = 2
        , lwd         = 0.2
        , show.legend = "line"
      ) +
      geom_sf(
          data        = prot
        , mapping     = aes(color = "National Parks")
        , fill        = NA
        , lty         = 1
        , lwd         = 0.5
        , show.legend = "line"
      ) +
      labs(
          x        = NULL
        , y        = NULL
        , col      = NULL
        , title    = "Heatmap"
        , subtitle = "After 2'000 Steps"
      ) +
      scale_color_manual(
        , breaks = c("National Parks", "Country Borders", "KAZA-TFCA Borders")
        , values = c("gray30", "black", "black")
        , guide  = guide_legend(
            override.aes = list(
              linetype = c(1, 2, 1)
            , shape    = c(NA, NA, NA)
            , lwd      = c(0.2, 0.2, 0.5)
          )
        )
      ) +
      coord_sf(
          crs    = 4326
        , xlim   = c(min(r$x), max(r$x))
        , ylim   = c(min(r$y), max(r$y))
        , expand = F
      ) +
      theme(
          legend.position       = c(0.20, 0.90)
        , legend.box            = "vertical"
        , legend.background     = element_blank()
        , legend.box.background = element_rect(fill  = "white")
        , legend.margin         = margin(0, 8, 2, 6)
        , legend.text           = element_text(color = "black")
        , legend.key            = element_blank()
        , legend.key.size       = unit(0.8, "lines")
        , legend.key.width      = unit(1.2, "lines")
        , panel.background      = element_blank()
      )

    # Extract legend
    legend <- get_legend(p2)

    # Put into main plot
    p3 <- p1 + annotation_custom(
          grob = legend
        , xmin = 18.35
        , xmax = 21
        , ymin = -13
        , ymax = -13.5
      )
  } else {
    p3 <- p1
  }

  # Return it
  return(p3)

}

# Let's apply the function to get all desired plots
labels <- paste0(rep(c("a", "b"), each = 6), rep(1:3, each = 2))
maps <- lapply(1:nlayers(heatmaps), function(x){
  subtitle <- paste0("After ", format(rasterized$steps[x], big.mark = "'"), " Steps")
  map <- plotHeatmap(
      x        = heatmaps[[x]]
    , subtitle = subtitle
    , legend   = F
    , barwidth = 7
  )
  map <- map + annotate("text"
    , x        = 18.8
    , y        = -13.3
    , label    = labels[x]
    , col      = "black"
    , fontface = 2
    , size     = 10
  )
  return(map)
})

# Arrange plots nicely
p <- ggarrange(maps[[2]], maps[[8]], maps[[4]], maps[[10]], maps[[6]], maps[[12]], ncol = 2, nrow = 3)

# Store the arranged plot
ggsave("04_Manuscript/99_HeatmapsIndividual.png"
  , plot   = p
  , scale  = 2
  , height = 9
  , width  = 6
  , bg     = "white"
)
