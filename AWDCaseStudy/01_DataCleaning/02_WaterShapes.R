################################################################################
#### Preparation of Major Water Areas
################################################################################
# Description: Preparation of a shapefile of all major water areas. This might
# be handy for a nice plot

# Clean environment
rm(list = ls())

# Load required packages
library(rgdal)  # To load spatial data

# Load water shapefiles
water <- readOGR("03_Data/01_RawData/GEOFABRIK/Water.shp")

# Specify the areas that we want to keep for plotting later
object1 <- water[grepl(water@data$name, pattern = "Okavango.*Delta"), ][2, ]
object2 <- water[grepl(water@data$name, pattern = "Linyanti.*Delta"), ]
object3 <- water[grepl(water@data$name, pattern = "Garangwe.*Pan"), ][2, ]
object4 <- water[grepl(water@data$name, pattern = "Lake.*Kariba"), ]
object5 <- water[grepl(water@data$name, pattern = "Lake.*Ngami"), ]
object6 <- water[grepl(water@data$name, pattern = "Barotse"), ]
object7 <- water[grepl(water@data$name, pattern = "Lukanga"), ]
object8 <- water[grepl(water@data$name, pattern = "Kafue"), ]

# Put all objects together
water <- rbind(object1, object2, object3, object4, object5, object6, object7, object8)

# Visualize them
plot(water, col = "lightblue", border = "black", lwd = 0.4)

# Store the file
writeOGR(
    water
  , dsn       = "03_Data/02_CleanData"
  , layer     = "03_LandscapeFeatures_MajorWaters_GEOFABRIK"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
