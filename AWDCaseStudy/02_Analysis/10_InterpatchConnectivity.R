################################################################################
#### Number of Connections Between Protected Areas
################################################################################
# Description: In this script, we analyze the simulated dispersal trajectories
# and identify the number of connections between protected areas

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(pbmcapply)      # To run on multiple cores with progress bar
library(rgeos)          # For spatial ananylsis
library(davidoff)       # Custom functions
library(rgdal)          # To load spatial data
library(sf)             # To plot spatial stuff
library(igraph)         # To plot networks
library(ggpubr)         # To put multiple plots together

################################################################################
#### Rasterize National Parks
################################################################################
# Load protected areas and subset to national parks only
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
prot <- subset(prot, Desig == "National Park")

# Load reference raster and rasterize all national parks
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")
prot$ID <- 1:nrow(prot)
prot_r <- raster(terra::rasterize(x = vect(prot), y = rast(r), field = "ID"))

# Visualize them
plot(prot_r, main = "National Parks", horizontal = T)
plot(prot, add = T)
text(prot, "Name", cex = 0.5, halo = T)

################################################################################
#### Prepare Simulations
################################################################################
# Load simulations. We'll focus on simulations initiated in the main study area
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")
sims <- subset(sims, Area == "Main")
sims <- dplyr::select(sims, x, y, TrackID, StepNumber)

# Identify the origin (source area) of each trajectory. For this, we create
# SpatialPoints from the first location of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )
plot(first, pch = ".")

# Assess from which protected area each trajectory leaves and add this
# information to the simulated tracks
first$SourceArea <- unlist(over(first, prot)["ID"], use.names = F)
sims$SourceArea <- first$SourceArea[match(sims$TrackID, first$TrackID)]

# Remove tracks that start outside national parks
sims <- subset(sims, !is.na(SourceArea))

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# Identify through which national parks the dispersers moved
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , Prot       = raster::extract(prot_r, sims)
)

# Calculate for each step the distance to the first coordinate. We'll use this
# to determine how far an individual had to disperse before reaching another
# national park.
visits <- visits %>%
  nest(data = -TrackID) %>%
  mutate(data = pbmclapply(data
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {

      # Project coordinates
      coords <- reprojCoords(
          xy   = x[, c("x", "y")]
        , from = CRS("+init=epsg:4326")
        , to   = CRS("+init=epsg:32734")
      )

      # Compute distance to first coordinate
      first <- coords[1, ]
      distance <- sqrt((coords[, 1] - first[1]) ** 2 + (coords[, 2] - first[2]) ** 2)
      x$DistanceFromFirst <- distance

      # Return the resulting object
      return(x)
  })) %>%
  unnest(data)

# Add the information on the reached parks and distance traveled to the
# simulations
sims$CurrentPark       <- visits$Prot
sims$DistanceFromFirst <- visits$DistanceFromFirst

# Convert simulations to regular dataframe
sims <- as.data.frame(sims, xy = T)
sims$xy <- NULL

# Each national park belongs to a country, let's assign the repsective country
# to the simulations too
sims$SourceAreaCountry <- as.character(prot$Country[match(sims$SourceArea, prot$ID)])
sims$CurrentParkCountry <- as.character(prot$Country[match(sims$CurrentPark, prot$ID)])

################################################################################
#### Identify Relative Number of Trajectories Reaching another National Park
################################################################################
# We want to determine the relative number of trajectories that successfully
# move into another national park. We do this by country. Moreover, we change
# the minimal distance considered to see how the number of reached national
# parks decreases as the minimal distance is increased.

# Compute how many dispersers were intiated within in each national parks
number_simulated_park <- sims %>%
  dplyr::select(TrackID, SourceArea) %>%
  distinct() %>%
  count(SourceArea) %>%
  setNames(c("SourceArea", "NumberSimulations"))
number_simulated_park$FromPark     <- prot$Name[match(number_simulated_park$SourceArea, prot$ID)]
number_simulated_park$FromCountry  <- prot$Country[match(number_simulated_park$SourceArea, prot$ID)]

# Summarize by country
number_simulated_country <- number_simulated_park %>%
  group_by(FromCountry) %>%
  summarize(NumberSimulations = sum(NumberSimulations))

# Function to determine how many of the simulated individuals (in percent)
# reached another national park, as well as the average dispersal duration
# required for this
getConnections <- function(min_distance) {
  reached <- sims %>%
      subset(DistanceFromFirst >= min_distance & CurrentPark != SourceArea) %>%
      group_by(TrackID, SourceAreaCountry, CurrentPark) %>%
      summarize(
          StepNumber = min(StepNumber)
        , .groups    = "drop"
      ) %>%
      group_by(SourceAreaCountry) %>%
      summarize(
          MeanStepNumber = mean(StepNumber)
        , SDStepNumber   = sd(StepNumber)
        , Frequency      = length(unique(TrackID))
        , .groups        = "drop"
      ) %>%
      left_join(number_simulated_country, by = c("SourceAreaCountry" = "FromCountry")) %>%
      mutate(PercentReachedOtherNationalPark = Frequency / NumberSimulations) %>%
      dplyr::select(Country = SourceAreaCountry, PercentReachedOtherNationalPark, MeanStepNumber, SDStepNumber)
  return(reached)
}

# Try it
getConnections(0)

# Run the function for different distances to see how the success rate
# deacreases if one only consideres more distant parks
results <- tibble(MinDistance = seq(0, 600 * 1000, length.out = 25))
results$Results <- pbmclapply(results$MinDistance
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {
      conns <- getConnections(x)
    return(conns)
})

# Convert meters to km
results$MinDistance <- results$MinDistance / 1000

# Unnest and convert fractions to percentages
reached <- unnest(results, Results)
reached$PercentReachedOtherNationalPark <- reached$PercentReachedOtherNationalPark * 100

# Plot the share of trajectories reaching other national parks in relation to
# the distance considered
p1 <- ggplot(reached, aes(x = MinDistance, y = PercentReachedOtherNationalPark, col = Country)) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(name = "Country of Origin") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 600, by = 100)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  xlab("Minimum Distance Considered (km)") +
  ylab("% Trajectories Reaching\nanother National Park") +
  theme(
      panel.grid.major = element_line(colour = "gray90", size = 0.1)
    , panel.grid.minor = element_line(colour = "gray90", size = 0.1)
    , legend.position = c(0.85, 0.7)
  )

# Also plot the average duration required to make those connections
p2 <- ggplot(reached, aes(x = MinDistance, y = MeanStepNumber, col = Country)) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(name = "Country of Origin") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 600, by = 100)) +
  scale_y_continuous(breaks = seq(0, 2000, by = 250), labels = function(x){format(x, big.mark = "'")}) +
  xlab("Minimum Distance Considered (km)") +
  ylab("Mean Dispersal Duration (Steps)\nbeforeReaching another National Park") +
  theme(
      panel.grid.major = element_line(colour = "gray90", size = 0.1)
    , panel.grid.minor = element_line(colour = "gray90", size = 0.1)
    , legend.position = "none"
  )

# Put the plots together
p <- ggarrange(p1, p2, nrow = 2, labels = c("a", "b"), label.y = 1.02)
p

# Store the plot
ggsave("04_Manuscript/99_AreasReached.png", plot = p, width = 8, height = 6)

################################################################################
#### Identify Direct Connections between National Parks
################################################################################
# Note: Here we are going to look at all direct connections between national
# parks. That is, if a trajectory moves from A through B to C, we generate an
# edge between A and B and between A and C, but not between B and C.

# Identify how long it takes on average to reach the different areas
visits <- sims %>%
  rename(From = SourceArea, To = CurrentPark) %>%
  group_by(TrackID, From, To) %>%
  summarize(
      StepNumber        = min(StepNumber)
    , DistanceFromFirst = min(DistanceFromFirst)
    , .groups           = "drop"
  ) %>%
  subset(!is.na(From) & !is.nan(To) & !is.na(To)) %>%
  arrange(TrackID, StepNumber)

# Replace national park IDs with proper park names. Also identify the country
# for each of the national parks.
visits$FromPark    <- as.character(prot$Name[match(visits$From, prot$ID)])
visits$ToPark      <- as.character(prot$Name[match(visits$To, prot$ID)])
visits$FromCountry <- as.character(prot$Country[match(visits$From, prot$ID)])
visits$ToCountry   <- as.character(prot$Country[match(visits$To, prot$ID)])

# Store visits to file
write_rds(visits, "03_Data/03_Results/99_InterpatchConnectivity.rds")

# Also store the number of simulated individuals
write_rds(number_simulated_country, "03_Data/03_Results/99_NumberSimulatedCountry.rds")
write_rds(number_simulated_park, "03_Data/03_Results/99_NumberSimulatedPark.rds")
