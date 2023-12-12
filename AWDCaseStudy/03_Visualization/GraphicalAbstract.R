################################################################################
#### Figures for the Graphical Abstract
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse) # To wrangle data
library(raster)    # To handle spatial data
library(terra)     # To handle spatial data
library(rgeos)     # To manipulate spatial data
library(igraph)    # To create networks
library(NLMR)      # To create random landscapes
library(gdistance) # To calculate least-cost paths
library(pbmcapply) # For multicore abilities with progress bar
library(spatstat)  # To quickly rasterize lines
library(maptools)  # To quickly rasterize lines
library(viridis)   # For nice colors
library(sf)        # To plot spatial stuff with ggplot
library(ggpubr)    # To arrange multiple plots
library(ggnetwork) # To plot networks with ggplot
library(ggridges)  # For ridgeline plot
library(lemon)     # For capped coordinate system

################################################################################
#### Required Functions
################################################################################
# Function to sample locations uniformly within a circle
simCircle <- function(n = 20, x = 0, y = 0, r = 1) {
  sl <- runif(n, 0, r ** 2)
  ta <- runif(n, 0, 2 * pi)
  x_ <- sqrt(sl) * cos(ta) + x
  y_ <- sqrt(sl) * sin(ta) + y
  df <- data.frame(x = x_, y = y_)
  return(df)
}

# Function to sample locations uniformly on a sine
simSine <- function(n = 100, var = 0.1, reverse = F) {
  x <- sort(runif(n, 0, 2 * pi))
  y <- sin(x)
  y <- y + rnorm(n, 0, var)
  if (reverse) {
    x <- rev(x)
    y <- rev(y)
  }
  df <- data.frame(x = x, y = y)
  return(df)
}

# Function to simulate a path using the above two functions
simPath <- function(n_circle = 20, n_sine = 20, reverse = F, var = 0.05) {
  l1 <- simCircle(n = n_circle, x = -0.25, y = 0, r = 0.5)
  l2 <- simSine(n = n_sine, var = var, reverse = reverse)
  l3 <- simCircle(n = n_circle, x = 2 * pi + 0.25, y = 0, r = 0.5)
  if (reverse) {
    l <- rbind(l3, l2, l1)
  } else {
    l <- rbind(l1, l2, l3)
  }
  return(l)
}

# Function to interpolate spatial coordinates (we use a slightly different
# interpolation function in the true analysis)
interpolatePath <- function(x, y, eps = 0.1) {
  ppp <- psp(
      x0 = x[-length(x)]
    , x1 = x[-1]
    , y0 = y[-length(y)]
    , y1 = y[-1]
    , window = owin(xrange = range(x), yrange = range(y))
  )
  ppp_int <- pointsOnLines(ppp, eps = eps, shortok = F)
  ppp_int <- as.data.frame(ppp_int)
  return(ppp_int)
}

# Function to generate the visitation history
visitHist <- function(path, r, singlecount = F, along_line = F, interpolate = F, eps = NULL) {

  # Interpolate path if desired
  if (interpolate) {
    path <- interpolatePath(path$x, path$y, eps = eps)
  }

  # Create spatial line from coordinates if desired, then extract raster values
  if (along_line) {
      path_sp <- spLines(cbind(path$x, path$y))
      extracted <- raster::extract(r, y = path_sp, along = T)
      extracted <- unname(extracted[[1]])
    } else {
      extracted <- raster::extract(r, y = cbind(path$x, path$y))
  }

  # Count cell transitions
  visits <- extracted %>%
    data.frame(from = lag(.), to = .) %>%
    na.omit() %>%
    subset(from != to) %>%
    group_by(from, to) %>%
    summarize(count = n(), .groups = "drop")

  # Adjust count if necessary
  if (singlecount) {
    visits$count <- 1
  }

  # Return the visits
  return(visits)
}

# Function to compute betweenness map
betweennessMap <- function(paths, r, id = NULL, singlecount = F, along_line = F, interpolate = F, eps = NULL) {
  if (is.null(id)) {
    paths$id <- "dummy"
  }
  paths_nested <- paths %>% nest(Data = -all_of(id))
  histories <- lapply(paths_nested$Data, function(x) {
    x %>%
      visitHist(r, singlecount, along_line, interpolate = interpolate, eps = eps) %>%
      group_by(from, to) %>%
      summarize(count = n(), .groups = "drop")
  })
  history <- histories %>%
    do.call(rbind, .) %>%
    group_by(from, to) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    mutate(weight = mean(count) / count)
  net <- graph_from_data_frame(history, vertices = unique(values(r)))
  bet <- r
  values(bet) <- betweenness(net)
  return(bet)
}

# Function to compute heatmap
heatMap <- function(paths, r, id = NULL) {
  if (is.null(id)) {
    paths$id <- "dummy"
  }
  lines <- lapply(unique(paths$id), function(x) {
    sub <- paths[paths$id == x, ]
    lines <- SpatialPoints(as.matrix(sub[, c("x", "y")]))
    lines <- spLines(lines)
    lines$ID <- x
    return(lines)
  })
  lines <- do.call(rbind, lines)
  values(r) <- 0
  im <- as.im.RasterLayer(r)
  summed <- im
  for (y in 1:length(lines)){
    line    <- as.psp(lines[y, ], window = im)
    line    <- as.mask.psp(line)
    line_r  <- as.im.owin(line, na.replace = 0)
    summed  <- Reduce("+", list(summed, line_r))
  }
  return(raster(summed))
}

################################################################################
#### Simulating Paths and Generating Maps
################################################################################
# Let's simulate some paths, going both directions
set.seed(123)
paths <- lapply(1:50, function(x) {
  reverse <- rbernoulli(n = 1, p = 0.2)
  path    <- simPath(n_circle = 20, n_sine = 10, var = 0.1, reverse = reverse)
  path$id <- x
  return(path)
}) %>% do.call(rbind, .)

# Generate a raster
n <- 50
r <- raster(
    xmn  = min(paths$x) - 0.5
  , xmx  = max(paths$x) + 0.5
  , ymn  = min(paths$y) - 0.5
  , ymx  = max(paths$y) + 0.5
  , ncol = n * 3
  , nrow = n
)
r[] <- 1:ncell(r)

# Compute heatmap and betweenness map
heat <- heatMap(paths, r = r, id = "id")
betw <- betweennessMap(paths, r = r, id = "id", singlecount = T, interpolate = T, eps = 0.01)

# Apply some transformations for nicer plots
heat <- heat %>%
  disaggregate(method = "bilinear", fact = 2) %>%
  aggregate(fact = 2, fun = max)
betw <- betw %>%
  disaggregate(method = "bilinear", fact = 2) %>%
  aggregate(fact = 2, fun = max) %>%
  sqrt()

# Plot of heatmap
p1 <- ggplot() +
  geom_raster(
      data    = as.data.frame(heat, xy = T)
    , mapping = aes(x = x, y = y, fill = layer)
  ) +
  scale_fill_gradientn(
      colours = c("white", magma(100))
    , guide   = guide_colorbar(
      , title          = expression("Traversal Frequency")
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  coord_equal() +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
  ) +
  theme_void() +
  theme(legend.position  = "none")

# Plot of betweenness
p2 <- ggplot() +
  geom_raster(
      data    = as.data.frame(betw, xy = T)
    , mapping = aes(x = x, y = y, fill = layer)
  ) +
  scale_fill_gradientn(
      colours = c("white", magma(100))
    # , trans   = "sqrt"
    , guide   = guide_colorbar(
      , show.limits = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
    , na.value = "transparent"
  ) +
  coord_equal() +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
  ) +
  theme_void() +
  theme(legend.position  = "none")
ggarrange(p1, p2, ncol = 1)


################################################################################
#### Plot of Arbitrary Movement Model
################################################################################
# Generate some arbitrary data to plot
x <- seq(-0.2, 0.2, by = 0.001)
y <- dnorm(x, mean = 0, sd = 0.05)
ci_90 <- qnorm(c(0.050, 0.950), mean = 0, sd = 0.05)
ci_95 <- qnorm(c(0.025, 0.975), mean = 0, sd = 0.05)
ci_99 <- qnorm(c(0.005, 0.995), mean = 0, sd = 0.05)

# Shift data to three different means
dat1 <- data.frame(Covariate = "v", Mean = -0.2, x = x - 0.2, y = y)
dat2 <- data.frame(Covariate = "w", Mean = 0.4, x = x + 0.4, y = y)
dat3 <- data.frame(Covariate = "x", Mean = -0.3, x = x - 0.3, y = y)
dat4 <- data.frame(Covariate = "y", Mean = -0.2, x = x - 0.2, y = y)
dat5 <- data.frame(Covariate = "z", Mean = 0.2, x = x + 0.2, y = y)

# Put all together
dat <- rbind(dat1, dat2, dat3, dat4, dat5)

# Calculate summaries
averaged <- data.frame(
    Covariate = c("v", "w", "x", "y", "z")
  , Mean = c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , LCI_90 = ci_90[1] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , UCI_90 = ci_90[2] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , LCI_95 = ci_95[1] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , UCI_95 = ci_95[2] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , LCI_99 = ci_99[1] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , UCI_99 = ci_99[2] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
)

# Generate plot
p3 <- ggplot(dat, aes(x = x, y = Covariate, height = y, group = Covariate)) +
  geom_errorbarh(
      data        = averaged
    , inherit.aes = F
    , height      = 0
    , size        = 0.5
    , alpha       = 1
    , col         = "orange"
    , aes(y = Covariate, xmin = LCI_99, xmax = UCI_99)
  ) +
  geom_errorbarh(
      data        = averaged
    , inherit.aes = F
    , height      = 0
    , size        = 1
    , alpha       = 1
    , col         = "orange"
    , aes(y = Covariate, xmin = LCI_95, xmax = UCI_95)
  ) +
  geom_errorbarh(
      data        = averaged
    , inherit.aes = F
    , height      = 0
    , size        = 2
    , alpha       = 0.5
    , col         = "orange"
    , aes(y = Covariate, xmin = LCI_90, xmax = UCI_90)
  ) +
  geom_point(
      data        = averaged
    , inherit.aes = F
    , size        = 3
    , col         = "orange"
    , aes(x = Mean, y = Covariate)
  ) +
  geom_vline(xintercept = 0, lty = 2, col = "gray50") +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , xlim   = c(-0.5, 0.5)
  ) +
  labs(x = expression(beta*"-Coefficient"))

# Show plots
p1 <- p1 +
  labs(title = "Heatmap", subtitle = "(Traversal Frequency)") +
  theme(plot.title = element_text(face = 2, hjust = 0.5), plot.subtitle = element_text(face = 3, hjust = 0.5))
p2 <- p2 +
  labs(title = "Betweenness", subtitle = "(Bottlenecks & Dispersal Corridors)") +
  theme(plot.title = element_text(face = 2, hjust = 0.5), plot.subtitle = element_text(face = 3, hjust = 0.5))
p3 <- p3 +
  labs(title = "Integrated Step\nSelection Model") +
  theme(plot.title = element_text(face = 2, hjust = 0.5), plot.subtitle = element_text(face = 3, hjust = 0.5))

# Store the plots
ggsave("test.png", plot = p1)
ggsave("test2.png", plot = p2)
ggsave("test3.png", plot = p4, height = 3, width = 4)
