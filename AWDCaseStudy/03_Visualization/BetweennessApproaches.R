################################################################################
#### Betweenness
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)
library(tidyverse)
library(igraph)
library(spatstat)
library(davidoff)
library(ggpubr)
library(terra)
library(microbenchmark)
library(latex2exp)

# Set seed
set.seed(123)

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
simSine <- function(n = 100, var = 0.1) {
  x <- sort(runif(n, 0, 2 * pi))
  y <- sin(x)
  y <- y + rnorm(n, 0, var)
  df <- data.frame(x = x, y = y)
  return(df)
}

# Function to simulate a path using the above two functions
simPath <- function(n_circle = 20, n_sine = 20) {
  l1 <- simCircle(n = n_circle, x = -0.25, y = 0, r = 0.5)
  l2 <- simSine(n = n_sine, var = 0.05)
  l3 <- simCircle(n = n_circle, x = 2 * pi + 0.25, y = 0, r = 0.5)
  l <- rbind(l1, l2, l3)
  return(l)
}

# Function to interpolate spatial coordinates (we use a slightly different
# interpolation function in the true analysis)
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

################################################################################
#### Example Using a Single Path
################################################################################
# Let's simulate a path
path <- simPath(n_circle = 5, n_sine = 5)

# Generate a raster
r   <- raster(xmn = -1.5, xmx = 7.5, ymn = -1.2, ymx = 1.2, ncol = 50, nrow = 50)
r[] <- 1:ncell(r)

# Visualize it
ggplot(path, aes(x = x, y = y)) +
  geom_raster(data = as.data.frame(r, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_path(size = 0.2) +
  geom_point(size = 0.2) +
  scale_fill_gradientn(colours = terrain.colors(100)) +
  theme_minimal() +
  coord_equal()

# Let's get the betweenness map for this path
par(mfrow = c(2, 2))
plot(betweennessMap(path, r, singlecount = T, along_line = F))
lines(y ~ x, path)
plot(betweennessMap(path, r, singlecount = T, along_line = T))
plot(betweennessMap(path, r, singlecount = T, along_line = F, interpolate = T, eps = 0.01))

################################################################################
#### Example Using Multiple Paths
################################################################################
# Let's simulate some paths
paths <- lapply(1:5, function(x) {
  path <- simPath(n_circle = 10, n_sine = 20)
  path$id <- x
  return(path)
}) %>% do.call(rbind, .)

# Generate a raster
r1   <- raster(xmn = -1.5, xmx = 7.5, ymn = -1.2, ymx = 1.2, ncol = 10, nrow = 5)
r1[] <- 1:ncell(r1)
r2   <- raster(xmn = -1.5, xmx = 7.5, ymn = -1.2, ymx = 1.2, ncol = 30, nrow = 10)
r2[] <- 1:ncell(r2)

# Let's get the betweenness maps for the paths using different approaches
map1 <- betweennessMap(paths, id = "id", r1, singlecount = T, along_line = F)
map2 <- betweennessMap(paths, id = "id", r1, singlecount = T, along_line = T)
map3 <- betweennessMap(paths, id = "id", r1, singlecount = T, along_line = F, interpolate = T, eps = 0.01)
map4 <- betweennessMap(paths, id = "id", r2, singlecount = T, along_line = F)
map5 <- betweennessMap(paths, id = "id", r2, singlecount = T, along_line = T)
map6 <- betweennessMap(paths, id = "id", r2, singlecount = T, along_line = F, interpolate = T, eps = 0.01)

# Compare computational speeds of the different approaches
results <- microbenchmark(
    map1  = betweennessMap(paths, id = "id", r1, singlecount = T, along_line = F)
  , map2  = betweennessMap(paths, id = "id", r1, singlecount = T, along_line = T)
  , map3  = betweennessMap(paths, id = "id", r1, singlecount = T, along_line = F, interpolate = T, eps = 0.01)
  , map4  = betweennessMap(paths, id = "id", r2, singlecount = T, along_line = F)
  , map5  = betweennessMap(paths, id = "id", r2, singlecount = T, along_line = T)
  , map6  = betweennessMap(paths, id = "id", r2, singlecount = T, along_line = F, interpolate = T, eps = 0.01)
  , times = 10
)

# Check the results
summary(results)

# Compare mean and standard deviations
comptimes <- results %>%
  mutate(time = time / 1e6) %>%
  group_by(expr) %>%
  summarize(
    mean = mean(time)
  , sd   = sd(time)
)

# Plot all
p1 <- ggplot(as.data.frame(r1, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  geom_text(aes(label = layer), size = 1.5) +
  geom_path(data = paths, inherit.aes = F, aes(x = x, y = y, col = as.factor(id)), size = 0.5) +
  geom_point(data = paths, inherit.aes = F, aes(x = x, y = y, col = as.factor(id)), size = 1.5) +
  scale_fill_gradientn(colors = terrain.colors(100)) +
  scale_color_discrete() +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle("Simulated movement + coarse network grid")
p2 <- ggplot(as.data.frame(map1, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(parse(text = TeX(paste0(
      "Point extraction ("
    , round(comptimes$mean[[1]])
    , " \\pm "
    , round(comptimes$sd[[1]])
    , " ms)"
  ))))
p3 <- ggplot(as.data.frame(map2, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(parse(text = TeX(paste0(
      "Line extraction ("
    , round(comptimes$mean[[2]])
    , " \\pm "
    , round(comptimes$sd[[2]])
    , " ms)"
  ))))
p4 <- ggplot(as.data.frame(map3, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(parse(text = TeX(paste0(
      "Interpolated point extraction ("
    , round(comptimes$mean[[3]])
    , " \\pm "
    , round(comptimes$sd[[3]])
    , " ms)"
  ))))
p5 <- ggplot(as.data.frame(r2, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  geom_text(aes(label = layer), size = 1.5) +
  geom_path(data = paths, inherit.aes = F, aes(x = x, y = y, col = as.factor(id)), size = 0.5) +
  geom_point(data = paths, inherit.aes = F, aes(x = x, y = y, col = as.factor(id)), size = 1.5) +
  scale_fill_gradientn(colors = terrain.colors(100)) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle("Simulated movement + detailed network grid")
p6 <- ggplot(as.data.frame(map4, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(parse(text = TeX(paste0(
      "Point extraction ("
    , round(comptimes$mean[[4]])
    , " \\pm "
    , round(comptimes$sd[[4]])
    , " ms)"
  ))))
p7 <- ggplot(as.data.frame(map5, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(parse(text = TeX(paste0(
      "Line extraction ("
    , round(comptimes$mean[[5]])
    , " \\pm "
    , round(comptimes$sd[[5]])
    , " ms)"
  ))))
p8 <- ggplot(as.data.frame(map6, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  ggtitle(parse(text = TeX(paste0(
      "Interpolated point extraction ("
    , round(comptimes$mean[[6]])
    , " \\pm "
    , round(comptimes$sd[[6]])
    , " ms)"
  ))))

# Also generate a common legend
legend <- ggplot(as.data.frame(map1, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  # scale_fill_viridis_c(option = "magma") +
  scale_fill_viridis_c(option = "magma"
    , breaks = c(min(values(map1)), max(values(map1)))
    , labels = c("Low", "High")
    , guide   = guide_colorbar(
      , title          = expression("Betweenness Score")
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  theme(legend.position = "bottom", legend.box = "vertical")
legend <- get_legend(legend)

# Arrange plots
p9  <- ggarrange(p1, p2, p3, p4, ncol = 1, align = "hv", widths = c(-1, -1, -1, -1), labels = c("a1", "a2", "a3", "a4"))
p10 <- ggarrange(p5, p6, p7, p8, ncol = 1, align = "hv", widths = c(-1, -1, -1, -1), labels = c("b1", "b2", "b3", "b4"))
p   <- ggarrange(p9, p10, ncol = 2, widths = c(-1, -1))
p   <- ggarrange(p, legend, heights = c(0.9, 0.1), ncol = 1)

# Store graph to file
ggsave(p, filename = "04_Manuscript/99_BetweennessApproaches.png", scale = 2, bg = "white", width = 5, height = 3.5)

################################################################################
#### Comparing two slightly different methods to interpolate points
################################################################################
# Function 1 to interpolate
interpolatePath1 <- function(x, y, eps = 0.1) {
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

# Function 2 to interpolate
interpolatePath2 <- function(x, y, eps = 0.1) {
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

# Simulate a path and visualize the two approaches
path <- simPath(n_circle = 10, n_sine = 10)
par(mfrow = c(2, 1))
plot(y ~ x, path, type = "l")
points(interpolatePath1(path$x, path$y, eps = 1), type = "o", col = "red", pch = 20)
plot(y ~ x, path, type = "l")
points(interpolatePath2(path$x, path$y, eps = 1), type = "o", col = "red", pch = 20)

# Clearly the second approach is more desirable
