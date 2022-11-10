library(sf)
library(raster)

## grids used in the National Park Hoge Kempen:
# total number of sites in the 300 m * 300 m finest resolution layer
n_sites <- 36^2 
grid_extent <- 
  as(as(raster::extent(0, 10800, 0, 10800), "SpatialPolygons"), "sf")

# grid 36*36 (resolution 300 m * 300 m)
grid_36x36 <- 
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(n_sites)),
  site_id = as.integer(factor(1:n_sites)))

# grid 18*18 (resolution 600 m * 600 m)
grid_18x18 <- 
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(n_sites/4)),
  site_id = as.integer(factor(1:(n_sites/4))))

# grid 12*12 (resolution 900 m * 900 m)
grid_12x12 <- 
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(n_sites/9)),
  site_id = as.integer(factor(1:(n_sites/9))))

# grid 9*9 (resolution 1200 m * 1200 m)
# use for small home range size
grid_9x9 <- 
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(n_sites/16)),
  site_id = as.integer(factor(1:(n_sites/16))))

# grid 6*6 (resolution 1800 m * 1800 m)
# use for small home range size
grid_6x6 <- 
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(sqrt(n_sites))),
  site_id = as.integer(factor(1:(sqrt(n_sites)))))

# grid 5*5 (resolution 2120 * 2160 m)
# use for normal home range size
grid_5x5 <- 
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(n_sites/ 51.84)),
  site_id = as.integer(factor(1:(n_sites/51.84))))

# grid 4*4 (resolution 2700 m * 2700 m)
# use for large home range size
grid_4x4 <- 
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(n_sites/ 81)),
  site_id = as.integer(factor(1:(n_sites/81))))

# grid 3*3 (resolution 3600 m * 3600 m)
# use for large home range size
grid_3x3 <-
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(n_sites/ 144)),
    site_id = as.integer(factor(1:(n_sites/144))))

# grid 2*2 (resolution 5400 m * 5400 m)
# use for large home range size
grid_2x2 <-
  st_as_sf(
    st_make_grid(grid_extent, square = T, n = sqrt(n_sites/ 324)),
    site_id = as.integer(factor(1:(n_sites/324))))

grid.list <- list(
  g36 <- grid_36x36,
  g18 <- grid_18x18,
  g12 <- grid_12x12,
  g9 <- grid_9x9,
  g6 <- grid_6x6,
  g5 <- grid_5x5,
  g4 <- grid_4x4,
  g3 <- grid_3x3,
  g2 <- grid_2x2,
  g1 <- grid_extent
)

names(grid.list) <- c("g36", "g18", "g12", "g9", "g6", 
                      "g5", "g4", "g3", "g2", "g1")

# save list of grids
saveRDS(grid.list, "data/grids.rds")

