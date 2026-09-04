# reach_slope(): slope of a short surveyed reach from the DEM.
#
# Takes the transect as an sf LINESTRING. Both ends are snapped to the nearest
# DEM stream cell (within snap_m) so a few meters of GPS error don't put an
# endpoint on the bank, then the DEM is sampled along the line.
# slope = elevation drop / length (m/m; x100 for percent).

library(sf)
library(terra)

reach_slope <- function(reach, dem, streams, snap_m = 10) {

  streams <- rast(streams)
  coords <- st_coordinates(st_transform(st_geometry(reach), crs(dem)))[, 1:2]

  snap <- function(xy) {
    circle <- vect(st_buffer(st_sfc(st_point(xy), crs = crs(dem)), snap_m))
    cells <- as.data.frame(mask(crop(streams, circle), circle), xy = TRUE, na.rm = TRUE)
    cells <- cells[cells[, 3] > 0, ]
    if (nrow(cells) == 0) return(xy)
    i <- which.min((cells$x - xy[1])^2 + (cells$y - xy[2])^2)
    c(cells$x[i], cells$y[i])
  }
  coords[1, ] <- snap(coords[1, ])
  coords[nrow(coords), ] <- snap(coords[nrow(coords), ])

  line <- st_sfc(st_linestring(coords), crs = crs(dem))
  len <- as.numeric(st_length(line))
  frac <- seq(0, 1, length.out = ceiling(len / res(dem)[1]) + 1)
  pts <- st_cast(st_line_sample(line, sample = frac), 'POINT')
  z <- extract(dem, st_coordinates(pts))[, 1]
  profile <- data.frame(dist_m = frac * len, z_m = z)

  list(stats = data.frame(length_m = len, drop_m = abs(z[1] - z[length(z)]),
                          slope = abs(z[1] - z[length(z)]) / len,
                          downstream_end = if (z[length(z)] <= z[1]) 'end' else 'start'),
       profile = profile, line = line)
}
