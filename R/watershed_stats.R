# watershed_stats(): area, relief ratio, mean slope.
#
# "Watershed slope" as elevation change / watershed length is a relief ratio:
#   (max elevation - outlet elevation) / longest flow path
# Mean terrain slope is a different thing. It depends on grid size (1 m lidar
# sees every boulder and reads steeper), so it's computed on a 10 m grid.

library(sf)
library(terra)
library(whitebox)

watershed_stats <- function(ws, outlet, dem, files) {

  f <- function(x) file.path(files$work_dir, x)
  dem_ws <- mask(crop(dem, vect(ws)), vect(ws))

  z_outlet <- extract(dem, vect(outlet))[1, 2]
  z_max <- global(dem_ws, 'max', na.rm = TRUE)[1, 1]

  # longest flow path (whitebox needs the watershed as a raster on the DEM grid)
  writeRaster(rasterize(vect(ws), rast(files$breached), field = 1), f('basin.tif'), overwrite = TRUE)
  wbt_longest_flowpath(files$breached, f('basin.tif'), f('lfp.shp'))
  lfp <- st_set_crs(st_read(f('lfp.shp'), quiet = TRUE), crs(dem))
  lfp <- lfp[which.max(st_length(lfp)), ]
  lfp_m <- as.numeric(st_length(lfp))

  slope <- terrain(aggregate(dem_ws, 10 / res(dem)[1]), 'slope', unit = 'degrees')

  list(stats = data.frame(area_ha = as.numeric(st_area(ws)) / 1e4,
                          z_outlet_m = z_outlet, z_max_m = z_max,
                          relief_m = z_max - z_outlet, lfp_m = lfp_m,
                          relief_ratio = (z_max - z_outlet) / lfp_m,
                          mean_slope_deg = global(slope, 'mean', na.rm = TRUE)[1, 1]),
       lfp = lfp)
}
