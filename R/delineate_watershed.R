# delineate_watershed(): DEM -> flow network -> snapped outlet -> watershed
#
#   1. Condition the DEM so water can flow everywhere: fill single-cell pits,
#      then least-cost breaching (carves through road embankments where the
#      culverts are, instead of ponding behind them).
#   2. D8 flow direction and accumulation; cells draining more than
#      stream_thresh_ha count as "stream".
#   3. Snap the site to the nearest stream cell within search_m, delineate.
#
# stream_thresh_ha is the knob that matters: too low and the site snaps to a
# rill or side channel; too high and the surveyed stream isn't a stream.
# Always look at the map.
#
# whitebox works on files, so intermediate rasters go in work_dir.

library(sf)
library(terra)
library(whitebox)
wbt_verbose(FALSE)

delineate_watershed <- function(pt, dem, work_dir,
                                search_m = 50,             # how far to look for a stream
                                stream_thresh_ha = 0.5) {  # minimum drainage area that makes a stream

  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  f <- function(x) file.path(work_dir, x)
  cell_m <- res(dem)[1]
  pt <- st_transform(st_geometry(pt), crs(dem))

  # flow network
  writeRaster(dem, f('dem.tif'), overwrite = TRUE)
  wbt_fill_single_cell_pits(f('dem.tif'), f('dem_nopits.tif'))
  wbt_breach_depressions_least_cost(f('dem_nopits.tif'), f('dem_breached.tif'),
                                    dist = ceiling(100 / cell_m), fill = TRUE)
  wbt_d8_pointer(f('dem_breached.tif'), f('d8.tif'))
  wbt_d8_flow_accumulation(f('d8.tif'), f('facc.tif'), out_type = 'cells', pntr = TRUE)
  wbt_extract_streams(f('facc.tif'), f('streams.tif'),
                      threshold = stream_thresh_ha * 1e4 / cell_m^2)

  # snap the site to the nearest stream cell, then delineate
  st_write(st_sf(id = 1, geometry = pt), f('site.shp'), append = FALSE, quiet = TRUE)
  wbt_jenson_snap_pour_points(f('site.shp'), f('streams.tif'), f('outlet.shp'),
                              snap_dist = search_m)
  outlet <- st_geometry(st_read(f('outlet.shp'), quiet = TRUE))
  snap_m <- as.numeric(st_distance(pt, outlet))
  if (snap_m > search_m) stop('no stream within ', search_m, ' m; raise search_m or lower stream_thresh_ha')
  message('snapped ', round(snap_m, 1), ' m to a stream cell')

  wbt_watershed(f('d8.tif'), f('outlet.shp'), f('watershed.tif'))
  ws <- st_union(st_as_sf(as.polygons(ifel(rast(f('watershed.tif')) > 0, 1, NA))))
  ws <- st_sf(area_ha = as.numeric(st_area(ws)) / 1e4, geometry = ws)

  if (any(abs(st_bbox(ws) - st_bbox(dem)) < 2 * cell_m)) {
    warning('watershed touches the DEM edge: refetch with a larger buffer_m')
  }

  list(watershed = ws, outlet = outlet, site = pt, snap_m = snap_m,
       files = list(dem = f('dem.tif'), breached = f('dem_breached.tif'),
                    d8 = f('d8.tif'), streams = f('streams.tif'), work_dir = work_dir))
}

# quick-look map. ALWAYS look at it.
plot_delineation <- function(d, dem, main = '') {
  bb <- st_bbox(st_buffer(st_union(st_geometry(d$watershed), d$site), 100))
  e <- ext(bb[c('xmin', 'xmax', 'ymin', 'ymax')])
  dem_c <- crop(dem, e)
  hs <- shade(terrain(dem_c, 'slope', unit = 'radians'), terrain(dem_c, 'aspect', unit = 'radians'))
  plot(hs, col = grey(0:100 / 100), legend = FALSE, main = main)
  plot(ifel(crop(rast(d$files$streams), e) > 0, 1, NA), col = 'dodgerblue', add = TRUE, legend = FALSE)
  plot(st_geometry(d$watershed), border = 'red', lwd = 2, add = TRUE)
  plot(d$site, pch = 21, bg = 'orange', cex = 1.5, add = TRUE)
  plot(d$outlet, pch = 4, col = 'red', lwd = 2, cex = 1.5, add = TRUE)
  mtext('blue: DEM streams   red: watershed   orange: site   x: snapped outlet',
        side = 3, line = 0.3, cex = 0.8)
}
