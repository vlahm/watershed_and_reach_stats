# get_dem(): fetch USGS 3DEP 1 m lidar around a point, straight from The
# National Map's S3 bucket. GDAL's /vsicurl reads only the window we need,
# so a 3 km square arrives in ~10 s instead of downloading whole 300 MB tiles.
# Output is in NAD83 / UTM (meters), which is what whitebox and slope math want.

library(sf)
library(terra)
library(jsonlite)

get_dem <- function(pt, buffer_m = 1500) {

  # UTM zone from longitude
  ll <- st_coordinates(st_transform(st_geometry(pt), 4326))[1, ]
  epsg <- 26900 + floor((ll[1] + 180) / 6) + 1
  bb_utm <- st_bbox(st_buffer(st_transform(st_geometry(pt), epsg), buffer_m))
  bb_ll <- st_bbox(st_transform(st_as_sfc(bb_utm), 4326))

  # ask The National Map which 1 m tiles cover the box
  url <- paste0('https://tnmaccess.nationalmap.gov/api/v1/products',
                '?datasets=', URLencode('Digital Elevation Model (DEM) 1 meter', reserved = TRUE),
                '&bbox=', paste(bb_ll[c('xmin', 'ymin', 'xmax', 'ymax')], collapse = ','),
                '&outputFormat=JSON&max=100')
  tiles <- fromJSON(url)$items$downloadURL
  if (is.null(tiles)) stop('no 1 m lidar here')

  # several lidar projects can overlap; keep tiles in our UTM zone from the
  # project with the most tiles
  tiles <- tiles[grepl(paste0('_', epsg - 26900, '_x'), tiles)]
  project <- basename(dirname(dirname(tiles)))
  tiles <- tiles[project == names(which.max(table(project)))]
  message('DEM: ', length(tiles), ' tile(s) from ', unique(project[project == names(which.max(table(project)))]))

  e <- ext(bb_utm[c('xmin', 'xmax', 'ymin', 'ymax')])
  pieces <- lapply(tiles, function(u) {
    r <- rast(paste0('/vsicurl/', u))
    if (is.null(intersect(ext(r), e))) NULL else crop(r, e)
  })
  pieces <- Filter(Negate(is.null), pieces)
  dem <- if (length(pieces) > 1) do.call(merge, pieces) else pieces[[1]]
  names(dem) <- 'elev_m'
  dem
}
