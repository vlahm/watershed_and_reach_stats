
library(tidyverse)
library(sf)
library(terra)
library(mapview)

source('R/get_dem.R')
source('R/delineate_watershed.R')
source('R/watershed_stats.R')
source('R/reach_slope.R')

#arbitrary headwater point within HBEF
s <- tibble(lat = 43.92801851113091, lon = -71.79193358276926) %>%
    st_as_sf(coords = c(y = 'lon', x = 'lat'), crs = 4326)

work_dir <- '/tmp'

dem <- get_dem(s, buffer_m = 1500)

d <- delineate_watershed(s, dem, work_dir = work_dir,
                         search_m = 50, stream_thresh_ha = 5)

st <- watershed_stats(d$watershed, d$outlet, dem, d$files)

mapview(d$watershed, map.types = "OpenTopoMap")


## reach slope demo: last 50 m of longest flow path
L <- as.numeric(st_length(st$lfp))
reach <- lwgeom::st_linesubstring(st_geometry(st$lfp), 1 - 50 / L, 1)  # line ends at the outlet
rs <- reach_slope(reach, dem, streams = d$files$streams)
print(rs$stats, digits = 3)
plot(rs$profile, type = 'b', pch = 16, xlab = 'distance along reach (m)',
     ylab = 'elevation (m)', main = 'W3: 50 m reach, 1 m lidar')
