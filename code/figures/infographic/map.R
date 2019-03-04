
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(ggplot2)
library(countrycode)
library(rgdal)
library(sf)
library(raster)
library(doParallel)
library(maptools)
library(mapdata)

# Directories
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures/infographic"

# Projections
moll84 <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Get world map
data(wrld_simpl)
# world <- map_data("world2Hires")
world <- spTransform(wrld_simpl, wgs84)

# Plot maps
png(file.path(plotdir, "background_map_hiRes.png"), width=16, height=8, units="in", res=600)
plot(world, col="grey80", border="white", lwd=0.6)
dev.off()

