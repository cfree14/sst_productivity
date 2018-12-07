
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(maps)
library(mapdata)
library(rgdal)
library(plyr)
library(dplyr)
library(freeR)
library(reshape2)
library(RColorBrewer)
library(fields) # colorbar.plot()

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/tables"
faodir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/fao_stat_areas/gis_data"

# Read data
stocks <- read.csv(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.csv", sep="/"), as.is=T)

# Read FAO data
fao_areas <- readOGR(dsn=faodir, laye="FAO_AREAS_NOCOASTLINE", verbose=F)
fao_areas <- subset(fao_areas, F_LEVEL=="MAJOR")
# plot(fao_areas)


# Plot data
################################################################################

# Setup figure
figname <- "SFig8_theta_map.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3, units="in", res=600)
layout(matrix(c(1,2), ncol=2, byrow=T), widths=c(0.9,0.1))
par(mar=c(0.1,0.1,0.1,0.1), xpd=NA)

# Parameters
# xlim <- c(-160, 170)
# ylim <- c(-50, 80)
xlim <- c(-180, 180)
ylim <- c(-90, 90)

# A. Theta map
####################################

# Randomize order for plotting
set.seed(1)
stocks <- stocks[sample(1:nrow(stocks)),]

# Add colors
summary(stocks$betaT)
betaT_breaks <- seq(-0.6,0.6,0.1)
stocks$betaT_bin <- cut(stocks$betaT, breaks=betaT_breaks)
colors <- colorpal(brewer.pal(11, "RdBu"), nlevels(stocks$betaT_bin))
colors_tr <- rgb(t(col2rgb(colors))/255, alpha=0.7)
stocks$betaT_bin_color <- colors_tr[stocks$betaT_bin]

# Plot FAO stat areas
plot(fao_areas, border="grey70", lty=3, xlim=xlim, ylim=ylim)

# Plot world countries
map("world", col="grey85", fill=T, border="white", lwd=0.3, 
    xlim=xlim, ylim=ylim, add=T)

# Add betaT values
stocks_look <- select(stocks, stockid, long_dd, lat_dd)
stocks_on_land <- c("SAILWATL", "SKJWATL", "SAILEATL", "SKJEATL")
stocks_plot <- subset(stocks, !stockid %in% stocks_on_land)
points(jitter(stocks_plot$long_dd,120), jitter(stocks_plot$lat_dd,120), pch=21, cex=1.6,
       bg=stocks$betaT_bin_color, col="grey30")

# Legend
plot(1:10, 1:10, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="")
colorbar.plot(x=0, y=1, adj.x=0, adj.y=0, horiz=F, col=colors,
              strip=seq(-1.5,1.5,length.out=nlevels(stocks$betaT_bin)), 
                        strip.width=0.3, strip.length=2, xpd=NA)
text(x=2, y=c(1.1, 2.95, 4.85), pos=4, labels=c("-0.6", "0.0", "0.6"), cex=0.8)
text(x=-0.5, y=5.5, pos=4, labels="SST\ninfluence", font=2, cex=0.8)

# Off
dev.off()
graphics.off()

