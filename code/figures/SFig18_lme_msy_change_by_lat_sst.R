
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/tables"

# Read data
data <- read.csv(paste(tabledir, "STable6_lme_msy_hindcast.csv", sep="/"), as.is=T)


# Plot data
################################################################################

# Setup figure
figname <- "SFig18_lme_msy_change_by_lat_sst.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3, units="in", res=600)
par(mfrow=c(1,3), mar=c(4,1,1,1), mgp=c(2.3, 0.5, 0), oma=c(0,4,0,0))

# % MSY change ~ latitute
plot(pdiff ~ abs(lat_dd), data, bty="n", las=1, xpd=NA, pch=21, bg="grey80", cex=1.3, cex.lab=1.2,
     xlim=c(0,80), ylim=c(-80,40), xlab="Latitude (° absolute)", ylab="% difference in MSY\n(1930-39 to 2000-10)")
lines(x=c(0,80), y=c(0,0), lty=3)

# % MSY change ~ SST mean
plot(pdiff ~ sst_c_mean, data, bty="n", las=1, pch=21, bg="grey80", cex=1.3,
     xlim=c(0,30), ylim=c(-80,40), xlab="SST mean (°C)", ylab="", cex.lab=1.2)
lines(x=c(0,30), y=c(0,0), lty=3)

# % MSY change ~ SST trend
plot(pdiff ~ sst_c_slope_10yr, data, bty="n", las=1, pch=21, bg="grey80", cex=1.3,
     xlab="SST trend (°C/10yr)", ylab="", ylim=c(-80,40), cex.lab=1.2)
lines(x=c(-0.1,0.15), y=c(0,0), lty=3)

# Off
dev.off()
graphics.off()








