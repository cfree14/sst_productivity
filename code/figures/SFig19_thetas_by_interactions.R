
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(tidyverse)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)


# Explore data
################################################################################

# Add BetaT colors
summary(data$betaT)
betaT_breaks <- seq(-1.25,1.25,0.25)
data$betaT_bin <- cut(data$betaT, breaks=betaT_breaks)
colors <- colorpal(brewer.pal(11, "RdBu"), nlevels(data$betaT_bin))
colors_tr <- rgb(t(col2rgb(colors))/255, alpha=0.7)
data$betaT_bin_color <- colors[data$betaT_bin]

# Plot F/FMSY * SST trend interaction
plot(ffmsy_avg ~ sst_c_trend, data, bty="n", las=1,
     xlim=c(-0.02, 0.08), ylim=c(0,5),
     xlab="SST trend (°C/1yr)", ylab=expression("F/F"["MSY"]*" mean"), 
     pch=21, bg=data$betaT_bin_color, cex=1.2)
lines(x=c(-0.02, 0.08), y=c(1, 1), lty=3)
lines(x=c(0,0), y=c(0,5), lty=3)

# Plot F/FMSY * SST mean interaction
plot(ffmsy_avg ~ sst_c_avg, data, bty="n", las=1,
     ylim=c(0,5),
     xlab="SST mean (°C)", ylab=expression("F/F"["MSY"]*" mean"), 
     pch=21, bg=data$betaT_bin_color, cex=1.2)
lines(x=c(0, 25), y=c(1, 1), lty=3)

# Plot F/FMSY * maximum age interaction
plot(ffmsy_avg ~ tmax_yr, data, bty="n", las=1,
     xlim=c(0, 100), ylim=c(0,5),
     xlab="Maximum age (yr)", ylab=expression("F/F"["MSY"]*" mean"), 
     pch=21, bg=data$betaT_bin_color, cex=1.2)
lines(x=c(0, 100), y=c(1, 1), lty=3)


# Plot SST trend * SST mean interaction
plot(sst_c_trend ~ sst_c_avg, data, bty="n", las=1,
     xlab="SST mean (°C)", ylab="SST trend (°C/yr)", 
     pch=21, bg=data$betaT_bin_color, cex=1.2)
lines(x=c(0, 25), y=c(0, 0), lty=3)


# Plot data
################################################################################

# Plot F/FMSY * SST trend interaction
plot(ffmsy_avg ~ sst_c_trend, data, bty="n", las=1,
     xlim=c(-0.02, 0.08), ylim=c(0,5),
     xlab="SST trend (°C/1yr)", ylab=expression("F/F"["MSY"]*" mean"), 
     pch=21, bg=data$betaT_bin_color, cex=1.2)
lines(x=c(-0.02, 0.08), y=c(1, 1), lty=3)
lines(x=c(0,0), y=c(0,5), lty=3)

# Build new dataset
data1 <- data %>% 
  filter(!is.na(ffmsy_avg)) %>% 
  select(stockid, ffmsy_avg, sst_c_trend, betaT) %>% 
  mutate(ffmsy_bin=cut(ffmsy_avg, breaks=seq(0,5,0.5)),
         trend_bin=cut(sst_c_trend, breaks=seq(-0.02,0.08,0.01)))

# Reshape to grid
data2 <- reshape2::dcast(data1, ffmsy_bin ~ trend_bin, value.var="betaT", fun.aggregate=mean)
data3 <- t(as.matrix(data2[,2:ncol(data2)]))

# Setup figure
figname <- "SFig19_thetas_overfishing_warming_interaction.png"
png(paste(plotdir, figname, sep="/"), width=5, height=5, units="in", res=600)
par(mar=c(3.5,3.5,0.5,0.8), mgp=c(2.2,0.7,0))

# Plot interaction
image(y=seq(0,5,0.5), x=seq(-0.02,0.08,0.01), z=data3, las=1,
      zlim=c(-1.25,1.25), col=colors,
      xlab="SST trend (°C/1yr)", ylab=expression("F/F"["MSY"]*" mean"))
lines(x=c(-0.02, 0.08), y=c(1, 1), lty=3)
lines(x=c(0,0), y=c(0,5), lty=3)
points(data1$sst_c_trend, data1$ffmsy_avg, col="grey40")


# Off
dev.off()
graphics.off()


