

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)
library(reshape2)
library(RColorBrewer)
library(mblm) # For Thiel-Sen slope

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"

# Read hindcast data
load(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme_msy_hindcast_10000traj_mvnorm.Rdata", sep="/"))
msy_trajs_df$year <- as.numeric(msy_trajs_df$year)

# Build data
################################################################################

# Stats
yrs2 <- 2001:2010
msy2 <- mean(msy_trajs_df$msy[msy_trajs_df$year%in%yrs2])

# Setup dataframe
stats <- data.frame(year=as.numeric(1850:1990), slope=NA, msy_avg=NA, pdiff=NA)
stats$year <- as.numeric(stats$year)

# Loop through years and build data
for(i in 1:nrow(stats)){
  
  # Subset data
  yr1 <- stats$year[i]
  sdata <- subset(msy_trajs_df, year%in%yr1:2010)
  
  # Fit Thiel-Sen slope
  tsfit <- mblm(msy ~ year, sdata, repeated=F)
  stats$slope[i] <- coef(tsfit)[2]
  
  # Percent change in 
  msy_avg_10yr <- mean(sdata$msy[1:10])
  stats$msy_avg[i] <- msy_avg_10yr
  stats$pdiff[i] <- (msy2-msy_avg_10yr)/msy2*100

}

# Plot data
################################################################################

# Hindcast year
hind_yr <- 1930
years <- 1930:2010
hdata <- subset(msy_trajs_df, year%in%years)
tsfit <- mblm(msy ~ year, hdata, repeated=F)

# Setup figure
figname <- "SFig17_hindcast_sensitivity.png"
png(paste(plotdir, figname, sep="/"), width=4, height=7, units="in", res=600)
par(mfrow=c(4,1), mar=c(0.5, 6, 0.5, 0.5), oma=c(3,0,0,0), mgp=c(2.8,0.8,0))

# A. Plot SST
plot(sst_c_sd ~ year, sst, type="n", bty="n", las=1,
     xaxt="n", xlim=c(1850,2020), ylim=c(-0.6, 0.6), cex.lab=1.2, cex.axis=1.1,
     xlab="", ylab="SST anomaly (°C)", col="grey30")
rect(xleft=1930, xright=2010, ybottom=-0.3, ytop=0.35, col="grey90", border=F)
text(x=1928, y=0.28, pos=4, labels="Hindcast window", cex=0.9, col="grey50")
lines(sst$year, sst$sst_c_sd, col="grey30")
axis(1, at=seq(1850,2020,10), labels=F, las=2)
lines(x=c(1850,2020), y=c(0,0), lty=3, col="grey50")
text(x=1875, y=-0.01, pos=3, labels="1971-2000 average", cex=1, col="grey50")
mtext("A", side=3, adj=0.05, line=-1.5, cex=0.8, font=2)

# B. Plot hindcast
plot(msy ~ year, msy_trajs_df, type="n", bty="n", cex.lab=1.2, cex.axis=1.1, las=1,
     xaxt="n", xlim=c(1850,2020), xlab="", ylab="MSY (millions of mt)", ylim=c(30,40))
axis(1, at=seq(1850, 2020, 10), las=2, labels=F)
polygon(x=c(msy_trajs_df$year, rev(msy_trajs_df$year)),
        y=c(msy_trajs_df$msy_lo, rev(msy_trajs_df$msy_hi)),
        col="grey80", border=F)
lines(x=msy_trajs_df$year, y=msy_trajs_df$msy, col="grey30")
lines(x=c(1850,2020), y=rep(msy_avg,2), lty=3, lwd=1.2)
mtext("B", side=3, adj=0.05, line=-1.5, cex=0.8, font=2)
curve(coef(tsfit)[1]+coef(tsfit)[2]*x, from=1930, to=2010, n=101, add=T, lwd=1.4)

# C. Plot Thiel-Sen slope
hind_yr_slope <- stats$slope[stats$year==hind_yr]*10
slope_label <- paste0(roundf(abs(hind_yr_slope), 2), " million mt\nper decade")
plot(slope*10 ~ year, stats, type="l", bty="n", las=1, xaxt="n", cex.lab=1.2, cex.axis=1.1,
     xlim=c(1850, 2020), ylim=c(-1,0), xlab="", ylab="Δ MSY (millions of mt)\nper decade")
axis(1, at=seq(1850,2020,10), las=2, labels=F)
points(x=hind_yr, y=hind_yr_slope, pch=16, cex=1.5)
text(x=hind_yr, y=hind_yr_slope*0.7, pos=2, labels=slope_label)
mtext("C", side=3, adj=0.05, line=-1.5, cex=0.8, font=2)

# D. Plot percent difference
hind_yr_pdiff <- stats$pdiff[stats$year==hind_yr]
pdiff_label <- paste0(roundf(abs(hind_yr_pdiff), 1), "% decline")
plot(pdiff ~ year, stats, type="l", bty="n", las=1, xaxt="n", cex.lab=1.2, cex.axis=1.1,
     xlim=c(1850, 2020), ylim=c(-12,0), xlab="", ylab="Percent difference in \nMSY relative to 2001-2010")
axis(1, at=seq(1850,2020,10), las=2, labels=T, cex.axis=1.1)
points(x=hind_yr, y=hind_yr_pdiff, pch=16, cex=1.5)
text(x=hind_yr, y=hind_yr_pdiff, pos=2, labels=pdiff_label)
mtext("D", side=3, adj=0.05, line=-1.5, cex=0.8, font=2)

# Off
dev.off()
graphics.off()

