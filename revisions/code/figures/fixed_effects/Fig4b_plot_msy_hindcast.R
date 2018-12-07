

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

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"


# Plot figure
################################################################################

# Setup figure
figname <- "FE_Fig4_msy_hindcast_fixed_random.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3.5, units="in", res=600)
par(mfrow=c(1,2), mar=c(3.5,2,1.5,0.5), mgp=c(2.2,0.8,0), oma=c(0,2,0,0))

# A. Plot MSY overall - fixed effects
##########################################

# Read hindcast data
load(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_fixed_hindcast_10000traj_unorm.Rdata", sep="/"))

# Stats
yrs1 <- 1930:1939
yrs2 <- 2001:2010
msy1 <- mean(msy_trajs_df$msy[msy_trajs_df$year%in%yrs1])
msy2 <- mean(msy_trajs_df$msy[msy_trajs_df$year%in%yrs2])
perc_diff <- (msy2-msy1)/msy2*100

# Subset data to use
hindcast_yrs <- c(1930:2010)
msy_df <- subset(msy_trajs_df, year%in%hindcast_yrs)
msy_pos_df <- subset(msy_trajs_pos_df, year%in%hindcast_yrs)
msy_neg_df <- subset(msy_trajs_neg_df, year%in%hindcast_yrs)
msy_non_df <- subset(msy_trajs_non_df, year%in%hindcast_yrs)
yr1 <- min(hindcast_yrs)
yr2 <- max(hindcast_yrs)

# Plot data
min(msy_df$msy_lo)
max(msy_df$msy_hi)
plot(msy ~ year, msy_df, type="n", bty="n", cex.lab=1.2, cex.axis=1.05, las=1,
     xaxt="n", xlab="", ylab="MSY (millions of mt)", ylim=c(30,50), xpd=NA, main="Fixed effects")
axis(1, at=seq(yr1, yr2, 10), las=2)
polygon(x=c(msy_df$year, rev(msy_df$year)),
        y=c(msy_df$msy_lo, rev(msy_df$msy_hi)),
        col="grey80", border=F)
lines(x=msy_df$year, y=msy_df$msy, col="black")
lines(x=c(yr1,yr2), y=rep(msy_avg,2), lty=3, lwd=1.2)

# Add text
ptext <- paste0(roundf(abs(perc_diff),1), "% decline")
text(x=2010, y=50, labels=ptext, pos=2, cex=0.8)


# B. Plot MSY overall - random effects
##########################################

# Read hindcast data
load(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_hindcast_10000traj_unorm.Rdata", sep="/"))

# Stats
yrs1 <- 1930:1939
yrs2 <- 2001:2010
msy1 <- mean(msy_trajs_df$msy[msy_trajs_df$year%in%yrs1])
msy2 <- mean(msy_trajs_df$msy[msy_trajs_df$year%in%yrs2])
perc_diff <- (msy2-msy1)/msy2*100

# Subset data to use
hindcast_yrs <- c(1930:2010)
msy_df <- subset(msy_trajs_df, year%in%hindcast_yrs)
msy_pos_df <- subset(msy_trajs_pos_df, year%in%hindcast_yrs)
msy_neg_df <- subset(msy_trajs_neg_df, year%in%hindcast_yrs)
msy_non_df <- subset(msy_trajs_non_df, year%in%hindcast_yrs)
yr1 <- min(hindcast_yrs)
yr2 <- max(hindcast_yrs)

# Plot data
min(msy_df$msy_lo)
max(msy_df$msy_hi)
plot(msy ~ year, msy_df, type="n", bty="n", cex.lab=1.2, cex.axis=1.05, las=1,
     xaxt="n", xlab="", ylab="", ylim=c(30,50), main="Random effects")
axis(1, at=seq(yr1, yr2, 10), las=2)
polygon(x=c(msy_df$year, rev(msy_df$year)),
        y=c(msy_df$msy_lo, rev(msy_df$msy_hi)),
        col="grey80", border=F)
lines(x=msy_df$year, y=msy_df$msy, col="black")
lines(x=c(yr1,yr2), y=rep(msy_avg,2), lty=3, lwd=1.2)

# Add text
ptext <- paste0(roundf(abs(perc_diff),1), "% decline")
text(x=2010, y=50, labels=ptext, pos=2, cex=0.8)

# Off
dev.off()
graphics.off()


