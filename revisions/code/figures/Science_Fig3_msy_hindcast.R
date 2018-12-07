

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

# Read hindcast data
load(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme_msy_hindcast_10000traj_mvnorm.Rdata", sep="/"))


# Stats for MS
################################################################################

# Stats
yrs1 <- 1930:1939
yrs2 <- 2001:2010
msy1 <- mean(msy_trajs_df$msy[msy_trajs_df$year%in%yrs1]) / 1E6
msy2 <- mean(msy_trajs_df$msy[msy_trajs_df$year%in%yrs2]) / 1E6
perc_diff <- (msy2-msy1)/msy2*100

# What happens if we only hindcasted 1970-2010?
yrs1b <- 1970:1979
msy1b <- mean(msy_trajs_df$msy[msy_trajs_df$year%in%yrs1b])
perc_diff2 <- (msy2-msy1b)/msy2*100

# How much SST change have we observed since 1970 compared to 1930?
sst$sst_c[sst$year==2015] - sst$sst_c[sst$year==1970] # 0.52
sst$sst_c[sst$year==2015] - sst$sst_c[sst$year==1930] # 0.59

# Confidence intervals for percent MSY change
stats <- data.frame(id=1:1000, msy1=NA, msy2=NA, pdiff=NA)
for(i in 1:nrow(stats)){
  msy_ts <- msy_trajs_sum[,i]
  msy1 <- mean(msy_ts[which(names(msy_ts)%in%1930:1939)])
  msy2 <- mean(msy_ts[which(names(msy_ts)%in%2001:2010)])
  pdiff <- (msy2-msy1)/msy2*100
  stats$msy1[i] <- msy1
  stats$msy2[i] <- msy2
  stats$pdiff[i] <- pdiff
}
stats$msy_diff <- stats$msy2-stats$msy1
summary(stats$pdiff)
quantile(stats$pdiff, probs=c(0.025,0.5, 0.975))

#  Stats for MS
mean(stats$pdiff)
mean(stats$msy1) / 1E6
mean(stats$msy2) / 1E6
mean(stats$msy_diff) / 1E6

# Proportion negative
sum(stats$pdiff<0) / nrow(stats)
sum(stats$pdiff>0) / nrow(stats)



# Plot data
################################################################################

# Subset data to use
hindcast_yrs <- c(1930:2010)
msy_df <- subset(msy_trajs_df, year%in%hindcast_yrs)
msy_pos_df <- subset(msy_trajs_pos_df, year%in%hindcast_yrs)
msy_neg_df <- subset(msy_trajs_neg_df, year%in%hindcast_yrs)
msy_non_df <- subset(msy_trajs_non_df, year%in%hindcast_yrs)
yr1 <- min(hindcast_yrs)
yr2 <- max(hindcast_yrs)

# Colors
display.brewer.pal(11, "RdBu")
colors <- brewer.pal(11, "RdBu")

# Sample size
n_inf <- table(stocks$betaT_inf)

# Setup figure
# figname <- "Science_Fig3_msy_hindcast.png"
# png(paste(plotdir, figname, sep="/"), width=6.5, height=5, units="in", res=600)
figname <- "Science_Fig3_msy_hindcast.pdf"
pdf(paste(plotdir, figname, sep="/"), width=6.5, height=5)
layout(matrix(c(1,3,
                1,4,
                2,5), byrow=T, ncol=2), widths=c(0.55, 0.45)) 
par(mar=c(1,4,0.5,0.5), mgp=c(2.4,0.7,0), oma=c(3,0,0,0))

# A. Plot MSY overall
##########################################

# Purple options
# purple.col <- brewer.pal(4, "Set1")[4]
purple.col <- "mediumorchid4"
purple.col <- "purple4"

# Plot data
min(msy_df$msy_lo)
max(msy_df$msy_hi)
plot(msy ~ year, msy_df, type="n", bty="n", cex.lab=1.2, cex.axis=1.05, las=1,
     xaxt="n", xlab="", ylab="MSY (millions of mt)", ylim=c(30,40))
axis(1, at=seq(yr1, yr2, 10))
polygon(x=c(msy_df$year, rev(msy_df$year)),
        y=c(msy_df$msy_lo, rev(msy_df$msy_hi)),
        col="grey80", border=F)
lines(x=msy_df$year, y=msy_df$msy, col="black")
mtext("A", side=3, adj=0.05, line=-2, cex=0.8, font=2)
text(x=1934, y=39.64, pos=4, labels=paste(sum(n_inf), "populations"), cex=1.05)
lines(x=c(yr1,yr2), y=rep(msy_avg,2), lty=3, lwd=1.2)

# E. Plot SST data
##########################################

# Plot SST data
plot(sst_c_sd ~ year, sst, type="n", bty="n", las=1, cex.lab=1.2, cex.axis=1.05,
     xaxt="n", xlim=c(1850,2020), ylim=c(-0.6, 0.6),
     xlab="", ylab="SST anomaly (°C)", col="grey30")
rect(xleft=1930, xright=2010, ybottom=-0.3, ytop=0.35, col="grey90", border=F)
text(x=1928, y=0.28, pos=4, labels="Hindcast window", cex=0.9, col="grey50")
lines(sst$year, sst$sst_c_sd, col="grey30")
axis(1, at=seq(1850,2020,10), labels=T, las=2)
mtext("E", side=3, adj=0.05, line=-1.5, cex=0.8, font=2)
lines(x=c(1850,2020), y=c(0,0), lty=3, col="grey50")
text(x=1880, y=-0.01, pos=3, labels="1971-2000 average", cex=0.9, col="grey50")

# B. Plot MSY significant positive
##########################################

# New par
par(mgp=c(2.2,0.7,0))

# Plot data
min(msy_pos_df$msy_lo)
max(msy_pos_df$msy_hi)
plot(msy ~ year, msy_pos_df, type="n", bty="n", las=1, cex.axis=1.05,
     xaxt="n", xlab="", ylab="", ylim=c(0.1,0.4), yaxt="n")
axis(1, at=seq(yr1, yr2, 10), labels=F)
axis(2, at=seq(0.1, 0.4, 0.1), cex.axis=1.05, las=1)
polygon(x=c(msy_pos_df$year, rev(msy_pos_df$year)),
        y=c(msy_pos_df$msy_lo, rev(msy_pos_df$msy_hi)),
        col=colors[9], border=F)
lines(x=msy_pos_df$year, y=msy_pos_df$msy, col=colors[11])
lines(x=c(yr1,yr2), y=rep(msy_avg_pos,2), lty=3, lwd=1.2)
mtext("B", side=3, adj=0.05, line=-1.5, cex=0.8, font=2)
text(x=1934, y=0.4-(0.4-0.1)*0.07, pos=4, labels=paste(n_inf["positive"], "populations"), cex=1.05)

# C. Plot MSY significant negative
##########################################

# Plot data
min(msy_neg_df$msy_lo)
max(msy_neg_df$msy_hi)
plot(msy ~ year, msy_neg_df, type="n", bty="n", las=1, cex.axis=1.05, 
     xaxt="n", xlab="", ylab="MSY (millions of mt)", ylim=c(1,6), yaxt="n",cex.lab=1.2)
axis(1, at=seq(yr1, yr2, 10), labels=F)
axis(2, at=1:6, cex.axis=1.05, las=1)
polygon(x=c(msy_neg_df$year, rev(msy_neg_df$year)),
        y=c(msy_neg_df$msy_lo, rev(msy_neg_df$msy_hi)),
        col=colors[3], border=F)
lines(x=msy_neg_df$year, y=msy_neg_df$msy, col=colors[1])
lines(x=c(yr1,yr2), y=rep(msy_avg_neg,2), lty=3, lwd=1.2)
mtext("C", side=3, adj=0.05, line=-1.5, cex=0.8, font=2)
text(x=1934, y=6-(6-1)*0.07, pos=4, labels=paste(n_inf["negative"], "populations"), cex=1.05)

# D. Plot MSY non-significant
##########################################

# Plot data
min(msy_non_df$msy_lo)
max(msy_non_df$msy_hi)
plot(msy ~ year, msy_non_df, type="n", bty="n", las=1, cex.axis=1.05,
     xaxt="n", xlab="", ylab="", ylim=c(28,36))
axis(1, at=seq(yr1, yr2, 10), las=2, cex.axis=1.05)
polygon(x=c(msy_non_df$year, rev(msy_non_df$year)),
        y=c(msy_non_df$msy_lo, rev(msy_non_df$msy_hi)),
        col="grey80", border=F)
lines(x=msy_non_df$year, y=msy_non_df$msy)
lines(x=c(yr1,yr2), y=rep(msy_avg_non,2), lty=3, lwd=1.2)
mtext("D", side=3, adj=0.05, line=-1.5, cex=0.8, font=2)
text(x=1934, y=36-(36-28)*0.07, pos=4, labels=paste(n_inf["none"], "populations"), cex=1.05)

# Off
dev.off()
graphics.off()


