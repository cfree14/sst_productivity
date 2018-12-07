
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)
library(RColorBrewer)
library(quantreg)
library(mblm) # For Thiel-Sen slope

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.csv", sep="/"), as.is=T)


# Stats for MS
################################################################################

table(data$class)


# Plot data
################################################################################

# Setup figure
# figname <- "Science_Fig2_theta_drivers.png"
# png(paste(plotdir, figname, sep="/"), width=6.5, height=2.5, units="in", res=600)
figname <- "Science_Fig2_theta_drivers.pdf"
pdf(paste(plotdir, figname, sep="/"), width=6.5, height=2.5)
par(mfrow=c(1,3), mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.2,0.8,0))

# Colors
##########################################

# Add colors
summary(data$betaT)
betaT_breaks <- round(seq(-0.6,0.6,0.1),1)
data$betaT_bin <- cut(data$betaT, breaks=betaT_breaks)
colors <- colorpal(brewer.pal(11, "RdBu"), nlevels(data$betaT_bin))
colors_tr <- rgb(t(col2rgb(colors))/255, alpha=0.7)
data$betaT_bin_color <- colors[data$betaT_bin]
data$col <- revalue(data$betaT_inf, c("positive"="blue", "negative"="red", "none"="grey60"))
data$sst_c_trend <- data$sst_c_trend * 10

# F/FMSY*warming interaction
##########################################

# Plot F/FMSY * SST trend interaction
plot(ffmsy_avg ~ sst_c_trend, data, bty="n", las=1,
     cex.axis=1.07, cex.lab=1.1,
     xlim=c(-0.2, 0.8), ylim=c(0,5),
     xlab="Temperature trend (°C/10yr)", ylab=expression("Mean F/F"["MSY"]), 
     pch=21, bg=data$betaT_bin_color, cex=1.2)
lines(x=c(-0.2, 0.8), y=c(1, 1), lty=3, lwd=1.2)
lines(x=c(0,0), y=c(0,5), lty=3, lwd=1.2)
text(x=(-0.2)-(0.8--0.2)*0.28, y=5, labels="A", font=2, cex=1.4, xpd=NA)

# Add sample size
n1 <- nrow(filter(data, !is.na(ffmsy_avg) & !is.na(sst_c_trend)))
text(labels=paste0("n=", n1), x=-0.2+(0.8--0.2)*1.05, y=0, pos=2, cex=1, xpd=NA)

# Maximum age
##########################################

# Plot data
plot(betaT ~ tmax_yr, data, bty="n", las=1, xpd=NA,
     cex.axis=1.07, cex.lab=1.1,
     xaxt="n", xlim=c(0, 100), ylim=c(-1,1), 
     xlab="Maximum age (yr)", ylab="",
     col=col, pch=1, cex=1.2)
axis(1, at=seq(0, 100, 20), cex.axis=1.07)
lines(x=c(0, 100), y=c(0,0), lty=3, lwd=1.2)
text(x=(0)-(100-0)*0.325, y=1, labels="B", font=2, cex=1.4, xpd=NA)

# Fit and plot quantile regression
qrfit <- rq(betaT ~ tmax_yr, data, tau=0.5)
qrfit_lo <- rq(betaT ~ tmax_yr, data,  tau=0.025)
qrfit_hi <- rq(betaT ~ tmax_yr, data,  tau=0.975)
curve(coef(qrfit)[1]+coef(qrfit)[2]*x, from=0, to=100, n=50, add=T, lwd=1.2)
curve(coef(qrfit_lo)[1]+coef(qrfit_lo)[2]*x, from=0, to=100, n=50, add=T, lwd=1.2, lty=2)
curve(coef(qrfit_hi)[1]+coef(qrfit_hi)[2]*x, from=0, to=100, n=50, add=T, lwd=1.2, lty=2)

# Add axis label
mtext("Temperature influence", side=2, adj=0.5, cex=0.8, line=2.8)

# Add sample size
n2 <- nrow(filter(data, !is.na(tmax_yr)))
text(labels=paste0("n=", n2), x=-0+(100-0)*1.05, y=-1, pos=2, cex=1, xpd=NA)

# Position in thermal niche
##########################################

# Species
spp <- c("Gadus morhua", "Clupea harengus")
cols <- brewer.pal(4, "Set1")[3:4]

# Setup empty
plot(1:10, 1:10, type="n", bty="n", las=1,
     cex.axis=1.07, cex.lab=1.1,
     xlim=c(0, 20), ylim=c(-1, 1), 
     xlab="Mean temperature (°C)", ylab="")
text(x=(0)-(20-0)*0.325, y=1, labels="C", font=2, cex=1.4, xpd=NA)

# Loop through species
for(i in 1:length(spp)){
  
  # Subset data
  sci_name <- spp[i]
  sdata <- subset(data, species==sci_name)
  
  # Add points
  points(betaT ~ sst_c_avg, sdata, pch=16, col=cols[i], cex=1.2)
  
  # Fit and plot Thiel-Sen slope
  tsfit <- mblm(betaT ~ sst_c_avg, sdata, repeated=F)
  # pvalue <- roundf(summary(tsfit)$coefficients[2,4],2)
  # lty <- ifelse(pvalue<0.1, 1, 2)
  curve(coef(tsfit)[1] + coef(tsfit)[2]*x, from=0, to=20, n=100, add=T, lty=1, col=cols[i], lwd=1.1)

}

# Add axis label
mtext("Temperature influence", side=2, adj=0.5, cex=0.8, line=2.8)

# Add legend
legend("topright", bty="n", col=cols, pch=16, lty=1, 
       legend=c("Atlantic cod (n=12)", "Atlantic herring (n=10)"))

# Add 0 influence line
lines(x=c(0, 20), y=c(0,0), lty=3, lwd=1.2)


# Off
dev.off()
graphics.off()

