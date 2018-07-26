
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(mblm) # For Thiel-Sen slope
library(plyr)
library(dplyr)
library(freeR)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output/fixed_effects"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures/fixed_effects"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_pella_0.20_cobe_merged.csv", sep="/"), as.is=T)


# Plot data
################################################################################

# Fit and plot quantile regression
# x <- data$ffmsy_avg; y <- betaT; xlim <- c(0,5)
fit_plot_qr <- function(x, y, xlim){
  qrfit <- rq(y ~ x, tau=0.5)
  qrfit_lo <- rq(y ~ x, tau=0.025)
  qrfit_hi <- rq(y ~ x, tau=0.975)
  curve(coef(qrfit)[1]+coef(qrfit)[2]*x, from=xlim[1], to=xlim[2], n=50, add=T, lwd=1.2)
  curve(coef(qrfit_lo)[1]+coef(qrfit_lo)[2]*x, from=xlim[1], to=xlim[2],n=50, add=T, lwd=1.2, lty=2)
  curve(coef(qrfit_hi)[1]+coef(qrfit_hi)[2]*x, from=xlim[1], to=xlim[2], n=50, add=T, lwd=1.2, lty=2)
}


# Models
models <- c("fixed", "random")
ylims <- matrix(c(-4, 8,
                  -1, 1.5), ncol=2, byrow=T)
zlims <- matrix(c(-8,8,
                  -1.25,1.25), ncol=2, byrow=T)

# Setup figure
figname <- "Fig6_drivers_fixed_vs_random_effects.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3.8, units="in", res=600)
par(mfrow=c(2,4), mar=c(2.5,2.5,0.5,0.5), mgp=c(2,0.7,0), oma=c(1,2,0,0))

# Loop and plot
for(i in 1:length(models)){
  
  # Reset par
  par(mgp=c(2,0.7,0))
  
  # Get data
  model <- models[i]
  if(model=="random"){
    betaT <- data$betaT_r
    colors <- revalue(data$betaT_inf_r, c("positive"="blue", "negative"="red", "none"="grey60"))
    xlabs <- c(expression("F/F"["MSY"]*" mean"), "Maximum age (yr)", "SST trend (°C/yr)")
    title <- "Random effects"
  }else{
    betaT <- data$betaT_f
    colors <- revalue(data$betaT_inf_f, c("positive"="blue", "negative"="red", "none"="grey60"))
    xlabs <- rep("", 3)
    title <- "Fixed effects"
  }
  
  # Plot BetaT ~ overfishing
  plot(x=data$ffmsy_avg, y=betaT, bty="n", las=1, col=colors, cex.axis=0.8,
       xlim=c(0,5), ylim=ylims[i,], xlab=xlabs[1], ylab="", xpd=NA)
  lines(x=c(0, 5), y=c(0,0), lty=3)
  mtext(title, side=3, adj=0.1, line=-1.3, cex=0.7, font=2)
  fit_plot_qr(x=data$ffmsy_avg, y=betaT, xlim=c(0,5))
  
  # Plot BetaT ~ maximum age
  plot(x=data$tmax_yr, y=betaT, bty="n", las=1, col=colors, cex.axis=0.8,
       xlim=c(0,100), ylim=ylims[i,], xlab=xlabs[2], ylab="", xpd=NA)
  lines(x=c(0, 100), y=c(0,0), lty=3)
  fit_plot_qr(x=data$tmax_yr, y=betaT, xlim=c(0,100))
  
  # Plot BetaT ~ SST trend
  plot(x=data$sst_c_trend, y=betaT, bty="n", las=1, col=colors, cex.axis=0.8,
       xlim=c(-0.02, 0.08), ylim=ylims[i,], xlab=xlabs[3], ylab="", xpd=NA)
  lines(x=c(-0.02, 0.08), y=c(0,0), lty=3)
  fit_plot_qr(x=data$sst_c_trend, y=betaT, xlim=c(-0.02,0.08))
  
  # Plot BetaT ~ overfishing*warming interaction
  #####################################################

  # Build new dataset
  data1 <- data %>% 
    filter(!is.na(ffmsy_avg)) %>% 
    select(stockid, ffmsy_avg, sst_c_trend, betaT_r, betaT_f) %>% 
    mutate(ffmsy_bin=cut(ffmsy_avg, breaks=seq(0,4.5,0.5)),
           trend_bin=cut(sst_c_trend, breaks=seq(-0.03,0.07,0.01)))
  value_var <- ifelse(model=="random", "betaT_r", "betaT_f")
  data2 <- reshape2::dcast(data1, ffmsy_bin ~ trend_bin, value.var=value_var, fun.aggregate=mean)
  data3 <- t(as.matrix(data2[,2:ncol(data2)]))
  
  # Plot interaction
  par(mgp=c(1.6, 0.7, 0))
  colors <- brewer.pal(11, "RdBu")
  image(y=seq(0,4.5,0.5), x=seq(-0.03,0.07,0.01), z=data3, las=1,
        zlim=zlims[i,], col=colors, cex.axis=0.8,
        xlab="SST trend (°C/1yr)", ylab=expression("F/F"["MSY"]*" mean"), xpd=NA)
  lines(x=c(-0.03, 0.07), y=c(1, 1), lty=3)
  lines(x=c(0,0), y=c(0,4.5), lty=3)
  points(data1$sst_c_trend, data1$ffmsy_avg, col="grey40")
  
  
}

# Axis labels
mtext(expression("SST influence (θ"["i"]*")"), outer=T, side=2, adj=0.5, line=0, cex=0.8)

# Off
dev.off()
graphics.off()




