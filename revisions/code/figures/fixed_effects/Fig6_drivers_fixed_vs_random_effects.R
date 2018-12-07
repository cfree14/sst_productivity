
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
library(quantreg)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_re_fe_merged.csv", sep="/"), as.is=T)
data$sst_c_trend <- data$sst_c_trend*10

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
range(data$betaT_f)
range(data$betaT_r)
models <- c("fixed", "random")
ylims <- matrix(c(-5, 5,
                  -1, 1), ncol=2, byrow=T)

# Setup figure
figname <- "FE_Fig6_drivers_fixed_vs_random_effects.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3.8, units="in", res=600)
par(mfrow=c(2,3), mar=c(2.5,2.5,0.5,0.5), mgp=c(2,0.7,0), oma=c(1,2,0,0))

# Loop and plot
for(i in 1:length(models)){
  
  # Reset par
  # par(mgp=c(2,0.7,0))
  
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
  
  # Overfishing
  #####################################################
  
  # Plot BetaT ~ overfishing
  plot(x=data$ffmsy_avg, y=betaT, bty="n", las=1, col=colors, cex.axis=0.8,
       xlim=c(0,5), ylim=ylims[i,], xlab=xlabs[1], ylab="", xpd=NA)
  lines(x=c(0, 5), y=c(0,0), lty=3)
  mtext(title, side=3, adj=0.1, line=-1.3, cex=0.7, font=2)
  fit_plot_qr(x=data$ffmsy_avg, y=betaT, xlim=c(0,5))
  # if(i==1){mtext("A", side=3, adj=0.95, line=-2, cex=0.8, font=2)}
  
  # Maximum age
  #####################################################
  
  # Plot BetaT ~ maximum age
  plot(x=data$tmax_yr, y=betaT, bty="n", las=1, col=colors, cex.axis=0.8,
       xlim=c(0,100), ylim=ylims[i,], xlab=xlabs[2], ylab="", xpd=NA)
  lines(x=c(0, 100), y=c(0,0), lty=3)
  fit_plot_qr(x=data$tmax_yr, y=betaT, xlim=c(0,100))

  
  # Plot thermal niche
  #####################################################
  
  
  # Species
  spp <- c("Gadus morhua", "Clupea harengus")
  cols <- brewer.pal(4, "Set1")[3:4]
  
  # Setup empty
  plot(1:10, 1:10, type="n", bty="n", las=1,
       cex.axis=0.8, xpd=NA,
       xlim=c(0, 20), ylim=ylims[i,], 
       xlab=c("","Mean temperature (°C)")[i], ylab="")
  lines(x=c(0, 20), y=c(0,0), lty=3)
  
  if(i==2){
    legend("topright", bty="n", col=cols, pch=16, lty=1,  cex=0.7,
           legend=c("Atlantic cod (n=12)", "Atlantic herring (n=10)"))
  }
  
  # Loop through species
  for(j in 1:length(spp)){
    
    # Subset data
    sci_name <- spp[j]
    sdata <- subset(data, species==sci_name)
    
    # Add points
    betaT_cols <- c("betaT_f", "betaT_r")
    points(sdata$sst_c_avg, sdata[,betaT_cols[i]], pch=16, col=cols[j], cex=1.2)
    
    # Fit and plot Thiel-Sen slope
    if(i==1){
      tsfit <- mblm(betaT_f ~ sst_c_avg, sdata, repeated=F)
    }else{
      tsfit <- mblm(betaT_r ~ sst_c_avg, sdata, repeated=F)
    }
    # pvalue <- roundf(summary(tsfit)$coefficients[2,4],2)
    # lty <- ifelse(pvalue<0.1, 1, 2)
    curve(coef(tsfit)[1] + coef(tsfit)[2]*x, from=0, to=20, n=100, add=T, lty=1, col=cols[j], lwd=1.1)
    
  }
  
  
  
  # Plot overfishing*warming interaction
  #####################################################
  
  # # Add fixed colors
  # summary(data$betaT_f)
  # betaT_breaks <- c(-16,-8,-4,-3,-2,-1,0,1,2,3,4,8,16)
  # data$betaT_bin_f <- cut(data$betaT_f, breaks=betaT_breaks)
  # colors <- colorpal(brewer.pal(11, "RdBu"), nlevels(data$betaT_bin_f))
  # colors_tr <- rgb(t(col2rgb(colors))/255, alpha=0.7)
  # data$betaT_bin_color_f <- colors[data$betaT_bin_f]
  # 
  # # Add random colors
  # summary(data$betaT_r)
  # betaT_breaks <- seq(-0.75,0.75,0.25)
  # data$betaT_bin_r <- cut(data$betaT_r, breaks=betaT_breaks)
  # colors <- colorpal(brewer.pal(11, "RdBu"), nlevels(data$betaT_bin_r))
  # colors_tr <- rgb(t(col2rgb(colors))/255, alpha=0.7)
  # data$betaT_bin_color_r <- colors[data$betaT_bin_r]
  # 
  # # Assign colors
  # if(i==1){
  #   colors <- data$betaT_bin_color_f
  # }else{
  #   colors <- data$betaT_bin_color_r
  # }
  # 
  # # F/FMSY*warming interaction
  # ##########################################
  # 
  # # Plot F/FMSY * SST trend interaction
  # plot(ffmsy_avg ~ sst_c_trend, data, bty="n", las=1,
  #      cex.axis=1.1, cex.lab=1.1,
  #      xlim=c(-0.2, 0.8), ylim=c(0,5),
  #      xlab="Temperature trend (°C/10yr)", ylab=expression("Mean F/F"["MSY"]), 
  #      pch=21, cex=1.2, bg=colors, xpd=NA)
  # lines(x=c(-0.2, 0.8), y=c(1, 1), lty=3, lwd=1.2)
  # lines(x=c(0,0), y=c(0,5), lty=3, lwd=1.2)
  # 
  # # Add sample size
  # n1 <- nrow(filter(data, !is.na(ffmsy_avg) & !is.na(sst_c_trend)))
  # text(labels=paste0("n=", n1), x=-0.2+(0.8--0.2)*1.05, y=0, pos=2, cex=1, xpd=NA)
  
  
  
  
}

# Axis labels
mtext(expression("SST influence (θ"["i"]*")"), outer=T, side=2, adj=0.5, line=0, cex=0.8)

# Off
dev.off()
graphics.off()







