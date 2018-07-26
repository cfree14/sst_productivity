
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(quantreg)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)

# Format data
data <- data %>%
  mutate(log_area_sqkm=log(area_sqkm), 
         log_tb_avg=log(tb_avg))

# Plot data
################################################################################

# Parameters
traits <- c("lat_dd_abs", "sst_c_avg2", "sst_c_trend2", "log_area_sqkm", 
            "log_tb_avg", "tb_sd_trend",  "bbmsy_avg",  "ffmsy_avg", "n")
xlabels <- c("Latitude (° absolute)", "SST mean (°C)", "SST trend (°C/yr)", 
             expression("log(stock area (km"^2*"))"), "log(biomass mean (MT))", "Scaled biomass trend", 
             expression("B/B"["MSY"]*" mean"), expression("F/F"["MSY"]*" mean"), "Number of years")

# Specify xmax's
apply(data[,traits], 2, min, na.rm=T)
apply(data[,traits], 2, max, na.rm=T)
xmins <- c(0, 0, -0.015, 5, 5, -0.04, 0, 0, 0)
xmaxs <- c(80, 25, 0.015, 20, 20, 0.06, 7, 5, 150)

# Setup figure
figname <- "SFig13_thetas_by_stock_characteristics.png"
png(paste(plotdir, figname, sep="/"), width=6, height=6, units="in", res=600)
par(mfrow=c(3,3), mar=c(3.5, 2, 0.5, 0.5), mgp=c(2,0.8,0), oma=c(0,3,0,0))

# Loop through traits
for(i in 1:length(traits)){
  
  # Trait
  trait <- traits[i]
  
  # Subset data
  sdata <- data[, c("stockid", trait, "betaT", "betaT_inf")]
  colnames(sdata) <- c("stockid", "trait", "betaT", "betaT_inf")
  sdata <- sdata %>% 
    filter(!is.na(trait)) %>% 
    mutate(col=revalue(betaT_inf, c("positive"="blue", "negative"="red", "none"="grey60")))
  
  # Plot data
  xmin <- xmins[i]
  xmax <- xmaxs[i]
  plot(betaT ~ trait, sdata, bty="n", las=1, xpd=NA,
       xlim=c(xmin, xmax), ylim=c(-1.5,1.5), xlab=xlabels[i], ylab="",
       col=col, pch=1)
  lines(x=c(xmin, xmax), y=c(0,0), lty=3)
  
  # Add sample size text (if not 201)
  n <- nrow(sdata)
  if(n!=nrow(data)){text(x=xmax+(xmax-xmin)*0.05, y=-1.47, pos=2, labels=paste0("n=", n), font=2)}
  
  # Fit and plot lm
  # lmfit <- lm(betaT ~ trait, sdata)
  # curve(coef(lmfit)[1]+coef(lmfit)[2]*x, from=xmin, to=xmax, n=50, add=T, lwd=1.2)
  
  # Fit and plot quantile regression
  qrfit <- rq(betaT ~ trait, sdata, tau=0.5)
  qrfit_lo <- rq(betaT ~ trait, sdata, tau=0.025)
  qrfit_hi <- rq(betaT ~ trait, sdata, tau=0.975)
  curve(coef(qrfit)[1]+coef(qrfit)[2]*x, from=xmin, to=xmax, n=50, add=T, lwd=1.2)
  curve(coef(qrfit_lo)[1]+coef(qrfit_lo)[2]*x, from=xmin, to=xmax, n=50, add=T, lwd=1.2, lty=2)
  curve(coef(qrfit_hi)[1]+coef(qrfit_hi)[2]*x, from=xmin, to=xmax, n=50, add=T, lwd=1.2, lty=2)

}

# Add y-axis label
xlabel <- expression("SST influence (θ"["i"]*")")
mtext(xlabel, outer=T, side=2, adj=0.53, line=1, cex=0.9)

# Off
dev.off()
graphics.off()

