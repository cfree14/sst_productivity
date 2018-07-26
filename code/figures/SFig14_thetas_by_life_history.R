
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)


# Plot data
################################################################################

# Parameters
traits <- c("k", "linf_cm", "winf_kg", 
            "m", "tmax_yr", "tmat_yr", "lmat_cm", "troph", "depth_m_mid")
xlabels <- c("Brody growth coefficient (K)", expression("L"["inf"]*" (cm)"), expression("W"["inf"]*" (kg)"), 
             "Natural mortality (1/yr)", expression("T"["max"]*" (yr)"), expression("T"["mat"]*" (yr)"), 
             expression("L"["mat"]*" (cm)"), "Trophic level", "Depth (m)")

# Specify xmax's
apply(data[,traits], 2, min, na.rm=T)
apply(data[,traits], 2, max, na.rm=T)
xmins <- c(0,0,0,0,0,0,0,2.5,0)
xmaxs <- c(5, 500, 500, 1.5, 100, 20, 200, 4.5,1600)
xbins <- c(0.2, 50, 50, 0.5, 20, 5, 50, 0.5, 400)

# Setup figure
figname <- "SFig14_thetas_by_life_history.png"
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

