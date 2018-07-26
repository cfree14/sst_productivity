

# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Read data
data_orig <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)


# LME
################################################################################

# LME stats
lme_n <- data_orig %>% 
  group_by(lme_name) %>% 
  summarize(betaT=median(betaT)) %>% 
  arrange(desc(betaT))
data_orig$lme_name <- factor(data_orig$lme_name, levels=lme_n$lme_name)

# Setup figure
figname <- "SFigX_thetas_by_lme.png"
png(paste(plotdir, figname, sep="/"), width=6, height=5, units="in", res=600)
par(mfrow=c(1,1), mar=c(3, 12, 0.8, 0.5), mgp=c(2,0.8,0))

# Plot data
boxplot(betaT ~ lme_name, data_orig, las=2, cex.axis=0.7, horizontal=T, lty=1,
        xlab=expression("SST influence (θ"["i"]*")"), cex=0.8, lwd=0.8, col="grey80")
abline(v=0, lty=3)

# Off
dev.off()

# Family
################################################################################

# Family stats
fam_n <- data_orig %>% 
  group_by(family) %>% 
  summarize(betaT=median(betaT)) %>% 
  arrange(desc(betaT))
data_orig$family <- factor(data_orig$family, levels=fam_n$family)

# Setup figure
figname <- "SFigX_thetas_by_family.png"
png(paste(plotdir, figname, sep="/"), width=6, height=5, units="in", res=600)
par(mfrow=c(1,1), mar=c(3, 12, 0.8, 0.5), mgp=c(2,0.8,0))

# Plot data
boxplot(betaT ~ family, data_orig, las=2, cex.axis=0.7, horizontal=T, lty=1,
        xlab=expression("SST influence (θ"["i"]*")"), cex=0.8, lwd=0.8, col="grey80")
abline(v=0, lty=3)

# Off
dev.off()


# LME
################################################################################

# LME stats
lme_n <- data_orig %>% 
  group_by(lme_name) %>% 
  summarize(betaT=median(betaT)) %>% 
  arrange(desc(betaT))
data_orig$lme_name <- factor(data_orig$lme_name, levels=lme_n$lme_name)

# Setup figure
figname <- "SFigX_thetas_by_lme.png"
png(paste(plotdir, figname, sep="/"), width=6, height=5, units="in", res=600)
par(mfrow=c(1,1), mar=c(3, 12, 0.8, 0.5), mgp=c(2,0.8,0))

# Plot data
boxplot(betaT ~ lme_name, data_orig, las=2, cex.axis=0.7, horizontal=T, lty=1,
        xlab=expression("SST influence (θ"["i"]*")"), cex=0.8, lwd=0.8, col="grey80")
abline(v=0, lty=3)

# Off
dev.off()

# FAO area
################################################################################

# FAO area stats
fao_n <- data_orig %>% 
  group_by(fao_area) %>% 
  summarize(betaT=median(betaT)) %>% 
  arrange(desc(betaT))
data_orig$fao_area <- factor(data_orig$fao_area, levels=fao_n$fao_area)

# Setup figure
figname <- "SFigX_thetas_by_fao_area.png"
png(paste(plotdir, figname, sep="/"), width=6, height=5, units="in", res=600)
par(mfrow=c(1,1), mar=c(3, 12, 0.8, 0.5), mgp=c(2,0.8,0))

# Plot data
boxplot(betaT ~ fao_area, data_orig, las=2, cex.axis=0.7, horizontal=T, lty=1,
        xlab=expression("SST influence (θ"["i"]*")"), cex=0.8, lwd=0.8, col="grey80")
abline(v=0, lty=3)

# Off
dev.off()


# Order
################################################################################

# Family stats
ord_n <- data_orig %>% 
  group_by(order) %>% 
  summarize(betaT=median(betaT)) %>% 
  arrange(desc(betaT))
data_orig$order <- factor(data_orig$order, levels=ord_n$order)

# Setup figure
figname <- "SFigX_thetas_by_order.png"
png(paste(plotdir, figname, sep="/"), width=6, height=5, units="in", res=600)
par(mfrow=c(1,1), mar=c(3, 12, 0.8, 0.5), mgp=c(2,0.8,0))

# Plot data
boxplot(betaT ~ order, data_orig, las=2, cex.axis=0.7, horizontal=T, lty=1,
        xlab=expression("SST influence (θ"["i"]*")"), cex=0.8, lwd=0.8, col="grey80")
abline(v=0, lty=3)

# Off
dev.off()




