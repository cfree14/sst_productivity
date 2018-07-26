

# Clear workspace
rm(list = ls())

# Turn off scientific notation
options(scipen=999)

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Read SST data
cobe <- read.csv(paste(sstdir, "ramldb_sst_yearly_cobe.csv", sep="/"), as.is=T)
ersst <- read.csv(paste(sstdir, "ramldb_sst_yearly_ersst.csv", sep="/"), as.is=T)
hadsst <- read.csv(paste(sstdir, "ramldb_sst_yearly_hadsst.csv", sep="/"), as.is=T)

# Read model results
stocks <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)


# Build data
################################################################################

# Build data
data <- cobe %>% 
  rename(cobe_sst_c=sst_c) %>% 
  # Add ERSST data
  left_join(ersst, by=c("assessid", "year")) %>% 
  rename(er_sst_c=sst_c) %>% 
  # Add HadISST data
  left_join(hadsst, by=c("assessid", "year")) %>% 
  rename(had_sst_c=sst_c) %>% 
  # Filter to stocks of interest
  filter(assessid %in% stocks$assessid) %>% 
  filter(year<=2015) %>% 
  # Add stockid and rearrange
  left_join(select(stocks, assessid, stockid), by="assessid") %>% 
  select(assessid, stockid, year, everything())
  


# Plot data
################################################################################

# For y-axis label
top.i <- seq(1, nrow(stocks), 24)

# Setup figure
figname <- "AppendixD_sst_time_series.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 4), mar=c(2.5, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,6,3,3), lwd=0.8)

# Loop through stocks: i <- 1
for(i in 1:nrow(stocks)){
  
  # Subset data
  stock <- stocks$stockid[i]
  sdata <- subset(data, stockid==stock)
  
  # Plot HadSST data
  tmin <- floor(min(sdata[,4:6], na.rm=T))
  tmax <- ceiling(max(sdata[,4:6], na.rm=T))
  plot(had_sst_c ~ year, sdata, bty="n", type="l", las=1,
    ylim=c(tmin, tmax), xlab="", ylab="", main=stock, col="red", lwd=0.6)
  
  # Add ERSST data
  lines(sdata$year, sdata$er_sst_c, col="blue", lwd=0.6)
  
  # Add COBE data
  lines(sdata$year, sdata$cobe_sst_c, col="green", lwd=0.6)
  
  # Add y-axis label each new page
  if(i%in%top.i){mtext("SST (°C)", outer=T, side=2, adj=0.5, line=0.8, cex=1)}
  
}

# Off
dev.off()
graphics.off()



