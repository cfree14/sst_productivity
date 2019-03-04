
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures/infographic"

# Read model data
load(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.Rdata", sep="/"))
spdata <- data; stocks <- results.wide
rm(data, hess, input.data, model, nstocks, output, params, problem.stocks, sd, results.df, results.wide)

# Read data
data_orig <- read.csv(paste(datadir, "msy_hindcast_time_series_using_lme_model.csv", sep="/"), as.is=T)


# Format data
################################################################################

# Hindcast period
yr1 <- 1930
yr2 <- 2010

# Format data
data <- data_orig %>%
  filter(year>=yr1 & year<=yr2) %>%
  mutate(msy_avg=msy_avg/1000,
         msy_real=msy_real/1000,
         msy_real_lo=msy_real_lo/1000,
         msy_real_hi=msy_real_hi/1000,
         msy_real_ma=zoo::rollapply(msy_real, width=10, FUN=mean, partial=T),
         sst_c_ma=zoo::rollapply(sst_c, width=10, FUN=mean, partial=T))

# Stats
stats <- data %>%
  filter(stockid %in% c("BSBASSMATLC", "CODIS")) %>% 
  group_by(stockid) %>% 
  summarize(msy1=mean(msy_real[year%in%1930:1939]),
            msy2=mean(msy_real[year%in%2001:2010]),
            msy_diff=(msy2-msy1)/msy1*100)


# Plot MSY hindcasts
################################################################################

# Colors
display.brewer.pal(11, "RdBu")
colors <- brewer.pal(11, "RdBu")

# Stocks to plot
stock <- "BSBASSMATLC"
stock <- "CODIS"

# Setup figure
figname <- paste0(stock, ".png")
png(file.path(plotdir, figname), width=5, height=4, units="in", res=600)
layout(mat=matrix(data=c(1,2), nrow=2), heights=c(0.6, 0.4))
par(mar=c(0.7, 5, 0.5, 0.5), mgp=c(2.5,0.8,0), xpd=NA, oma=c(2.5,0,0,0))


# Plot MSY over time
######################################

# Subset data
sdata <- subset(data, stockid==stock)
betaT <- format(round(unique(stocks$betaT[stocks$stockid==stock]),2), nsmall=2)
betaT_sig <- stocks$betaT_inf[stocks$stockid==stock]
msy_avg <- unique(sdata$msy_avg)

# Colors
if(betaT_sig=="none"){
  scolor <- ifelse(betaT<0, colors[5], colors[7])
  lcolor <- ifelse(betaT<0, colors[3], colors[9])
}else{
  scolor <- ifelse(betaT<0, colors[3], colors[9])
  lcolor <- ifelse(betaT<0, colors[1], colors[11])
}

# Plot MSY hindcast
ymin <- floor(min(sdata$msy_real_ma))
ymax <- ceiling(max(sdata$msy_real_ma))
plot(msy_real ~ year, sdata, type="n", bty="n", las=1,
     col="grey50", lwd=1, ylim=c(ymin, ymax), cex.lab=1.1,
     xaxt="n", yaxt="n", xlab="", ylab="Sustainable catch\n(1000s of tons)")
# axis(1, at=seq(yr1, yr2, 10), label=F, cex.axis=0.9, las=1)
axis(2, at=c(ymin, ymax), label=round(c(ymin, ymax), 1), cex.axis=1.1, las=1)

# Plot CI shading
# polygon(x=c(sdata$year, rev(sdata$year)),
#         y=c(sdata$msy_real_lo, rev(sdata$msy_real_hi)), col=scolor, border=F)

# Add MSY line
# lines(sdata$year, sdata$msy_real, col=lcolor)
lines(sdata$year, sdata$msy_real_ma, col=scolor, lwd=2)

# Add mean MSY line
lines(x=c(min(sdata$year), max(sdata$year)), y=rep(msy_avg,2), lty=3, lwd=0.8, col="grey40")


# Plot SST over time
######################################
  
# Add SST time series
tmin <- floor(min(sdata$sst_c_ma))
tmax <- ceiling(max(sdata$sst_c_ma))
sst_avg <- mean(sdata$sst_c_ma)
plot(sst_c_ma ~ year, sdata, type="l", col="black", lty=1, lwd=2, las=1, cex.lab=1.1,
     bty="n", xaxt="n", yaxt="n", xlab="", ylab="Ocean\ntemperature (°C)", ylim=c(tmin, tmax))
axis(1, at=seq(yr1, yr2, 10), label=T, cex.axis=1.1, las=2)
axis(2, at=c(tmin, tmax), label=T, las=1, cex.axis=1.1)
lines(x=c(min(sdata$year), max(sdata$year)), y=rep(sst_avg,2), lty=3, lwd=0.8, col="grey40")

# Off
dev.off()
graphics.off()

