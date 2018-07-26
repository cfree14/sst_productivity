

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(maps)
library(mapdata)
library(rgdal)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

# Directories
faodir1 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/trade_and_collapse/data/fao_landings/data"
faodir2 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/trade_and_collapse/data/fao_landings/predictions"
faodir3 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/fao_stat_areas/gis_data"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Read status predictions
stocks_orig <- read.csv(file.path(faodir1, "FAO_stocks_used_in_analysis.csv"), as.is=T)
centroids <- read.csv(file.path(faodir1, "FAO_ISO3_stock_centroids.csv"), as.is=T)
data <- read.csv(file.path(faodir2, "1950_2017_FAO_bbmsy_timeseries_cmsy13.csv"), as.is=T)

# Read SST data
sst <- read.csv(file.path(sstdir, "fao_area_sst_yearly_cobe.csv"), as.is=T)

# Read FAO data
fao_areas <- readOGR(dsn=faodir3, laye="FAO_AREAS_NOCOASTLINE", verbose=F)
fao_areas <- subset(fao_areas, F_LEVEL=="MAJOR")

# Build data
################################################################################

# Mean SST last 25 years (by area)
s_temp <- sst %>% 
  filter(!is.na(sst_c)) %>% 
  group_by(areaid) %>% 
  summarize(sst_c_avg=mean(sst_c[year>1990], na.rm=T),
            sst_c_trend=coef(lm(sst_c ~ year, subset=year>1990))[2]*10)

# Mean status last 25 years (by stock)
s_status <- data %>% 
  group_by(stockid) %>% 
  summarize(bbmsy_all=mean(bbmsy_q50),
            bbmsy_25yr=mean(bbmsy_q50[year>1990]),
            catch_25yr=mean(catch[year>1990]))

# Build stock data
stocks <- stocks_orig %>% 
  # Add area id
  mutate(areaid=paste(fao_code, iso3, sep="-")) %>% 
  # Add mean status
  left_join(s_status, by="stockid") %>% 
  # Add mean SST
  left_join(s_temp, by="areaid") %>% 
  # Add species SST percentile
  group_by(spp_code) %>% 
  mutate(spp_sst_q=ecdf(sst_c_avg)(sst_c_avg),
         spp_sst_bin=cut(spp_sst_q, breaks=seq(0,1,0.1)),
         spp_sst_col=rev(brewer.pal(nlevels(spp_sst_bin), "RdBu"))[spp_sst_bin]) %>% 
  # Mark "vulnerable stocks"
  mutate(vuln=bbmsy_25yr<0.5 & sst_c_trend>0 & spp_sst_q>0.8)


# Histograms of status
hist(s_status$bbmsy_all, breaks=seq(0,2,0.1), main="", xlab="B/BMSY mean", col="grey60", las=1)
hist(s_status$bbmsy_25yr, breaks=seq(0,2,0.1), main="", xlab="B/BMSY mean", col="grey60", las=1)
  

# Plot interaction
################################################################################

# Setup figure
figname <- "SFig20_vuln_stocks_interaction.png"
png(paste(plotdir, figname, sep="/"), width=4.5, height=4.5, units="in", res=600)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,0.5), mgp=c(2,0.8,0))

# Plot interaction
plot(bbmsy_25yr ~ sst_c_trend, stocks, bty="n", las=1,
     col=spp_sst_col,
     xlim=c(-0.4, 0.8), ylim=c(0,2), 
     xlab="Temperature trend (°C/10yr)", 
     ylab=expression("Mean B/B"["MSY"]))
lines(x=c(-0.4,0.8), y=c(0.5,0.5), lty=3, lwd=1.5)
lines(x=c(0,0), y=c(0,2), lty=3, lwd=1.5)

# Off
dev.off()
graphics.off()


# Plot vulnerability distribution
################################################################################

# Number of vulnerable stocks by country
stats <- stocks %>% 
  # Remove stocks that don't have B/BMSY predictions
  filter(!is.na(bbmsy_25yr)) %>% 
  # Number and catch of vulnerable stocks
  group_by(areaid) %>% 
  summarize(n=n(),
            n_vuln=sum(vuln),
            tc=sum(catch_25yr)/1000,
            tc_vuln=sum(catch_25yr[vuln==T])/1000) %>% 
  # Add area centroids
  left_join(centroids, by="areaid")

# Percent vulnerable
sum(stats$n_vuln) / sum(stats$n) * 100


# Setup figure
figname <- "SFig21_vuln_stocks_by_country.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3, units="in", res=600)
layout(matrix(c(1,2), ncol=2, byrow=T), widths=c(0.88,0.12))
par(mar=c(0.1,0.1,0.1,0.1), xpd=NA)

# Parameters
# xlim <- c(-160, 170)
# ylim <- c(-50, 80)
xlim <- c(-180, 180)
ylim <- c(-90, 90)

# Plot LME boundaries
# plot(lmes, col="grey60", lty=3)
plot(fao_areas, border="grey70", lty=3, xlim=xlim, ylim=ylim)

# Plot world countries
map("world", col="grey85", fill=T, border="white", lwd=0.3,
    xlim=xlim, ylim=ylim, add=T)

# Add points showing number/catch of vulnerable stocks
quantile(stats$tc_vuln)
# barplot(hist(stats$tc_vuln, breaks=c(0, 1, 5, 10, 50, 100, 500, 1000), plot=F)$counts)
tc_breaks <-  c(1, 5, 10, 50, 100, 500, 1000)
tcs <- cut(stats$tc_vuln, breaks=c(0, tc_breaks))
wts <- seq(1.8,4,length.out=nlevels(tcs))
cols <- freeR::tcolor(brewer.pal(length(wts), "YlOrRd"), 0.5)[tcs]
tcs_wts <- wts[tcs]
points(x=stats$long_dd, y=stats$lat_dd, 
       pch=21, bg=cols, col="grey40", cex=tcs_wts)

# Add number of stocks
text(x=stats$long_dd, y=stats$lat_dd, labels=stats$n_vuln, cex=0.5, font=2, col="black")

# MSY size legend
# points(x=rep(205, length(tc_breaks)), 
#        y=rep(40, length(tc_breaks)), cex=wts, lwd=0.6)
points(x=rep(205, length(tc_breaks)), pch=21, bg=rev(brewer.pal(length(wts), "YlOrRd")),
       y=rev(56+seq(0,5.5, length.out=length(tc_breaks))), cex=rev(wts), lwd=0.6)
text(x=190, y=80, pos=4, labels="Catch\n1000s of mt", font=2, xpd=NA, cex=0.75)
text(x=210, y=42+seq(0,27, length.out=length(tc_breaks)),
     pos=4, labels=tc_breaks, xpd=NA, cex=0.6)

# Legend
# plot(1:10, 1:10, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="")
# colorbar.plot(x=0, y=1.25, adj.x=0, adj.y=0, horiz=F, col=colors_msy,
#               strip=seq(-1.5,1.5, length.out=nlevels(data$perc_bin)), 
#               strip.width=0.3, strip.length=2, xpd=NA)
# text(x=2, y=c(1.25, 3.625, 6.05), pos=4, labels=c("-70%", "0%", "70%"), cex=0.7)
# text(x=-0.5, y=6.6, pos=4, labels="% change\nin MSY", font=2, cex=0.75)

# Off
dev.off()
graphics.off()





  



