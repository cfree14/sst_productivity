
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
library(freeR)
library(reshape2)
library(RColorBrewer)
library(fields) # colorbar.plot()

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/tables"
lmedir1 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme"
lmedir2 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme/lmes"
faodir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/fao_stat_areas/gis_data"

# Read data
data <- read.csv(paste(tabledir, "STable6_lme_msy_hindcast.csv", sep="/"), as.is=T) # FROM: AppendixH_msy_hindcasts_lme.R
stocks <- read.csv(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.csv", sep="/"), as.is=T)

# Read LME data
lme_key <- read.csv(paste(lmedir1, "lme_hsa_merge_key.csv", sep="/"), as.is=T)
lmes <- readOGR(dsn=lmedir2, laye="LME66_Offshore", verbose=F)

# Read FAO data
fao_areas <- readOGR(dsn=faodir, laye="FAO_AREAS_NOCOASTLINE", verbose=F)
fao_areas <- subset(fao_areas, F_LEVEL=="MAJOR")
# plot(fao_areas)


# Calculate MSY at median temperature for Reviewer 2
################################################################################

# Calculate MSY at average temperature and at median temperature
p <- 0.55
div <- (p+1)^((p+1)/p)
stocks <- stocks %>% 
  mutate(k1 = k * tb_max,
         msy_avg = (r * k1) / div * exp(betaT * (sst_c_avg-sst_c_avg)),
         msy_med = (r * k1) / div * exp(betaT * (sst_c_median-sst_c_avg)))

# Compares stock-level MSYs at average and median temperatures
plot(msy_med ~ msy_avg, stocks)

# Calculate LME-scale MSYs at median temperature
lme_msy <- stocks %>% 
  group_by(lme_name) %>% 
  summarize(msy_tmt_med=sum(msy_med)/1000,
            msy_tmt_avg=sum(msy_avg)/1000)

# Add LME-scale MSYs at median temperature to data
data1 <- data %>% 
  left_join(lme_msy, by=c("lme"="lme_name"))

# Compare LME-scale average MSYs calculated here with those calculated before to confirm that LME-scale median was calc'ed properly
plot(msy_tmt_avg ~ msy_tmt, data1)

# Setup figure
figname <- "msy_at_mean_median_temp_corr_reviewer2.png"
png(paste(plotdir, figname, sep="/"), width=4, height=4, units="in", res=600)
par(mar=c(3.5,3.5,0.5,0.5), mgp=c(2,0.8,0))

# Compare LME-scale MSY and average/median temps
options(scipen = 999)
plot(msy_tmt_med ~ msy_tmt_avg, data1, bty="n", 
     xlab="MSY at average temperature (1000s mt)", ylab="MSY at median temperature (1000s mt)")

# off
dev.off()

# Format data for plotting
################################################################################

# Add colors
summary(data1$pdiff)
hist(data1$pdiff, n=15)
data1$perc_bin <- cut(data1$pdiff, breaks=seq(-35,35,5))
colors_msy <- colorRampPalette(brewer.pal(11,"RdBu"))(nlevels(data1$perc_bin))
colors_msy_tr <- rgb(t(col2rgb(colors_msy))/255, alpha=0.7)
data1$perc_bin_color <- colors_msy_tr[data1$perc_bin]

# Add plotting position (pos) and plotting label to data
data1 <- data1 %>% 
  select(-lat_dd) %>% 
  left_join(select(lme_key, name, lat_dd, long_dd), by=c("lme"="name")) %>% 
  mutate(pos=revalue(lme, c("Agulhas Current"=4,
                            "Baltic Sea"=4,
                            "Barents Sea"=3,
                            "Bay of Biscay"=2, 
                            "Benguela Current"=2,
                            "California Current"=2,
                            "Canary Current"=4,
                            "Celtic-Biscay Shelf"=4,
                            "East Bering Sea"=1,
                            "East China Sea"=2,
                            "Faroe Plateau"=2,
                            "Greenland Sea"=3,
                            "Gulf of Alaska"=1,
                            "Gulf of Mexico"=1,
                            "Humboldt Current"=2,
                            "Iberian Coastal" =2,
                            "Iceland Shelf and Sea"=2,
                            "Indian Ocean"=2,
                            "Kuroshio Current"=4,
                            "Labrador - Newfoundland"=2,
                            "Labrador Sea"=2,
                            "Mediterranean Sea"=4,
                            "New Zealand Shelf"=3,
                            "North Atlantic Ocean"=4,
                            "North Brazil Shelf"=4,
                            "North Pacific Ocean"=2,
                            "North Sea"=2,
                            "Northeast U.S. Continental Shelf"=2,
                            "Norwegian Sea"=4,
                            "Patagonian Shelf"=4,
                            "Scotian Shelf"=2,
                            "Sea of Japan"=4,
                            "South Atlantic Ocean"=2,
                            "South Pacific Ocean"=2,
                            "South West Australian Shelf"=3,
                            "Southeast Australian Shelf"=2,
                            "Southeast U.S. Continental Shelf"=2,
                            "West Bering Sea"=3)),
         label=revalue(lme, c("Agulhas Current"="Agulhas Current",
                              "Baltic Sea"="Baltic Sea",
                              "Barents Sea"="Barents Sea",
                              "Benguela Current"="Benguela Current",
                              "California Current"="CA Current",
                              "Canary Current"="Canary Current",
                              "Celtic-Biscay Shelf"="Celtic-Biscay Shelf",
                              "East Bering Sea"="E Bering Sea",
                              "East China Sea"="E China Sea",
                              "Faroe Plateau"="Faroe Plateau",
                              "Greenland Sea"="Greenland Sea",
                              "Gulf of Alaska"="Gulf of Alaska",
                              "Gulf of Mexico"="Gulf of Mexico",
                              "Humboldt Current"="Humboldt Current",
                              "Iberian Coastal"="Iberian Coastal",
                              "Iceland Shelf and Sea"="Iceland Shelf & Sea",
                              "Indian Ocean"="Indian Ocean",
                              "Kuroshio Current"="Kuroshio Current",
                              "Labrador - Newfoundland"="Labrador-Newfoundland",
                              "Mediterranean Sea"="Mediterranean Sea",
                              "New Zealand Shelf"="NZ Shelf",
                              "North Atlantic Ocean"="N Atlantic Ocean",
                              "North Pacific Ocean"="N Pacific Ocean",
                              "North Sea"="North Sea",
                              "Northeast U.S. Continental Shelf"="NE US Shelf",
                              "Norwegian Sea"="Norwegian Sea",
                              "Patagonian Shelf"="Patagonian Shelf",
                              "Sea of Japan"="Sea of Japan",
                              "South Atlantic Ocean"="S Atlantic Ocean",
                              "South Pacific Ocean"="S Pacific Ocean",
                              "South West Australian Shelf"="SW Austr. Shelf",
                              "Southeast Australian Shelf"="SE Austr. Shelf",
                              "Southeast U.S. Continental Shelf"="SE US Shelf",
                              "West Bering Sea"="W Bering Sea")))
  
# LMEs without data
lmes_wout <- lme_key %>% 
  filter(!name %in% data1$lme)


# Plot data
################################################################################

# Setup figure
figname <- "Science_Fig4_msy_hindcast_map_for_reviewer2.png"
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

# Add MSY change
# barplot(hist(data1$msy_tmt_med, breaks=c(0, 100, 500, 1000, 2000, 10000, 20000), plot=F)$counts)
msy_breaks <-  c(100, 500, 1000, 2000, 10000, 20000)
msys <- cut(data1$msy_tmt_med, breaks=c(0, msy_breaks))
wts <- seq(1.8,4,length.out=nlevels(msys))
msy_wts <- wts[msys]
points(x=data1$long_dd, y=data1$lat_dd, 
       pch=21, bg=data1$perc_bin_color, col="grey40", cex=msy_wts)

# LME names
text(x=data1$long_dd, y=data1$lat_dd, pos=data1$pos, offset=0.6,
     labels=data1$label, cex=0.52, font=3)

# Add number of stocks
text(x=data1$long_dd, y=data1$lat_dd, labels=data1$n, cex=0.5, font=2, col="black")

# MSY size legend
# points(x=rep(205, length(msy_breaks)), 
#        y=rep(40, length(msy_breaks)), cex=wts, lwd=0.6)
points(x=rep(205, length(msy_breaks)), 
       y=56+seq(0,5.5, length.out=length(msy_breaks)), cex=wts, lwd=0.6)
text(x=190, y=80, pos=4, labels="MSY\n1000s of mt", font=2, xpd=NA, cex=0.75)
text(x=210, y=42+seq(0,27, length.out=length(msy_breaks)),
     pos=4, labels=msy_breaks, xpd=NA, cex=0.6)

# Legend
plot(1:10, 1:10, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="")
colorbar.plot(x=0, y=1.25, adj.x=0, adj.y=0, horiz=F, col=colors_msy,
              strip=seq(-1.5,1.5, length.out=nlevels(data1$perc_bin)), 
              strip.width=0.3, strip.length=2, xpd=NA)
text(x=2, y=c(1.25, 3.625, 6.05), pos=4, labels=c("-35%", "0%", "35%"), cex=0.7)
text(x=-0.5, y=6.6, pos=4, labels="% change\nin MSY", font=2, cex=0.75)

# Off
dev.off()
graphics.off()




