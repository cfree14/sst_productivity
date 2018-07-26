
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

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"
lmedir1 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme"
lmedir2 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme/lmes"
faodir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/fao_stat_areas/gis_data"

# Read data
data <- read.csv(paste(datadir, "msy_hindcast_time_series_using_lme_model.csv", sep="/"), as.is=T)
stocks <- read.csv(paste(datadir, "ramldb_v3.8_spsst_cobe_lme.csv", sep="/"), as.is=T)

# Read LME data
lme_key <- read.csv(paste(lmedir1, "lme_hsa_merge_key.csv", sep="/"), as.is=T)
lmes <- readOGR(dsn=lmedir2, laye="LME66_Offshore", verbose=F)

# Read FAO data
fao_areas <- readOGR(dsn=faodir, laye="FAO_AREAS_NOCOASTLINE", verbose=F)
fao_areas <- subset(fao_areas, F_LEVEL=="MAJOR")
# plot(fao_areas)


# Build data
################################################################################

# MSY time series (with LME)
data1 <- data %>% 
  left_join(select(stocks, stockid, lme_name), by="stockid") %>% 
  select(lme_name, stockid, year, msy_real) %>% 
  rename(lme=lme_name, msy=msy_real)

# LME-year MSY totals
lme_yr_msy <- data1 %>% 
  group_by(lme, year) %>% 
  summarize(nstocks=n_distinct(stockid),
            msy=sum(msy))

# Summarize data
years1 <- 1930:1940
years2 <- 2000:2010
lme_msy <- lme_yr_msy %>% 
  group_by(lme) %>% 
  summarize(nstocks=unique(nstocks),
            msy1avg=mean(msy[year%in%years1]), 
            msy2avg=mean(msy[year%in%years2]),
            change_mt=msy2avg-msy1avg,
            change_perc=(msy2avg-msy1avg)/msy1avg*100) %>% 
  left_join(lme_key, by=c("lme"="name"))


# Format data for plotting
################################################################################

# Add colors
summary(lme_msy$change_perc)
hist(lme_msy$change_perc, n=15)
lme_msy$perc_bin <- cut(lme_msy$change_perc, breaks=seq(-65,65,10))
colors <- colorRampPalette(brewer.pal(11,"RdBu"))(nlevels(lme_msy$perc_bin))
colors_tr <- rgb(t(col2rgb(colors))/255, alpha=0.7)
lme_msy$perc_bin_color <- colors_tr[lme_msy$perc_bin]

# Add plotting position (pos) and plotting label to data
lme_msy <- lme_msy %>% 
  mutate(pos=revalue(lme, c("Agulhas Current"=4,
                            "Baltic Sea"=4,
                            "Barents Sea"=3,
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
                            "Mediterranean Sea"=4,
                            "New Zealand Shelf"=3,
                            "North Atlantic Ocean"=4,
                            "North Pacific Ocean"=2,
                            "North Sea"=2,
                            "Northeast U.S. Continental Shelf"=2,
                            "Norwegian Sea"=4,
                            "Patagonian Shelf"=4,
                            "Sea of Japan"=4,
                            "South Pacific Ocean"=2,
                            "South West Australian Shelf"=2,
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
                              "New Zealand Shelf"="New Zealand Shelf",
                              "North Atlantic Ocean"="N Atlantic Ocean",
                              "North Pacific Ocean"="N Pacific Ocean",
                              "North Sea"="North Sea",
                              "Northeast U.S. Continental Shelf"="NE US Shelf",
                              "Norwegian Sea"="Norwegian Sea",
                              "Patagonian Shelf"="Patagonian Shelf",
                              "Sea of Japan"="Sea of Japan",
                              "South Pacific Ocean"="S Pacific Ocean",
                              "South West Australian Shelf"="SW Austr. Shelf",
                              "Southeast Australian Shelf"="SE Austr. Shelf",
                              "Southeast U.S. Continental Shelf"="SE US Shelf",
                              "West Bering Sea"="W Bering Sea")))
  
# LMEs without data
lmes_wout <- lme_key %>% 
  filter(!name %in% lme_msy$lme)


# Plot data
################################################################################

# Setup figure
figname <- "Fig6_theta_msy_hindcast_map.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6, units="in", res=600)
par(mfrow=c(2,1), mar=c(0.1,0.1,0.1,0.1), xpd=NA)

# A. Theta map
####################################

# Randomize order for plotting
set.seed(1)
stocks <- stocks[sample(1:nrow(stocks)),]

# Add colors
summary(stocks$betaT)
stocks$betaT_bin <- cut(stocks$betaT, breaks=seq(-1.5,1.5,0.20))
colors <- colorpal(brewer.pal(11, "RdBu"), nlevels(stocks$betaT_bin))
colors_tr <- rgb(t(col2rgb(colors))/255, alpha=0.7)
stocks$betaT_bin_color <- colors_tr[stocks$betaT_bin]

# Plot FAO stat areas
plot(fao_areas, border="grey70", lty=3)

# Plot world countries
map("world", col="grey80", fill=T, border="white", lwd=0.3, 
    xlim=c(-180,180), ylim=c(-90,90), add=T)

# Add betaT values
points(stocks$long_dd, stocks$lat_dd, pch=21, cex=1.3,
       bg=stocks$betaT_bin_color, col="grey30")

# B. Hindcast by LME
####################################

# Plot LME boundaries
xlim <- c(-160, 170)
ylim <- c(-50, 80)
# plot(lmes, col="grey60", lty=3)
plot(fao_areas, border="grey70", lty=3)

# Plot world countries
map("world", col="grey80", fill=T, border="white", lwd=0.3,
    xlim=c(-180,180), ylim=c(-90,90), add=T)

# Add MSY change
msys <- lme_msy$msy1avg
msy_wts <- 1.5 + 1.6 * (msys - min(msys))/(max(msys) - min(msys))
points(x=lme_msy$long_dd, y=lme_msy$lat_dd, 
       pch=21, bg=lme_msy$perc_bin_color, col="grey30", cex=msy_wts)

# Check that LME names are correct (comment out for figure)
text(x=lme_msy$long_dd, y=lme_msy$lat_dd, pos=lme_msy$pos,
     labels=lme_msy$label, cex=0.4, font=3)

# Add number of stocks
text(x=lme_msy$long_dd, y=lme_msy$lat_dd, labels=lme_msy$nstocks, cex=0.45, font=2, col="black")

# Off
dev.off()
graphics.off()




