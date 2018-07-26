
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(maps)
library(mapdata)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"
lmedir1 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme"
lmedir2 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme/lme_hsa_merge"

# Read data
data <- read.csv(paste(datadir, "msy_hindcast_time_series.csv", sep="/"), as.is=T)
stocks <- read.csv(paste(datadir, "ramldb_v3.8_spmodel_temp_dpdt_sst_yr_t.csv", sep="/"), as.is=T)

# Read LME data
lme_key <- read.csv(paste(lmedir1, "lme_hsa_merge_key.csv", sep="/"), as.is=T)
# lmes <- readOGR(dsn=lmedir2, laye="lme_hsa_merge", verbose=F)

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
years1 <- 1941:1945
years2 <- 2011:2015
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
lme_msy$perc_bin <- cut(lme_msy$change_perc, breaks=seq(-30,30,5))
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
                            "Sea of Okhotsk"=2,
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
                              "Sea of Okhotsk"="Sea of Okhotsk",
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
figname <- "Fig6_msy_hindcast_lme.png"
png(paste(plotdir, figname, sep="/"), width=6, height=3, units="in", res=600)
par(mar=c(0.5,0.5,0.5,0.5), xpd=NA)

# Plot LME boundaries
xlim <- c(-160, 170)
ylim <- c(-50, 80)
# plot(lmes, border="grey60", lty=3, xlim=xlim, ylim=ylim)

# Plot world countries
map("world", col="grey80", fill=T, border="white", lwd=0.3,
    xlim=xlim, ylim=ylim)

# Add MSY change
msys <- lme_msy$msy1avg
msy_wts <- 1.5 + 1.6 * (msys - min(msys))/(max(msys) - min(msys))
points(x=lme_msy$long_dd, y=lme_msy$lat_dd, 
       pch=21, bg=lme_msy$perc_bin_color, col="grey30", cex=msy_wts)

# Check that LME names are correct (comment out for figure)
text(x=lme_msy$long_dd, y=lme_msy$lat_dd, pos=lme_msy$pos,
     labels=lme_msy$label, cex=0.4, font=3)

# Add number of stocks
text(x=lme_msy$long_dd, y=lme_msy$lat_dd, labels=lme_msy$nstocks, cex=0.4, font=2, col="grey30")

# Off
dev.off()
graphics.off()




