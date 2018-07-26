
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
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/tables"
lmedir1 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme"
lmedir2 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme/lmes"
faodir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/fao_stat_areas/gis_data"

# Read data
data <- read.csv(paste(tabledir, "STable6_lme_msy_hindcast.csv", sep="/"), as.is=T)
stocks <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)

# Read LME data
lme_key <- read.csv(paste(lmedir1, "lme_hsa_merge_key.csv", sep="/"), as.is=T)
lmes <- readOGR(dsn=lmedir2, laye="LME66_Offshore", verbose=F)

# Read FAO data
fao_areas <- readOGR(dsn=faodir, laye="FAO_AREAS_NOCOASTLINE", verbose=F)
fao_areas <- subset(fao_areas, F_LEVEL=="MAJOR")
# plot(fao_areas)


# Build data
################################################################################

# # MSY time series (with LME)
# data1 <- data %>% 
#   left_join(select(stocks, stockid, lme_name), by="stockid") %>% 
#   select(lme_name, stockid, year, msy_real) %>% 
#   rename(lme=lme_name, msy=msy_real)
# 
# # LME-year MSY totals
# lme_yr_msy <- data1 %>% 
#   group_by(lme, year) %>% 
#   summarize(nstocks=n_distinct(stockid),
#             msy=sum(msy))
# 
# # Summarize data
# years1 <- 1930:1939
# years2 <- 2001:2010
# lme_msy <- lme_yr_msy %>% 
#   group_by(lme) %>% 
#   summarize(nstocks=unique(nstocks),
#             msy1avg=mean(msy[year%in%years1]), 
#             msy2avg=mean(msy[year%in%years2]),
#             change_mt=msy2avg-msy1avg,
#             change_perc=(msy2avg-msy1avg)/msy1avg*100) %>% 
#   left_join(lme_key, by=c("lme"="name"))
# 

# Format data for plotting
################################################################################

# Add colors
summary(data$pdiff)
hist(data$pdiff, n=15)
data$perc_bin <- cut(data$pdiff, breaks=seq(-70,70,10))
colors_msy <- colorRampPalette(brewer.pal(11,"RdBu"))(nlevels(data$perc_bin))
colors_msy_tr <- rgb(t(col2rgb(colors_msy))/255, alpha=0.7)
data$perc_bin_color <- colors_msy_tr[data$perc_bin]

# Add plotting position (pos) and plotting label to data
data <- data %>% 
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
                              "South Atlantic Ocean"="S Atlantic Ocean",
                              "South Pacific Ocean"="S Pacific Ocean",
                              "South West Australian Shelf"="SW Austr. Shelf",
                              "Southeast Australian Shelf"="SE Austr. Shelf",
                              "Southeast U.S. Continental Shelf"="SE US Shelf",
                              "West Bering Sea"="W Bering Sea")))
  
# LMEs without data
lmes_wout <- lme_key %>% 
  filter(!name %in% data$lme)


# Plot data
################################################################################

# Setup figure
figname <- "Fig6_theta_msy_hindcast_map.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6, units="in", res=600)
layout(matrix(c(1,2,3,4), ncol=2, byrow=T), widths=c(0.9,0.1))
par(mar=c(0.1,0.1,0.1,0.1), xpd=NA)

# Parameters
# xlim <- c(-160, 170)
# ylim <- c(-50, 80)
xlim <- c(-180, 180)
ylim <- c(-90, 90)

# A. Theta map
####################################

# Randomize order for plotting
set.seed(1)
stocks <- stocks[sample(1:nrow(stocks)),]

# Add colors
summary(stocks$betaT)
betaT_breaks <- seq(-1.25,1.25,0.25)
stocks$betaT_bin <- cut(stocks$betaT, breaks=betaT_breaks)
colors <- colorpal(brewer.pal(11, "RdBu"), nlevels(stocks$betaT_bin))
colors_tr <- rgb(t(col2rgb(colors))/255, alpha=0.7)
stocks$betaT_bin_color <- colors_tr[stocks$betaT_bin]

# Plot FAO stat areas
plot(fao_areas, border="grey70", lty=3, xlim=xlim, ylim=ylim)
text(x=-172, y=82, labels="A", font=2, cex=1.2)

# Plot world countries
map("world", col="grey85", fill=T, border="white", lwd=0.3, 
    xlim=xlim, ylim=ylim, add=T)

# Add betaT values
stocks_look <- select(stocks, stockid, long_dd, lat_dd)
stocks_on_land <- c("SAILWATL", "SKJWATL", "SAILEATL", "SKJEATL")
stocks_plot <- subset(stocks, !stockid %in% stocks_on_land)
points(jitter(stocks_plot$long_dd,120), jitter(stocks_plot$lat_dd,120), pch=21, cex=1.6,
       bg=stocks$betaT_bin_color, col="grey30")

# Legend
plot(1:10, 1:10, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="")
colorbar.plot(x=0, y=1, adj.x=0, adj.y=0, horiz=F, col=colors,
              strip=seq(-1.5,1.5,length.out=nlevels(stocks$betaT_bin)), 
                        strip.width=0.3, strip.length=2, xpd=NA)
text(x=2, y=c(1.1, 2.95, 4.85), pos=4, labels=c("-1.25", "0.0", "1.25"), cex=0.8)
text(x=-0.5, y=5.5, pos=4, labels="SST\ninfluence", font=2, cex=0.8)

# B. Hindcast by LME
####################################

# # Plot LME boundaries
# # plot(lmes, col="grey60", lty=3)
# plot(fao_areas, border="grey70", lty=3, xlim=xlim, ylim=ylim)
# text(x=-172, y=82, labels="B", font=2, cex=1.2)
# 
# # Plot world countries
# map("world", col="grey85", fill=T, border="white", lwd=0.3,
#     xlim=xlim, ylim=ylim, add=T)
# 
# # Add MSY change
# # barplot(hist(data$msy_tmt, breaks=c(0, 100, 500, 1000, 2000, 10000, 20000), plot=F)$counts)
# msy_breaks <-  c(100, 500, 1000, 2000, 10000, 20000)
# msys <- cut(data$msy_tmt, breaks=c(0, msy_breaks))
# wts <- seq(1.8,4,length.out=nlevels(msys))
# msy_wts <- wts[msys]
# points(x=data$long_dd, y=data$lat_dd, 
#        pch=21, bg=data$perc_bin_color, col="grey40", cex=msy_wts)
# 
# # LME names
# text(x=data$long_dd, y=data$lat_dd, pos=data$pos, offset=0.6,
#      labels=data$label, cex=0.52, font=3)
# 
# # Add number of stocks
# text(x=data$long_dd, y=data$lat_dd, labels=data$n, cex=0.5, font=2, col="black")
# 
# # MSY size legend
# # points(x=rep(205, length(msy_breaks)), 
# #        y=rep(40, length(msy_breaks)), cex=wts, lwd=0.6)
# points(x=rep(205, length(msy_breaks)), 
#        y=40+seq(0,4.6, length.out=length(msy_breaks)), cex=wts, lwd=0.6)
# text(x=190, y=60, pos=4, labels="MSY\n1000s of mt", font=2, xpd=NA, cex=0.8)
# text(x=210, y=23+seq(0,27, length.out=length(msy_breaks)),
#      pos=4, labels=msy_breaks, xpd=NA, cex=0.7)
# 
# # Legend
# plot(1:10, 1:10, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="")
# colorbar.plot(x=0, y=1, adj.x=0, adj.y=0, horiz=F, col=colors_msy,
#               strip=seq(-1.5,1.5, length.out=nlevels(data$perc_bin)), 
#               strip.width=0.3, strip.length=2, xpd=NA)
# text(x=2, y=c(1.1, 2.95, 4.85), pos=4, labels=c("-70%", "0%", "70%"), cex=0.8)
# text(x=-0.5, y=5.5, pos=4, labels="% change\nin MSY", font=2, cex=0.8)

# Off
dev.off()
graphics.off()




