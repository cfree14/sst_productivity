
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(mblm) # For Thiel-Sen slope
library(plyr)
library(dplyr)
library(freeR)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/tables"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"
lmedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme"

# Read LME model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
rm(hess, input.data, model, output, sd, nstocks, params, problem.stocks, results.wide, data)

# Read hindcast data
load(paste(datadir, "spsst_pella40%_cobe_lme_msy_hindcast_10000traj_mvnorm.Rdata", sep="/"))

# Read SST data
sst <- read.csv(paste(sstdir, "lme_hsa_sst_yearly_cobe.csv", sep="/"), as.is=T)

# Read LME/HSA key
lme_key <- read.csv(paste(lmedir, "lme_hsa_merge_key.csv", sep="/"), as.is=T)

# Format data
################################################################################

# Format LME model hyper-means
lme_infs <- results.df %>% 
  filter(param=="mu_group") %>% 
  select(-param) %>% 
  rename(lme=stockid, betaT=estimate, betaT_se=stderror)

# Format SST stats
sst_stats <- sst %>%
  filter(year %in% 1930:2010) %>% 
  group_by(lme) %>% 
  summarize(sst_c_mean=mean(sst_c),
            sst_c_slope_10yr=coef(lm(sst_c ~ year))[2]*10)

# Format MSY hindcast data
# Add plotting position (pos) and plotting label to data
msy_lme_df <- msy_lme_df %>% 
  mutate(label=revalue(lme, c("Iceland Shelf and Sea"="Iceland Shelf & Sea",
                              "Labrador - Newfoundland"="Labrador-Nwfld.",
                              "North Atlantic Ocean"="N Atlantic Ocean",
                              "North Pacific Ocean"="N Pacific Ocean",
                              "Northeast U.S. Continental Shelf"="NE US Contl. Shelf",
                              "South Atlantic Ocean"="S Atlantic Ocean",
                              "South Pacific Ocean"="S Pacific Ocean",
                              "South West Australian Shelf"="SW Australia Shelf",
                              "Southeast Australian Shelf"="SE Australia Shelf",
                              "Southeast U.S. Continental Shelf"="SE US Contl. Shelf")))


# Plot MSY hindcasts
################################################################################

# LMEs
lmes <- sort(unique(stocks$lme_name))

# Build stats data frame
stats <- stocks %>% 
  # LME sample size and MSY sum
  rename(lme=lme_name) %>% 
  group_by(lme) %>% 
  summarize(n=n(),
            msy_tmt=sum(msy_avg)/1000) %>% 
  # Add LME type
  left_join(select(lme_key, name, type, lat_dd), by=c("lme"="name")) %>% 
  # Add SST stats
  left_join(sst_stats, by="lme") %>% 
  # Add SST influence
  left_join(select(lme_infs, lme, betaT), by="lme") %>% 
  # Add MSY change stats
  mutate(ts_slope=NA,
         ts_pval=NA,
         ts_slope_sd=NA,
         ts_pval_sd=NA,
         pdiff=NA) %>% 
  # Rearrange columns
  select(type, lme, n, msy_tmt, lat_dd,
         sst_c_mean, sst_c_slope_10yr, betaT, 
         ts_slope, ts_pval, ts_slope_sd, ts_pval_sd, pdiff) %>% 
  ungroup()

# Colors
display.brewer.pal(11, "RdBu")
colors <- brewer.pal(11, "RdBu")

# For y-axis label
top.i <- seq(1, nrow(stocks), 24)

# Setup figure
figname <- "AppendixH_msy_hindcasts_lme.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 4), mar=c(1, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,6,3,3), lwd=0.8)

# Loop through stocks: i <- 1
for(i in 1:length(lmes)){

  # Subset data
  lme1 <- lmes[i]
  lme_label <- unique(msy_lme_df$label[msy_lme_df$lme==lme1])
  sdata <- msy_lme_df %>% 
    filter(lme==lme1 & year %in% 1930:2010) %>% 
    mutate(year=as.numeric(year),
           msy_sd=msy/max(msy))

  # Plot data
  ymin <- floor1(min(sdata$msy_lo), 0.1)
  ymax <- ceiling1(max(sdata$msy_hi), 0.1)
  plot(msy ~ year, sdata, bty="n", type="n", xlab="", ylab="", las=1, 
       xlim=c(1930,2010), ylim=c(ymin, ymax), xaxt="n")
  axis(1, at=seq(1930,2010,10), las=2)
  
  # Add MSY shading and line
  polygon(x=c(sdata$year, rev(sdata$year)),
          y=c(sdata$msy_lo, rev(sdata$msy_hi)),
          col="grey80", border=F)
  lines(x=sdata$year, y=sdata$msy, col="black")
  
  # Percent difference and Thiel-Sen slope
  tsfit <- mblm(msy ~ year, sdata, repeated=F)
  tsfit_sd <- mblm(msy_sd ~ year, sdata, repeated=F)
  stats$ts_slope[i] <- coef(tsfit)[2]*10*1E6 # mt per decade (convert millions of mt to just mt)
  stats$ts_pval[i] <- summary(tsfit)$coefficients[2,4]
  stats$ts_slope_sd[i] <- coef(tsfit_sd)[2]*10
  stats$ts_pval_sd[i] <- summary(tsfit_sd)$coefficients[2,4]
  msy1 <- mean(sdata$msy[sdata$year%in%1930:1939])
  msy2 <- mean(sdata$msy[sdata$year%in%2001:2010])
  stats$pdiff[i] <- (msy2-msy1)/msy2*100
  
  # Add labels
  n <- sum(stocks$lme_name==lme1)
  lme_text <- paste0(lme_label, " (n=", n, ")")
  stat_text <- paste0(roundf(stats$pdiff[i],1), "% change\n", roundf(stats$ts_slope[i], 1), " mt / 10 yr")
  text(x=1928, y=ymax-(ymax-ymin)*0.05, pos=4, labels=lme_text, font=2, xpd=NA)
  text(x=1928, y=ymax-(ymax-ymin)*0.25, pos=4, labels=stat_text, xpd=NA)
  
  # Add SST
  par(new = T)
  sst1 <- subset(sst, lme==lme1 & year%in%1930:2010)
  tmin <- min(sst1$sst_c)
  tmax <- max(sst1$sst_c)
  ymax <- tmax + (tmax-tmin)*0.4
  plot(sst_c ~ year, sst1, type="l", col=tcolor("red", 0.6), lwd=0.7,
       bty="n", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(tmin, ymax))
  
  
  # Add y-axis label
  if(i%in%top.i){mtext("MSY (millions mt)", outer=T, side=2, adj=0.5, line=0.8)}
  

}

# Off
dev.off()
graphics.off()

# Export table
################################################################################

# Export table
write.csv(stats, paste(tabledir, "STable6_lme_msy_hindcast.csv", sep="/"), row.names=F)





