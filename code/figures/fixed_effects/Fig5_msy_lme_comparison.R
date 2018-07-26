
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
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output/fixed_effects"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures/fixed_effects"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"
lmedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/lme"


# Function to build data
################################################################################

# Function to build data
build_data <- function(model_file, hindcast_file){

  # Read model
  load(paste(datadir, model_file, sep="/"))
  rm(hess, input.data, model, output, sd, nstocks, params, problem.stocks, results.wide, data)
  
  # Read hindcast data
  load(paste(datadir, hindcast_file, sep="/"))
  
  # Read SST data
  sst <- read.csv(paste(sstdir, "lme_hsa_sst_yearly_cobe.csv", sep="/"), as.is=T)
  
  # Read LME/HSA key
  lme_key <- read.csv(paste(lmedir, "lme_hsa_merge_key.csv", sep="/"), as.is=T)
  
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
    # Add MSY change stats
    mutate(ts_slope=NA,
           ts_pval=NA,
           ts_slope_sd=NA,
           ts_pval_sd=NA,
           pdiff=NA) %>% 
    # Rearrange columns
    select(type, lme, n, msy_tmt, lat_dd,
           sst_c_mean, sst_c_slope_10yr, 
           ts_slope, ts_pval, ts_slope_sd, ts_pval_sd, pdiff) %>% 
    ungroup()
  
  # Loop through stocks: i <- 1
  for(i in 1:length(lmes)){
  
    # Subset data
    lme1 <- lmes[i]
    lme_label <- unique(msy_lme_df$label[msy_lme_df$lme==lme1])
    sdata <- msy_lme_df %>% 
      filter(lme==lme1 & year %in% 1930:2010) %>% 
      mutate(year=as.numeric(year),
             msy_sd=msy/max(msy))
  
    
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
    
  }
  return(stats)
  
}


# Build data
################################################################################

# LME-scale change: fixed effects 
model_file <- "ramldb_v3.8_pella_0.20_cobe_fixed.Rdata"
hindcast_file <- "ramldb_v3.8_pella_0.20_cobe_fixed_msy_hindcast.Rdata"
stats_f <- build_data(model_file, hindcast_file)

# LME-scale change: random effects 
model_file <- "ramldb_v3.8_pella_0.20_cobe_random_normal.Rdata"
hindcast_file <- "ramldb_v3.8_pella_0.20_cobe_random_msy_hindcast.Rdata"
stats_r <- build_data(model_file, hindcast_file)

# Merge results
stats <- stats_f %>% 
  select(type, lme, n, msy_tmt, sst_c_mean, sst_c_slope_10yr, pdiff) %>% 
  rename(pdiff_f=pdiff) %>% 
  # Add random effects results
  left_join(select(stats_r, lme, pdiff), by="lme") %>% 
  rename(pdiff_r=pdiff)


# Plot data
################################################################################

# Setup figure
figname <- "Fig5_msy_lme_comparison_random_fixed.png"
png(paste(plotdir, figname, sep="/"), width=5, height=5, units="in", res=600)
par(mar=c(4,5,0.5,0.5), mgp=c(2.8,0.7,0))

# MSY weights
msy_breaks <-  c(100, 500, 1000, 2000, 10000, 20000)
msys <- cut(stats$msy_tmt, breaks=c(0, msy_breaks))
wts <- seq(1.8,4,length.out=nlevels(msys))
msy_wts <- wts[msys]

# Plot LME-scale changes: fixed vs. random effects
plot(pdiff_f ~ pdiff_r, stats, bty="n", las=1, 
     xlim=c(-80, 40), ylim=c(-250, 100), pch=21, bg="grey80", cex=msy_wts,
     xlab="Percent change in MSY\n(Random effects)",
     ylab="Percent change in MSY\n(Fixed effects)")

# Zero lines
lines(x=c(0,0), y=c(-250, 100), lty=3)
lines(x=c(-80,40), y=c(0,0), lty=3)

# One-to-one line
lines(x=c(-80,40), y=c(-80,40))

# Sample size
text(x=stats$pdiff_r, y=stats$pdiff_f, labels=stats$n, cex=0.8)

# Label outlier
out <- filter(stats, pdiff_f < -200)
text(x=out$pdiff_r, y=out$pdiff_f, labels=out$lme, pos=4, cex=0.8, offset=1)

# Off
dev.off()
graphics.off()
