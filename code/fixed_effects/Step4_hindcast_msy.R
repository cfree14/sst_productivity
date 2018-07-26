
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output/fixed_effects"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures/fixed_effects"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"

# Read model data
# It doesn't matter whether you read fixed or random -- same data
load(paste(datadir, "ramldb_v3.8_pella_0.20_cobe_fixed.Rdata", sep="/"))
spdata <- data
rm(data, hess, input.data, model, nstocks, output, params, problem.stocks, sd, stocks, results.df, results.wide)

# Read results
data_orig <- read.csv(paste(datadir, "ramldb_v3.8_pella_0.20_cobe_merged.csv", sep="/"), as.is=T)

# SST time series data
sst <- read.csv(paste(sstdir, "ramldb_sst_yearly_cobe.csv", sep="/"), as.is=T)


# Build MSY hindcasts
################################################################################

# Add K and MSY to results 
p <- 0.2
div <- (p+1)^((p+1)/p)
data <- data_orig %>%
  rename(k_orig_r = k_r,
         k_orig_f = k_f, 
         lme=lme_name) %>% 
  mutate(k_r = k_orig_r * tb_max,
         k_f = k_orig_f * tb_max,
         msy_avg_r = (r_r * k_r) / div,
         msy_avg_f = (r_f * k_f) / div) %>% 
  select(assessid, stockid, lme, sst_c_avg, betaT_r, r_r, k_r, msy_avg_r, betaT_f, r_f, k_f, msy_avg_f)
msy_avg_r <- sum(data$msy_avg_r)
msy_avg_f <- sum(data$msy_avg_f)

# Build realized MSY data
stockids <- sort(unique(data$stockid))
assessids <- sort(unique(data$assessid))
msys <- sst %>% 
  # Reduce to stocks in analysis
  filter(assessid %in% assessids & year<=2015) %>% 
  # Add stockid/r/k/theta/SSTavg from data
  left_join(select(data, assessid, stockid, sst_c_avg,
                   r_r, k_r, msy_avg_r, betaT_r,
                   r_f, k_f, msy_avg_f, betaT_f), by="assessid") %>% 
  # Scale SST time series and calculate realized MSY
  mutate(sst_c_sd=sst_c-sst_c_avg, 
         msy_real_r=r_r*k_r/div*exp(sst_c_sd*betaT_r),
         msy_real_f=r_f*k_f/div*exp(sst_c_sd*betaT_f)) %>% 
  # Rearrange columns
  select(assessid, stockid, sst_c_avg, 
         r_r, k_r, betaT_r, r_f, k_f, betaT_f, 
         year, sst_c, sst_c_sd, 
         msy_avg_r, msy_avg_f, msy_real_r, msy_real_f)

# Inspect completeness
apply(msys, 2, function(x) sum(is.na(x)))

# Summarize realized MSY over years
msy_ts <- msys %>% 
  group_by(year) %>% 
  summarize(sst=mean(sst_c),
            msy_r=sum(msy_real_r),
            msy_f=sum(msy_real_f))

# Plot quickly
par(mfrow=c(1,1))
plot(sst ~ year, msy_ts, type="l")
plot(msy_r/1E6 ~ year, msy_ts, type="l", ylim=c(30,50), col="blue")
lines(msy_ts$year, msy_ts$msy_f/1E6, col="red")

# Export MSY time series
write.csv(msys, paste(datadir, "msy_hindcast_time_series.csv", sep="/"), row.names=F)
write.csv(msy_ts, paste(datadir, "msy_hindcast_time_series_overall.csv", sep="/"), row.names=F)


