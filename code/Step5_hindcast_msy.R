
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"

# Read model data
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
spdata <- data
rm(data, hess, input.data, model, nstocks, output, params, problem.stocks, sd, stocks, results.df, results.wide)

# Read results
data_orig <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)

# SST time series data
sst <- read.csv(paste(sstdir, "ramldb_sst_yearly_cobe.csv", sep="/"), as.is=T)


# Build MSY hindcasts
################################################################################

# Add K and MSY to results 
p <- 0.2
div <- (p+1)^((p+1)/p)
data <- data_orig %>%
  rename(k_orig = k, lme=lme_name) %>% 
  mutate(k = k_orig * tb_max,
         msy_avg = (r * k) / div) %>% 
  select(assessid, stockid, lme, sst_c_avg, betaT, betaT_lo, betaT_hi, r, k, msy_avg)
msy_avg <- sum(data$msy_avg)

# Build realized MSY data
stockids <- sort(unique(data$stockid))
assessids <- sort(unique(data$assessid))
msys <- sst %>% 
  # Reduce to stocks in analysis
  filter(assessid %in% assessids & year<=2015) %>% 
  # Add stockid/r/k/theta/SSTavg from data
  left_join(select(data, assessid, stockid, r, k, msy_avg, betaT, betaT_lo, betaT_hi, sst_c_avg), by="assessid") %>% 
  # Scale SST time series and calculate realized MSY
  mutate(sst_c_sd=sst_c-sst_c_avg, 
         msy_real=r*k/div*exp(sst_c_sd*betaT),
         msy_real_lo=pmin(r*k/div*exp(sst_c_sd*betaT_lo), r*k/div*exp(sst_c_sd*betaT_hi)),
         msy_real_hi=pmax(r*k/div*exp(sst_c_sd*betaT_lo), r*k/div*exp(sst_c_sd*betaT_hi))) %>% 
  # Rearrange columns
  select(assessid, stockid, r, k, betaT, betaT_lo, betaT_hi, sst_c_avg, 
         year, sst_c, sst_c_sd, msy_avg, msy_real, msy_real_lo, msy_real_hi)

# Inspect completeness
apply(msys, 2, function(x) sum(is.na(x)))

# Summarize realized MSY over years
msy_ts <- msys %>% 
  group_by(year) %>% 
  summarize(sst=mean(sst_c),
            msy=sum(msy_real),
            msy_lo=sum(msy_real_lo),
            msy_hi=sum(msy_real_hi))

# Plot quickly
par(mfrow=c(1,1))
plot(sst ~ year, msy_ts, type="l")
plot(msy/1E6 ~ year, msy_ts, type="l", ylim=c(15,40))
lines(msy_ts$year, msy_ts$msy_lo/1E6, lty=2)
lines(msy_ts$year, msy_ts$msy_hi/1E6, lty=2)
abline(h=msy_avg/1E6)

# Export MSY time series
write.csv(msys, paste(datadir, "msy_hindcast_time_series_using_lme_model.csv", sep="/"), row.names=F)
write.csv(msy_ts, paste(datadir, "msy_hindcast_time_series_overall_using_lme_model.csv", sep="/"), row.names=F)


