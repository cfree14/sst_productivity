
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)
library(reshape2)

# Directories
inputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output/fixed_effects"
ramldbdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/ramldb/ramldb_v3.8"
boundarydir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/stock_boundaries/data"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/code"

# Source helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# Read fixed effects output
load(paste(outputdir, "ramldb_v3.8_pella_0.20_cobe_fixed.Rdata", sep="/"))
data_f <- results.wide
rm(hess, input.data, model, nstocks, output, params, problem.stocks, sd, stocks, results.df)

# Read random effects output
load(paste(outputdir, "ramldb_v3.8_pella_0.20_cobe_random_normal.Rdata", sep="/"))
data_r <- results.wide
rm(hess, input.data, model, nstocks, output, params, problem.stocks, sd, stocks, results.df)

# Read/format RAMLDB metadata
# Note: species and common names omitted, added in LH data
ramstocks <- read.csv(paste(inputdir, "all_stocks_in_ramldb.csv", sep="/"), as.is=T)
ramstocks <- select(ramstocks, assessid, stockid, country, region, area, 
                    method, method1, method2, species, comm_name, b0_true, msy_true)

# Read/format stock boundary centroid/area info
centroids <- read.csv(paste(boundarydir, "ramldb_v3.8_stock_boundary_centroids_areas_fixed.csv", sep="/"), as.is=T)
centroids <- centroids %>% 
  mutate(lat_dd_abs=abs(lat_dd)) %>% 
  select(assessid, lme_name, fao_area, long_dd, lat_dd, lat_dd_abs, area_sqkm)

# Read SST time series
sst_ts <- read.csv(paste(sstdir, "ramldb_sst_yearly_cobe.csv", sep="/"), as.is=T)

# Read/format life history data
lhdata <- read.csv(paste(inputdir, "ramldb_v3.8_life_history_traits.csv", sep="/"), as.is=T)
lhdata <- rename(lhdata, vonb_k=k)


# Calculate SST/TB stats
################################################################################

# SST stats
sst <- sst_ts %>%
  filter(year%in%1930:2010 & !is.na(sst_c)) %>% 
  group_by(assessid) %>% 
  summarize(sst_c_avg2=mean(sst_c),
            sst_c_trend2=coef(lm(sst_c~year))[2])

# Calculate SST mean/trend and TB mean/trend
stats <- data %>%
   group_by(assessid) %>% 
    summarize(n=n(),
              yr1=min(year),
              yr2=max(year),
              sst_c_avg=mean(cobe_sst_c),
              sst_c_trend=coef(lm(cobe_sst_c~year))[2],
              tb_max=max(tb),
              tb_avg=mean(tb),
              tb_trend=coef(lm(tb~year))[2], 
              tb_sd_trend=coef(lm(tb_sd~year))[2],
              catch_avg=mean(catch), 
              bbmsy_n=sum(!is.na(bbmsy)),
              bbmsy_avg=mean(bbmsy, na.rm=T),
              ffmsy_n=sum(!is.na(ffmsy)),
              ffmsy_avg=mean(ffmsy, na.rm=T)) %>% 
  left_join(sst, by="assessid")

# Check SST from 1930-2010 against time series range
plot(sst_c_avg ~ sst_c_avg2, stats)
plot(sst_c_trend ~ sst_c_trend2, stats)


# Build data
################################################################################

# Merge model output
params <- data_f %>% 
  left_join(data_r, by="stockid")
colnames(params) <- gsub(".x", "_f", colnames(params))
colnames(params) <- gsub(".y", "_r", colnames(params))

# Merge parameter estimates with RAMLDB metadata
final <- ramstocks %>%
  right_join(lhdata, by="species") %>%
  right_join(centroids, by="assessid") %>% 
  right_join(stats, by="assessid") %>% 
  right_join(params, by="stockid")

# Inspect completeness
complete(final)

# Fix empty method
sort(unique(final$method1))
subset(final, is.na(method))
final$method[final$method1=="unknown"] <- "Unknown"
final$method2[final$method1=="unknown"] <- "Unknown"

# Calculate variable completeness
sum(!is.na(final$bbmsy_avg)) / nrow(final)
sum(!is.na(final$ffmsy_avg)) / nrow(final)

# Finfish life history
fin <- subset(final, type=="finfish")
sum(!is.na(fin$troph)) / nrow(fin)
sum(!is.na(fin$habitat)) / nrow(fin)
sum(!is.na(fin$depth_m_mid)) / nrow(fin)

# Invertebrate life history
inv <- subset(final, type=="invertebrate")
sum(!is.na(inv$m)) / nrow(inv)
sum(!is.na(inv$k)) / nrow(inv)
sum(!is.na(inv$linf_cm)) / nrow(inv)
sum(!is.na(inv$winf_kg)) / nrow(inv)
sum(!is.na(inv$lmat_cm)) / nrow(inv)
sum(!is.na(inv$tmat_yr)) / nrow(inv)
sum(!is.na(inv$tmax_yr)) / nrow(inv)
sum(!is.na(inv$troph)) / nrow(inv)
sum(!is.na(inv$habitat)) / nrow(inv)
sum(!is.na(inv$depth_m_mid)) / nrow(inv)

# Export data
################################################################################

# Export data
write.csv(final, paste(outputdir, "ramldb_v3.8_pella_0.20_cobe_merged.csv", sep="/"), row.names=F)



