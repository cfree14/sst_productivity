
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/tables"

# Read data
data_orig <- read.csv(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.csv", sep="/"), as.is=T)

# Read MSY hindcasts
msy_ts <- read.csv(paste(datadir, "msy_hindcast_time_series_using_lme_model.csv", sep="/"), as.is=T)


# Stats for manuscript
################################################################################

# Sample size info
n_distinct(data$stockid)
n_distinct(data$species)
n_distinct(data$lme_name)


# Build data
################################################################################

# Years
yrs1 <- 1930:1939
yrs2 <- 2001:2010

# MSY stats
msy <- msy_ts %>% 
  group_by(stockid) %>% 
  summarise(msy1=mean(msy_real[year%in%yrs1]),
            msy2=mean(msy_real[year%in%yrs2]),
            msy_pdiff=(msy2-msy1)/msy1*100)

# Build data
data <- data_orig %>%
  mutate(comm_name=revalue(comm_name, c("common European sole"="Common sole",
                                        "Common seabream"="Red porgy",
                                        "Hake"="European hake",
                                        "Herring"="Atlantic herring",
                                        "Walleye pollock"="Alaska pollock", 
                                        "Pollock"="Saithe")),
         spp_name=paste0(comm_name, " (", species, ")"),
         betaT=round(betaT,2),
         sst_c_avg=round(sst_c_avg,1),
         sst_c_trend_10yr=round(sst_c_trend*10,2)) %>% 
  select(lme_name, stockid, spp_name, area, sst_c_avg, sst_c_trend_10yr, betaT, betaT_inf) %>% 
  left_join(msy, by="stockid") %>% 
  mutate(msy_pdiff=round(msy_pdiff,1))
  arrange(desc(betaT))
  
# Export data
################################################################################

# Export table
write.csv(data, paste(tabledir, "TableS3_all_stocks_msy_change_stats.csv", sep="/"), row.names=F)
