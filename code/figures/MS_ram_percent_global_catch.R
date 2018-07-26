

# Clear workspace
rm(list = ls())

# Turn off scientific notation
options(scipen=999)

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/input"
datadir1 <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
ramdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/ramldb/ramldb_v3.8"

# RAMLDB stock and time series data
ram.stocks <- read.csv(paste(datadir, "all_stocks_in_ramldb.csv", sep="/"), as.is=T)
ram.ts.vals <- read.csv(paste(ramdir, "ramldb_v38_timeseries_values_views.csv", sep="/"), as.is=T)
ram.ts.units <- read.csv(paste(ramdir, "ramldb_v38_timeseries_units_views.csv", sep="/"), as.is=T)
colnames(ram.ts.vals) <- tolower(colnames(ram.ts.vals))


# In whole RAM database
################################################################################

# Build data
data <- ram.ts.vals %>% 
  # Filter to 2000
  filter(year==2000) %>% 
  select(stockid, year, tc, tl) %>% 
  # Add units
  left_join(select(ram.ts.units, stockid, TC, TL), by="stockid") %>% 
  rename(tc_units=TC, tl_units=TL) %>% 
  # Overwrite E00 units
  mutate(tc=ifelse(tc_units=="E00", NA, tc),
         tl=ifelse(tl_units=="E00", NA, tl)) %>% 
  # Calculate "final" catch
  mutate(catch=pmax(tc, tl, na.rm=T)) %>% 
  filter(!is.na(catch))

# Total catch in RAM in 2000
catch_tot <- sum(data$catch) / 1E6
catch_tot / 86


# In our dataset
################################################################################

# Read data
load(paste(datadir1, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))

catch_us <- sum(data$catch[data$year==2000]) / 1E6
catch_us / 86






