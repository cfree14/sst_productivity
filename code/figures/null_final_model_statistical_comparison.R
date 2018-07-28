
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Load final model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd, results.df)
final <- results.wide

# Load primary null model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme_n1.Rdata", sep="/"))
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd, results.df)
null <- results.wide

# Merge results
data <- final %>% 
  select(stockid, betaT) %>% 
  rename(betaT_obs=betaT) %>% 
  left_join(select(null, stockid, betaT), by="stockid") %>% 
  rename(betaT_null=betaT)



