
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
  select(stockid, betaT, betaT_inf) %>% 
  rename(betaT_obs=betaT, betaT_obs_inf=betaT_inf) %>% 
  left_join(select(null, stockid, betaT, betaT_inf), by="stockid") %>% 
  rename(betaT_null=betaT, betaT_null_inf=betaT_inf)

# Params
n <- nrow(data)
nsig_obs <- sum(data$betaT_obs_inf != "none")
nsig_null <- sum(data$betaT_null_inf != "none")

# Binomial exact test
# Compare proportion significant in final vs. proportion significant in null
binom.test(x=nsig_obs, n=n, p=nsig_null/n)





