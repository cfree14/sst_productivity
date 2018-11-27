

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(TMB)
library(plyr)
library(dplyr)
library(devtools)
library(reshape2)
library(freeR)
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Directories
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Read data
load(file.path(outputdir, "ramldb_v3.8_sp.Rdata"))
resids_base <- model$report()

# Read data
load(file.path(outputdir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata"))
resids_final <- model$report()

# Build data
################################################################################

# Build data
df <- data %>% 
  ungroup() %>% 
  select(stockid, year) %>% 
  mutate(resids_base=resids_base[[1]],
         resids_final=resids_final[[1]])

# Setup plot
pdf(file.path(plotdir, "base_vs_final_model_resid_check.pdf"), width=8.5, height=11)
par(mfrow=c(6,4), oma=c(2,2,2,2))

# Loop through and plot
for(i in 1:nstocks){
  stock <- stocks[i]
  sdata <- filter(df, stockid==stock)
  plot(resids_base ~ year, sdata, type="l", bty="n", xlab="", ylab="SP residuals", las=1, main=stock)
  lines(resids_final ~ year, sdata, col="red", lty=3)
  abline(h=0, lty=3)
}

dev.off()





