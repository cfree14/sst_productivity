
# Clear workspace
rm(list = ls())

# Turn off scientific notation
options(scipen=999)

# Read data
################################################################################

# Packages
library(TMB)
library(plyr)
library(dplyr)
library(ggplot2)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output/simulations"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/code"
  
# Read model output and input
load(paste(datadir, "simulations_self.Rdata", sep="/"))

# Souce "fit model" code
source(paste(codedir, "fit_sp_sst_model.R", sep="/"))


# Simulation test
################################################################################

# Scenarios
scenarios <- unique(sim_key$scenario)

# Loop through scenarios and apply model to stocks
for(i in 1:length(scenarios)){
  
  # Subset data
  scenario_i <- scenarios[i]
  stocks <- subset(sim_key, scenario==scenario_i)
  sdata <- subset(sim_data, scenario==scenario_i)
  
  # Fit SP-SST model
  spfit <- fit_sp_sst(id_t=sdata$stockid, 
                      p_t=sdata$pt_sim,
                      b_t=sdata$bt_sim,
                      sst_t=sdata$sst_obs)
  
  # Save results
    
  
  
}


















