

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/tables"


# Which models to compare?
################################################################################

# PT model
sp <- "ramldb_v3.8_pella_0.20.Rdata"

# PT SST-linked (fixed effecst)
sp_sst_fixed <- "ramldb_v3.8_pella_0.20_cobe_fixed.Rdata"

# PT SST-linked (random effects - normal)
sp_sst_random_norm <- "ramldb_v3.8_pella_0.20_cobe_random_normal.Rdata"

# Merge models
models <- c(sp, sp_sst_fixed, sp_sst_random_norm)

# Model name
model_names <- c("SP", 
                 "SP-SST-Fixed", 
                 "SP-SST-Random")
length(models)==length(model_names)


# Build table
################################################################################

# Data frame
aic_df <- data.frame(model=model_names, k=NA, lik=NA, aic=NA, stringsAsFactors=F)

# Loop through models and calculate/record AIC value
for(i in 1:length(models)){
  
  # Load data
  load(paste(datadir, models[i], sep="/"))
  
  # Calculate/record AIC value
  k <- length(output[["par"]])
  lik <- output[["objective"]]
  aic_val <- TMBhelper::TMBAIC(output)
  aic_df$k[i] <- k
  aic_df$lik[i] <- lik
  aic_df$aic[i] <- aic_val
  
}

# Format table
aic_final <- aic_df %>% 
  arrange(aic) %>% 
  mutate(daic=aic-min(aic))
  

# Export table
################################################################################

# Export data
write.csv(aic_final, paste(tabledir, "Table1_model_aic_comparison.csv", sep="/"), row.names=F)

