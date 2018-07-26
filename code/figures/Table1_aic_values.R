

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
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/tables"


# Which models to compare?
################################################################################

# SP model
sp <- "ramldb_v3.8_sp.Rdata"

# SP-SST models
sp_sst <- c("ramldb_v3.8_spsst_cobe.Rdata",
            "ramldb_v3.8_spsst_er.Rdata",
            "ramldb_v3.8_spsst_had.Rdata")

# SP-SST PT models
sp_sst_pt <- c("ramldb_v3.8_spsst_pella_0.01_cobe.Rdata",
               "ramldb_v3.8_spsst_pella_0.20_cobe.Rdata",
               "ramldb_v3.8_spsst_pella_0.55_cobe.Rdata")

# SP-SST group models
sp_sst_group <- c("ramldb_v3.8_spsst_pella_cobe_order.Rdata",
                  "ramldb_v3.8_spsst_pella_cobe_family.Rdata",
                  "ramldb_v3.8_spsst_pella_cobe_fao_area.Rdata",
                  "ramldb_v3.8_spsst_pella_cobe_lme.Rdata",
                  "ramldb_v3.8_spsst_pella_cobe_method.Rdata",
                  "ramldb_v3.8_spsst_pella_cobe_method1.Rdata")

# SP-SST null models
sp_sst_null <- c("ramldb_v3.8_spsst_pella_cobe_method1_n1.Rdata", 
                 "ramldb_v3.8_spsst_pella_cobe_method1_n2.Rdata",
                 "ramldb_v3.8_spsst_pella_cobe_method1_n3.Rdata")

# Merge models
models <- c(sp, sp_sst, sp_sst_pt, sp_sst_group, sp_sst_null)

# Model name
model_names <- c("SP", 
                 "SP-SST-Cobe", "SP-SST-ERSST", "SP-SST-HadISST", 
                 "SP-SST-Pella 37%", "SP-SST-Pella 40%", "SP-SST-Pella 45%",
                 "SP-SST-Pella-Taxa-Order", "SP-SST-Pella-Taxa-Family",
                 "SP-SST-Pella-Region-FAO area", "SP-SST-Pella-Region-LME",
                 "SP-SST-Pella-SA-Method #1", "SP-SST-Pella-SA-Method #2",
                 "SP-SST-Pella-Region-Null 1", "SP-SST-Pella-Region-Null 2", "SP-SST-Pella-Region-Null 3")
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

