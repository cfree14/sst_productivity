

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


# Which models to compare?
################################################################################

# SP model
sp <- "ramldb_v3.8_sp.Rdata"

# SP-AR1 model
sp_ar1 <- "ramldb_v3.8_sp_ar1.Rdata"

# Pella-AR1 models
pella_ar1 <- c("ramldb_v3.8_sp_ar1_pella0.55.Rdata",
               "ramldb_v3.8_sp_ar1_pella0.20.Rdata",
               "ramldb_v3.8_sp_ar1_pella0.01.Rdata")

# Pella-AR1-SST models
pella_ar1_sst <- c("ramldb_v3.8_sp_ar1_pella0.55_cobe.Rdata",
                   "ramldb_v3.8_sp_ar1_pella0.55_er.Rdata",
                   "ramldb_v3.8_sp_ar1_pella0.55_had.Rdata")

# SP-SST group models
pella_ar1_sst_group <- c("ramldb_v3.8_sp_ar1_pella0.55_cobe_order.Rdata",
                         "ramldb_v3.8_sp_ar1_pella0.55_cobe_family.Rdata",
                         "ramldb_v3.8_sp_ar1_pella0.55_cobe_fao_area.Rdata",
                         "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.Rdata",
                         "ramldb_v3.8_sp_ar1_pella0.55_cobe_method.Rdata",
                         "ramldb_v3.8_sp_ar1_pella0.55_cobe_method1.Rdata")

# SP-SST null models
pella_ar1_sst_group_null <- c("ramldb_v3.8_sp_ar1_pella0.55_cobe_lme_n1.Rdata",
                              "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme_n2.Rdata",
                              "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme_n3.Rdata")

# Merge models
models <- c(sp, sp_ar1, pella_ar1,
            pella_ar1_sst,
            pella_ar1_sst_group,
            pella_ar1_sst_group_null)

# Model name
model_names <- c("SP", "SP-AR1",
                 "Pella-45%-AR1", "Pella-40%-AR1", "Pella-37%-AR1",
                 "Pella-45%-AR1-Cobe", "Pella-45%-AR1-ERSST", "Pella-45%-AR1-HadISST",
                 "Pella-45%-AR1-Cobe-Taxa-Order", "Pella-45%-AR1-Cobe-Taxa-Family",
                 "Pella-45%-AR1-Cobe-Region-FAO area", "Pella-45%-AR1-Cobe-Region-LME",
                 "Pella-45%-AR1-Cobe-SA-Method #1", "Pella-45%-AR1-Cobe-SA-Method #2",
                 "Pella-45%-AR1-Cobe-LME-Null 1", "Pella-45%-AR1-Cobe-LME-Null 2", "Pella-45%-AR1-Cobe-LME-Null 3")
length(models)==length(model_names)


# Build table
################################################################################

# Data frame
aic_df <- data.frame(model=model_names, k=NA, lik=NA, aic=NA, convergence=NA, gradient=NA, stringsAsFactors=F)

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
  
  # Convergence check
  # If fit using Optimize()
  if("Convergence_check"%in%names(output)){
    aic_df$convergence[i] <- output$Convergence_check
    aic_df$gradient[i] <- output$max_gradient
  # If fit using nlminb()    
  }else{
    c_code <- output$convergence
    c_message <- output$message
    c_combo <- paste(c_code, c_message, sep=" - ")
    aic_df$convergence[i] <- c_combo
  }

}

# Format table
aic_final <- aic_df %>% 
  arrange(aic) %>% 
  mutate(daic=aic-min(aic)) %>% 
  select(model:aic, daic, convergence, gradient)
  

# Export table
################################################################################

# Export data
write.csv(aic_final, paste(tabledir, "Table1_model_aic_comparison.csv", sep="/"), row.names=F)

