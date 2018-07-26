
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

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)


# Stats for manuscript
################################################################################

# Sample size info
n_distinct(data$stockid)
n_distinct(data$species)
n_distinct(data$lme_name)


# Plot data
################################################################################

# Build data
data1 <- data %>%
  mutate(comm_name=revalue(comm_name, c("common European sole"="Common sole",
                                        "Common seabream"="Red porgy",
                                        "Hake"="European hake",
                                        "Hawaiian morwong"="Tarakihi",
                                        "Herring"="Atlantic herring",
                                        "Walleye pollock"="Alaska pollock", 
                                        "Pollock"="Saithe")),
         spp_name=paste0(comm_name, " (", species, ")"),
         betaT=round(betaT,2)) %>% 
  filter(betaT_inf!="none") %>% 
  select(stockid, spp_name, area, betaT) %>% 
  arrange(desc(betaT))
  
# Export data
################################################################################

# Export table
write.csv(data1, paste(tabledir, "STable5_significant_stocks.csv", sep="/"), row.names=F)
