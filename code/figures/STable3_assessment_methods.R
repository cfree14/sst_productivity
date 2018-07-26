

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/tables"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)


# Build table
################################################################################

# Revalue countries
table(data$country)
table(data$region[data$country=="multinational"])
data$country1 <- data$country
data$country1[data$region%in%c("Europe (EU)", "Europe (non-EU)")] <- "Europe"
data$country1[grepl("Ocean", data$region)] <- "Tuna-RFMO"
data$country1[data$region=="Mediterranean-Black Sea"] <- "Tuna-RFMO"
data$country1[data$country=="multinational" & data$region=="South America"] <- "Chile"
data$country1[data$country=="multinational" & data$region=="US West Coast"] <- "Tuna-RFMO"
data$country1[data$country=="multinational" & data$region=="Canada East Coast"] <- "Canada"
data$country1[data$region=="West Africa"] <- "West Africa"
table(data$country1)

# Build table
methods <- data %>% 
  group_by(method, method1, method2) %>% 
  summarize(n=n(),
            regions=paste(sort(unique(country1)), collapse=", ")) %>% 
  arrange(method, desc(n)) %>% 
  rename(category=method, acronym=method1, name=method2)


# Export table
################################################################################

# Export table
write.csv(methods, paste(tabledir, "STable3_stock_assessment_methods.csv", sep="/"), row.names=F)


