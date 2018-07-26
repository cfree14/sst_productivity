
# Clear workspace
rm(list = ls())

# Turn off scientific notation
options(scipen=999)

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(rfishbase)

# Define directories
datadir <- "~/Dropbox/Chris/Rutgers/projects/productivity/data/ramldb/ramldb_v3.8"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"

# Read keys
lme_key <- read.csv(paste(datadir, "ramldb_v38_lme_key.csv", sep="/"), as.is=T)
taxa_key_ram <- read.csv(paste(datadir, "ramldb_v38_taxonomy.csv", sep="/"), as.is=T)
stock_key <- read.csv(paste(datadir, "ramldb_v38_stock.csv", sep="/"), as.is=T)
assessment_key <- read.csv(paste(datadir, "ramldb_v38_assessment.csv", sep="/"), as.is=T)
method_key <- read.csv(paste(datadir, "ramldb_v38_assessmethod.csv", sep="/"), as.is=T)
assessor_key <- read.csv(paste(datadir, "ramldb_v38_assessor.csv", sep="/"), as.is=T)
area_key <- read.csv(paste(datadir, "ramldb_v38_area.csv", sep="/"), as.is=T)

# Read data
values <- read.csv(paste(datadir, "ramldb_v38_timeseries_values_views.csv", sep="/"), as.is=T)
units <- read.csv(paste(datadir, "ramldb_v38_timeseries_units_views.csv", sep="/"), as.is=T)
bioparams_vals <- read.csv(paste(datadir, "ramldb_v38_bioparams_values_views.csv", sep="/"), as.is=T)
bioparams_units <- read.csv(paste(datadir, "ramldb_v38_bioparams_units_views.csv", sep="/"), as.is=T)

# FishBase/SeaLifeBase taxa keys
taxa_key_fb <- load_taxa()
taxa_key_slb <- sealifebase


# Format taxa info
################################################################################

# Merge FishBase/SeaLifeBase taxa keys
taxa_key <- taxa_key_fb %>% 
    bind_rows(taxa_key_slb) %>% 
    setNames(tolower(names(.))) %>% 
    mutate(sciname=paste(genus, species)) %>%
    select(class, order, family, genus, sciname) %>% 
    rename(species=sciname) %>% 
    unique()


# Format time series info
################################################################################

# Catch, abundance, and recruitment
# Catch: TL=total landings, TC=total catch
# Abundance: TB=total biomass, SSB=spawning stock biomass, TN=total number
# Recruitment: R=recruitment

# Time series units
ts_units <- units %>% 
  select(-c(stockid, stocklong, F, ER)) %>% 
  rename(tb_units=TB, ssb_units=SSB, tn_units=TN, r_units=R, 
         tc_units=TC, tl_units=TL)
ts_units[ts_units==""] <- NA

# Time series lengths and averages
ts_stats <- values %>%
    group_by(assessid) %>%
    summarize(tl_nyrs = sum(!is.na(TL)),
              tc_nyrs = sum(!is.na(TC)),
              tb_nyrs = sum(!is.na(TB)),
              ssb_nyrs = sum(!is.na(SSB)),
              tn_nyrs = sum(!is.na(TN)),
              r_nyrs = sum(!is.na(R)),
              tl_avg = mean(TL, na.rm=T),
              tc_avg = mean(TC, na.rm=T),
              tb_avg = mean(TB, na.rm=T),
              ssb_avg = mean(SSB, na.rm=T),
              tn_avg = mean(TN, na.rm=T),
              r_avg = mean(R, na.rm=T)) %>%
    mutate_all(funs(replace(., is.na(.), NA)))

# Merge time series info
ts_info <- ts_units %>% 
  left_join(ts_stats, by="assessid") %>% 
  mutate(catch_type=ifelse(tc_nyrs!=0 | tl_nyrs!=0, 
                           ifelse(tc_nyrs>=tl_nyrs, "tc", "tl"), "none"),
         catch_nyrs=ifelse(catch_type=="tc", tc_nyrs, tl_nyrs),
         catch_avg=ifelse(catch_type=="tc", tc_avg, tl_avg), 
         catch_units=ifelse(catch_type=="tc", tc_units, tl_units)) %>%  
  select(assessid,
         tl_nyrs, tl_avg, tl_units,
         tc_nyrs, tc_avg, tc_units,
         catch_type, catch_nyrs, catch_avg, catch_units,
         r_nyrs, r_avg, r_units, 
         ssb_nyrs, ssb_avg, ssb_units,
         tb_nyrs, tb_avg, tb_units,
         tn_nyrs, tn_avg, tn_units)

# Biological parameters data
bioparams <- bioparams_vals %>%
  select(assessid, MSY, SSB0) %>% 
  rename(msy_true=MSY, b0_true=SSB0) %>% 
  left_join(select(bioparams_units, assessid, MSY, SSB0), by="assessid") %>% 
  rename(msy_true_units=MSY, b0_true_units=SSB0)
bioparams[bioparams==""] <- NA


# Format assessment method data
################################################################################

# Remove the unknown BDM and SCAA models because they cause field duplication
# when merging the broad and long method names to the assessment key
# This introduces more generic "Unknowns" when more info is known but I'll correct these later
method_key1 <- method_key %>% 
  filter(!methodlong%in%c("Unknown biomass dynamics model", "Unknown statistical catch-at-age model"))

# Build merged data
################################################################################

# Build merged dataset
data <- assessment_key %>% 
  # Format assessment info
  select(assessid, stockid, assessorid, stocklong, assessyear, assessmethod) %>% 
  rename(assessor=assessorid, years=assessyear, method1=assessmethod) %>% 
  # Add more assessment method info 
  left_join(method_key1, by=c("method1"="methodshort")) %>%
  rename(method=category, method2=methodlong) %>% 
  mutate(method1=revalue(method1, c("Unknown"="unknown")),
       year=as.numeric(substr(years, 6, 9))) %>% 
  # Add accessor info: country and mgmt agency
  left_join(select(assessor_key, -assessorfull), by=c("assessor"="assessorid")) %>% 
  mutate(country=revalue(country, c("Multinational"="multinational"))) %>% 
  # Add stock info: region, area id, scientific name
  left_join(select(stock_key, -c(stocklong, tsn, inmyersdb, myersstockid, commonname)), by="stockid") %>% 
  rename(species_orig=scientificname) %>%
  mutate(species_orig=trimws(species_orig),
         region=revalue(region, c("Europe non EU"="Europe (non-EU)",
                                  "European Union"="Europe (EU)"))) %>% 
  # Add common name
  left_join(select(taxa_key_ram, scientificname, commonname1), by=c("species_orig"="scientificname")) %>% 
  rename(comm_name=commonname1) %>%       
  # Add area info: area name
  left_join(select(area_key, -c(country, alternateareaname, areatype, areacode)), by="areaid") %>%
  rename(area=areaname) %>% 
  # Format species names and add taxanomic info
  mutate(species=trimws(species_orig),
         species=revalue(species, c("Cervimunida Johni" = "Cervimunida johni",
                                    "Chrysophrys auratus" = "Pagrus auratus",
                                    # "Clupea bentincki" = "Strangomera bentincki",
                                    "Clupea pallasii" = "Clupea pallasii pallasii",
                                    "Epinephelus niveatus" = "Hyporthodus niveatus",
                                    "Epinephelus flavolimbatus" = "Hyporthodus flavolimbatus",
                                    "Etrumeus teres" = "Etrumeus sadina",
                                    "Loligo pealeii" = "Doryteuthis pealeii",
                                    # "Loligo reynaudii" = "Loligo vulgaris reynaudii",
                                    "Merluccius gayi" = "Merluccius gayi gayi",
                                    "Mullus barbatus" = "Mullus barbatus barbatus",
                                    "Neoplatycephalus richardsoni" = "Platycephalus richardsoni",
                                    "Psetta maxima" = "Scophthalmus maximus",
                                    "Reinhardtius stomias" = "Atheresthes stomias",
                                    "Sardinops melanostictus" = "Sardinops sagax",
                                    "Scomber australacius"="Scomber australasicus", 
                                    "Solea vulgaris" = "Solea solea",
                                    "Sprattus fuengensis" = "Sprattus fuegensis",
                                    "Tetrapturus albidus" = "Kajikia albida"))) %>% 
  left_join(taxa_key, by="species") %>% 
  # Add MSY/SSB0 true
  left_join(bioparams, by="assessid") %>% 
  # Add time series info
  left_join(ts_info, by="assessid") %>% 
  # Reorder columns
  select(assessid, stockid, mgmt, assessor, years, year, method, method1, method2,
         country, region, area, areaid, stocklong, 
         class, order, family, genus, species, species_orig, comm_name, everything())

# Fix problem species
prob_spp <- sort(unique(data$species_orig[is.na(data$class)]))
for(i in 1:length(prob_spp)){
  spp <- prob_spp[i]
  genus <- strsplit(spp, " ")[[1]][1]
  data$genus[data$species==spp] <- genus
  data$family[data$species==spp] <- unique(taxa_key$family[taxa_key$genus==genus])
  data$order[data$species==spp] <- unique(taxa_key$class[taxa_key$genus==genus])
  data$class[data$species==spp] <- unique(taxa_key$order[taxa_key$genus==genus])
}
data$comm_name[data$species=="Atheresthes stomias"] <- "Arrowtooth flounder"
data$comm_name[data$species=="Scomber australasicus"] <- "Blue mackerel"
data$comm_name[data$species=="Penaeus semisulcatus"] <- "Green tiger prawn"

# Inspect completness
apply(data, 2, function(x) sum(is.na(x)))


# Mark candidate stocks for analysis
################################################################################

# SURPLUS PRODUCTION MODEL CONDITIONS
# 1. Exclude salmon stocks
# 2. Exclude assessments with <20 years of catch (TC>TL) and abundance (TB) data in MT

# Salmon stocks
table(data$region)
salmon.regions <- c("US West Coast (Pacific Salmon)", "US Alaska (Pacific Salmon)",
                    "Canada West Coast (Pacific Salmon)", "Russia Japan (Pacific Salmon)")

# Assessment methods
table(data$method)

# Examine surplus production dataset sample size
nrow(subset(data, !(region%in%salmon.regions)))
nrow(subset(data, !(region%in%salmon.regions)
            & catch_units=="MT" & tb_units=="MT"))
nrow(subset(data, !(region%in%salmon.regions)
            & catch_units=="MT" & tb_units=="MT" & catch_nyrs>=20 & tb_nyrs >=20))

# Mark stocks for use
sp.model.stocks <- subset(data, !(region%in%salmon.regions)
                          & catch_nyrs>=20 & tb_nyrs >=20 & catch_units=="MT" & tb_units=="MT")
data$sp_model[data$assessid%in%sp.model.stocks$assessid] <- "yes"

# Final sample size
sum(data$sp_model=="yes", na.rm=T)


# Export data
################################################################################

# Export data
write.csv(data, paste(outputdir, "all_stocks_in_ramldb.csv", sep="/"), row.names=F)

