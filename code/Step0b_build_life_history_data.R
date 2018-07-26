
# Clear workspace
rm(list = ls())

# Turn off scientific notation
options(scipen=999)

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)
library(stringr)
library(FishLife)
library(rfishbase)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/input"
ramdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/ramldb/ramldb_v3.8"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"

# RAMLDB stocks
ramstocks <- read.csv(paste(ramdir, "all_stocks_in_ramldb.csv", sep="/"), as.is=T)


# Taxa key
################################################################################

# FishBase/SeaLifeBase taxa keys
taxa_key_fb <- load_taxa()
taxa_key_slb <- sealifebase

# Build taxa key
taxa_key <- taxa_key_fb %>% 
  bind_rows(taxa_key_slb) %>% 
  setNames(tolower(names(.))) %>% 
  mutate(sciname=paste(genus, species)) %>%
  select(class, order, family, genus, sciname) %>% 
  rename(species=sciname) %>% 
  unique()


# Build taxanomic info
################################################################################

# Build data
data <- ramstocks %>%
  # If unique species
  select(species) %>% 
  unique() %>% 
  # Add taxanomic info
  left_join(taxa_key, by="species") %>% 
  # Add taxanomic level
  mutate(taxa_level=ifelse(grepl("spp", species), "genus", "species")) %>% 
  select(everything(), species, taxa_level)

# Fill in taxanomic gaps (taxa identified to genus only)
genus_spp <- data$species[data$taxa_level=="genus"]
for(i in 1:length(genus_spp)){
  spp <- genus_spp[i]
  genus1 <- word(spp, 1)
  taxa_info <- taxa_key %>% 
    select(-species) %>% 
    filter(genus==genus1) %>% 
    unique()
  data[data$species==spp,2:5] <- taxa_info
}

# Fix "Loligo reynaudii" (which isn't in SeaLifeBase)
data[data$species=="Loligo reynaudii",2:5] <- taxa_key[taxa_key$species=="Loligo vulgaris",1:4]

# Add finfish/shellfish column and r
finfish_classes <- c("Actinopterygii", "Elasmobranchii")
data <- data %>% 
  mutate(type=ifelse(class%in%finfish_classes, "finfish", "shellfish")) %>% 
  select(type, class, order, family, genus, species, taxa_level)

# Are species names unique?
anyDuplicated(data$species)

# Inspect completness
apply(data, 2, function(x) sum(is.na(x)))


# Get finfish life history from FishLife
################################################################################

# Notes on FishLife
# FishLife reports predictions in log-space, except temperature
# FishLife reports lengths in cm, weights in g, and temps in C

# Setup container
finfish <- data$species[data$type=="finfish"]
fin_lh <- data.frame(species=finfish, linf=NA, k=NA, winf=NA, tmax=NA, tm=NA, 
                     m=NA, lm=NA, temp=NA, stringsAsFactors=F)

# Loop through species
for(i in 1:nrow(fin_lh)){
  
  # Get spp info
  sciname <- fin_lh$species[i]
  genus <- word(sciname, 1)
  nwords_in_spp <- length(strsplit(sciname, " ")[[1]])
  species <- word(sciname, start=2, end=nwords_in_spp)
  species <- ifelse(species=="spp", "predictive", species)
  
  # Get and plot life history info
  spp_info <- Plot_taxa(Search_species(Genus=genus, Species=species)$match_taxonomy)
  spp_lh_vals <- exp(spp_info[[1]]$Mean_pred)
  fin_lh[i,2:ncol(fin_lh)] <- spp_lh_vals
  
}

# Format life history data
fin_lh1 <- fin_lh %>% 
  mutate(temp=log(temp),
         winf=winf/1000) %>% 
  rename(linf_cm=linf, winf_kg=winf, tmax_yr=tmax, tmat_yr=tm, lmat_cm=lm, temp_c=temp)


# Get FB/SLB life history info
################################################################################

# Finfish
##################################################

# Get FB data
finfish <- data$species[data$type=="finfish"]
base_fb <- species(finfish) # for habitat and depth preferences - one per species
troph_fb <- ecology(finfish) # for trophic level - one per species

# Format FishBase info
base_fb1 <- base_fb %>% 
  select(sciname, DemersPelag, 
         DepthRangeShallow, DepthRangeDeep, 
         DepthRangeComShallow, DepthRangeComDeep) %>% 
  rename(habitat=DemersPelag,
         depth_m_min=DepthRangeShallow,
         depth_m_max=DepthRangeDeep,
         comm_depth_m_min=DepthRangeComShallow,
         comm_depth_m_max=DepthRangeComDeep) %>% 
  mutate(depth_m_mid=(depth_m_max+depth_m_min)/2) %>% 
  select(sciname, habitat, depth_m_mid, everything())
troph_fb1 <- troph_fb %>%
  select(sciname, FoodTroph, FeedingType) %>%
  rename(troph=FoodTroph, feeding=FeedingType)

# Add to finfish data
fin_lh2 <- fin_lh1 %>% 
  left_join(base_fb1, by=c("species"="sciname")) %>% 
  left_join(troph_fb1, by=c("species"="sciname"))


# Invertebrates
##################################################

# Set to SeaLifeBase
options(FISHBASE_API = "https://fishbase.ropensci.org/sealifebase")

# Get SLB data
inverts <- data$species[data$type=="shellfish"]
base_slb <- species(inverts) # for habitat and depth preferences - one per species
troph_slb <- ecology(inverts) # for trophic level - one per species
growth_slb <- popgrowth(inverts)
stocks_slb <- stocks(inverts)
maxs_slb <- popchar(inverts) # tmax=tmax_yr, Lmax=lmax_cm, Wmax=wmax_g
maturity_slb <- maturity(inverts) # tm=tmat_yr, Lm=lmat_cm

# Format FishBase info
base_slb1 <- base_slb %>% 
  select(sciname, DemersPelag, Length, Weight, LongevityWild,
         DepthRangeShallow, DepthRangeDeep, DepthRangeComShallow, DepthRangeComDeep) %>% 
  rename(habitat=DemersPelag, lmax_cm=Length, wmax_kg=Weight, tmax_yr=LongevityWild,
         depth_m_min=DepthRangeShallow, depth_m_max=DepthRangeDeep,
         comm_depth_m_min=DepthRangeComShallow, comm_depth_m_max=DepthRangeComDeep) %>% 
  mutate(wmax_kg = wmax_kg / 1000,
         depth_m_mid=(depth_m_max+depth_m_min)/2)
troph_slb1 <- troph_slb %>%
  select(sciname, FoodTroph, FeedingType) %>%
  rename(troph=FoodTroph, feeding=FeedingType)
maxs_slb1 <- maxs_slb %>% 
  group_by(sciname) %>%
  summarize(lmax_cm=mean(Lmax[Type=="TL"], na.rm=T),
            wmax_kg=mean(Wmax, na.rm=T)/1000,
            tmax_yr=mean(tmax, na.rm=T))
growth_slb1 <- growth_slb %>% 
  group_by(sciname) %>% 
  summarize(linf_cm=mean(Loo, na.rm=T),
            winf_kg=mean(Winfinity, na.rm=T)/1000,
            k=mean(K, na.rm=T),
            m=mean(M, na.rm=T))
# tmax_yr=mean(tmax, na.rm=T), 
# tmat_yr=mean(tm, na.rm=T), 
# lmat_cm=mean(Lm[TypeLm=="TL"], na.rm=T)) 
maturity_slb1 <- maturity_slb %>% 
  group_by(sciname) %>% 
  summarize(tmat_yr=mean(tm, na.rm=T),
            lmat_cm=mean(Lm[Type1=="TL"], na.rm=T)) 

# Merge FishBase life history data
inv_lh <- base_slb1 %>% 
  left_join(growth_slb1, by="sciname") %>% 
  left_join(troph_slb1, by="sciname") %>% 
  left_join(maturity_slb1, by="sciname") %>% 
  mutate_all(funs(replace(., is.nan(.), NA))) %>% 
  rename(species=sciname) %>% 
  select(species, linf_cm, k, winf_kg, tmax_yr, tmat_yr, m, lmat_cm, habitat, troph, feeding)


# Merge finfish and invertebrate data
##################################################

# Merge finfish/invert data
all_lh <- rbind.fill(fin_lh2, inv_lh)

# Merge taxanomic and life history data
data1 <- data %>% 
  left_join(all_lh, by="species") %>% 
  mutate(type=revalue(type, c("shellfish"="invertebrate")))

# Completeness
complete(data1)

# Export data
################################################################################

# Export data
save(base_fb, troph_fb,
     base_slb, troph_slb, growth_slb, stocks_slb, maxs_slb, maturity_slb,
     file=paste(datadir, "ramldb_v3.8_species_fb_slb_data.Rdata", sep="/"))
write.csv(data1, paste(datadir, "ramldb_v3.8_life_history_traits.csv", sep="/"), row.names=F)

