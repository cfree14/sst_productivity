
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"

# Read data
load(paste(datadir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))
data_orig <- data

# Which SST dataset and shape parameter?
# SST choices: cobe, er, had
# Shape parameters: 0.55 (45%), 0.20 (40%), 0.01 (37%)
p <- 0.2
sst <- "cobe"
sst_col <- paste0(sst, "_sst_c_sd")


# Helper functions
################################################################################

# Function to fit surplus production model
fit_pella_nls <- function(sp, tb, p, r_start=NA, k_start=NA){
  if(is.na(r_start)){r_start <- log(0.4)}
  if(is.na(k_start)){k_start <- log(max(tb) * 1.5)}
  spfit <- try(nls(sp ~ exp(r)/p*tb*(1-(tb/exp(k))^p),
                   start=list(r=r_start, k=k_start)))
  return(spfit)
}

# Function to fit SST-linked surplus production model
fit_pella_sst_nls <- function(sp, tb, sst, p, r_start=NA, k_start=NA){
  betaT_start <- 0
  if(is.na(r_start)){r_start <- log(0.4)}
  if(is.na(k_start)){k_start <- log(max(tb) * 1.5)}
  spfit <- try(nls(sp ~ exp(r)/p*tb*(1-(tb/exp(k))^p)*exp(sst*betaT),
                   start=list(r=r_start, k=k_start, betaT=betaT_start)))
  return(spfit)
}


# Fit surplus production model
################################################################################

# Data frame to record SP fits
spfits <- data %>%
    group_by(assessid, stockid) %>%
    summarize(tb_max=max(tb)) %>%
    mutate(r=NA, k=NA, 
           r2=NA, k2=NA, betaT=NA)

# Loop and fit
for(i in 1:nrow(spfits)){
  
  # Subset data
  stock <- spfits$stockid[i]
  sdata <- subset(data, stockid==stock)
  print(paste(i, stock))
  # plot(sp_sd ~ tb_sd, sdata)
  
  # Fit SP model
  spfit <- fit_pella_nls(sp=sdata$sp_sd, tb=sdata$tb_sd, p=p)
  if(!(inherits(spfit, "try-error"))){
    r <- exp(coef(spfit)["r"])
    k <- exp(coef(spfit)["k"])
    spfits$r[spfits$stockid==stock] <- r
    spfits$k[spfits$stockid==stock] <- k
    # curve(r/p*x*(1-(x/k)^p), from=0, to=1, lty=1, lwd=0.9, col="black")
  }
  
  # Fit SP model again
  r <- spfits$r[spfits$stockid==stock]
  k <- spfits$r[spfits$stockid==stock]
  spfit <- fit_pella_nls(sp=sdata$sp_sd, tb=sdata$tb_sd, p=p, r_start=r, k_start=k)
  if(!(inherits(spfit, "try-error"))){
    r <- exp(coef(spfit)["r"])
    k <- exp(coef(spfit)["k"])
    spfits$r[spfits$stockid==stock] <- r
    spfits$k[spfits$stockid==stock] <- k
    # curve(r/p*x*(1-(x/k)^p), from=0, to=1, lty=1, lwd=0.9, col="black")
  }
  
  # Fit SP-SST model
  r <- spfits$r[spfits$stockid==stock]
  k <- spfits$r[spfits$stockid==stock]
  spfit <- fit_pella_sst_nls(sp=sdata$sp_sd, tb=sdata$tb_sd, sst=as.matrix(sdata[,sst_col]), p=p, r_start=r, k_start=k)
  if(!(inherits(spfit, "try-error"))){
    r <- exp(coef(spfit)["r"])
    k <- exp(coef(spfit)["k"])
    betaT <- coef(spfit)["betaT"]
    spfits$r2[spfits$stockid==stock] <- r
    spfits$k2[spfits$stockid==stock] <- k
    spfits$betaT[spfits$stockid==stock] <- betaT
    # curve(r/p*x*(1-(x/k)^p), from=0, to=1, lty=1, lwd=0.9, col="black")
  }
  
}

# Check completeness
sum(is.na(spfits$r))
sum(is.na(spfits$betaT))

# Export data
write.csv(spfits, paste(outputdir, "ramldb_v3.8_pella_0.2_cobe_nls.csv", sep="/"), row.names=F)

