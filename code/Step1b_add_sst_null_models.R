

# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(forecast)
library(freeR)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/input"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Read data
load(paste(datadir, "ramldb_v3.8_production_data_clean.Rdata", sep="/"))
data <- final; rm(final)


# Build data
################################################################################

# Stocks
stocks <- sort(unique(data$stockid))

# ARIMA fit lists
arfit_list <- list()
arfit_n1_list <- list()
arima_matched <- rep(T, length(stocks))

# Add Null fields
data <- data %>%
  ungroup() %>% 
  mutate(cobe_sst_c_n1=NA,
         cobe_sst_c_n2=NA,
         cobe_sst_c_n3=NA)

# Loop through stocks
for(i in 1:length(stocks)){
  
  # Subset data
  stock <- stocks[i]
  sdata <- subset(data, stockid==stock)
  yrs <- sdata$year
  ssts <- sdata$cobe_sst_c
  print(paste(i, stock))
  
  # Analyze observed time series
  lmfit <- lm(cobe_sst_c ~ year, sdata)
  arfit <- auto.arima(ssts)
  
  # Extract statistics
  mu_obs <- mean(ssts)
  sd_obs <- sd(ssts)
  slope_obs <- coef(lmfit)[2]
  ar_pdq_obs <- arimaorder(arfit)
  sst_trend <- as.numeric(predict(lmfit, newdata=data.frame(year=sdata$year)))
  
  # Null 1: simulate ARIMA properties w/ trend
  # Simulate SST data using ARIMA fit until properties match observed properties
  counter <- 0
  counter_cutoff <- 500
  ar_pdq_n1 <- c(-999, -999, -999)
  while(sum(ar_pdq_n1==ar_pdq_obs)!=length(ar_pdq_obs)){
    counter <- counter + 1
    print(paste(i, stock, "- try", counter))
    # Simulate SST
    sst_ar <- as.numeric(arima.sim(model=as.list(coef(arfit)), n=length(ssts), sd=sqrt(arfit$sigma2)))
    sst_n1 <- sst_trend + sst_ar
    # Check simulated SST arima against observed SST arima
    arfit_n1 <- auto.arima(sst_n1)
    ar_pdq_n1 <- arimaorder(arfit_n1)
    # Check counter
    if(counter>=counter_cutoff){
      arima_matched[i] <- F
      break
    }
  }
  
  # Null 2: use ARIMA sim residuals
  sst_n2 <- mu_obs + sst_ar
  
  # Null 3: shuffle data
  sst_n3 <- sample(ssts)
  
  # Record results in data
  data$cobe_sst_c_n1[data$stockid==stock] <- sst_n1
  data$cobe_sst_c_n2[data$stockid==stock] <- sst_n2
  data$cobe_sst_c_n3[data$stockid==stock] <- sst_n3
  
  # Record model fits
  arfit_list[[i]] <- arfit
  arfit_n1_list[[i]] <- arfit_n1
  
}

# Inspect data
str(data)
apply(data, 2, function(x) sum(is.na(x)))


# Center simulated SST time series
################################################################################

# Center SST
final <- data %>%
  group_by(assessid) %>% 
  mutate(cobe_sst_c_n1_sd=cobe_sst_c_n1-mean(cobe_sst_c_n1),
         cobe_sst_c_n2_sd=cobe_sst_c_n2-mean(cobe_sst_c_n2),
         cobe_sst_c_n3_sd=cobe_sst_c_n3-mean(cobe_sst_c_n3))

# Confirm centered
data.frame(tapply(final$cobe_sst_c_sd, final$assessid,  mean, na.rm=T)) # should be 0
data.frame(tapply(final$cobe_sst_c_n1_sd, final$assessid,  mean, na.rm=T)) # should be 0
data.frame(tapply(final$cobe_sst_c_n2_sd, final$assessid,  mean, na.rm=T)) # should be 0
data.frame(tapply(final$cobe_sst_c_n3_sd, final$assessid,  mean, na.rm=T)) # should be 0

# Rename data
data <- final
complete(data)

# How many Null 1 models have different ARIMA order than observed?
sum(arima_matched==F)

# Export data
################################################################################

# Export data
save(data, spfits, 
     file=paste(datadir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))

# Export SST simulation fits
save(stocks, data, arfit_list, arfit_n1_list, arima_matched,
     file=paste(datadir, "ramldb_v3.8_sst_null_model_fits.Rdata", sep="/"))


