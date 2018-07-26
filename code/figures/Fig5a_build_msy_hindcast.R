
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/cobe"

# Read SST data
sst <- read.csv(paste(sstdir, "sst_mean_annual_cobe.csv", sep="/"), as.is=T)

# Read model data
stocks_orig <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)

# Read hindcast data
data_orig <- read.csv(paste(datadir, "msy_hindcast_time_series_using_lme_model.csv", sep="/"), as.is=T)


# Format data
################################################################################

# Hindcast window
# yr1 <- 1930
# yr2 <- 2010

# Format data
data <- data_orig 
# %>% filter(year>=yr1 & year<=yr2)

# Format stocks
p <- 0.2
div <- (p+1)^((p+1)/p)
stocks <- stocks_orig %>%
  rename(k_orig = k) %>% 
  mutate(k = k_orig * tb_max,
         msy_avg = (r * k) / div) 

# Format SST data
avg_yrs <- 1971:2000
sst_avg <- mean(sst$sst_c[sst$year %in% avg_yrs])
sst <- mutate(sst, sst_c_sd=sst_c-sst_avg)

# MSY total
msy_avg <- sum(stocks$msy_avg)/1E6
msy_avg_pos <- sum(stocks$msy_avg[stocks$betaT_inf=="positive"])/1E6
msy_avg_neg <- sum(stocks$msy_avg[stocks$betaT_inf=="negative"])/1E6
msy_avg_non <- sum(stocks$msy_avg[stocks$betaT_inf=="none"])/1E6


# Bootstrap MSY trajectories
################################################################################

# Parameters
years <- sst$year
nyrs <- length(years)
nstocks <- nrow(stocks)
ntrajs <- 1000

# Setup container
msy_trajs <- array(data=NA, dim=c(nyrs, ntrajs, nstocks), dimnames=list(years, 1:ntrajs, stocks$stockid))

# Loop through stocks
for(i in 1:nrow(stocks)){
  
  # Stock
  stock <- stocks$stockid[i]
  sdata <- subset(data, stockid==stock)
  
  # Draw 1000 thetas
  betaT <- stocks$betaT[i]
  betaT_se <- stocks$betaT_se[i]
  thetas <- rnorm(ntrajs, mean=betaT, sd=betaT_se)
  # hist(thetas)
  
  # MSY over time for each theta
  r <- stocks$r[i]
  k <- stocks$k[i]
  sst_sd <- sdata$sst_c_sd
  msys <- sapply(thetas, function(x) r * k / div * exp(sst_sd*x))
  # tgray <- rgb(t(col2rgb("grey40"))/255, alpha=0.4)
  # plot(x=1940:2010, y=msys[,1], type="l", ylim=c(min(msys), max(msys)), col=tgray, las=1)
  # for(j in 2:ncol(msys)){lines(x=1940:2010, y=msys[,j], type="l", col=tgray)}
  
  # Merge MSY trajectories of 1 stock in array of all stocks
  msy_trajs[,,i] <- msys
  
  
}


# Helper functions
################################################################################

# Function to format quantile data
format_data <- function(q_data){
  q_data_df <- data.frame(t(q_data)) %>% 
    rename(msy=X50.,
           msy_lo=X2.5.,
           msy_hi=X97.5.) %>% 
    mutate(year=row.names(t(q_data)),
           msy=msy/1E6,
           msy_lo=msy_lo/1E6,
           msy_hi=msy_hi/1E6) %>% 
    select(year, msy, msy_lo, msy_hi)
  return(q_data_df)
}


# Summarize by Large Marine Ecoregion (LME)
################################################################################

# LMEs and container
lmes <- sort(unique(stocks$lme_name))
msy_trajs_sum_lme <- array(data=NA, dim=c(nyrs, ntrajs, length(lmes)), dimnames=list(years, 1:ntrajs, lmes))

# Loop through LMEs: i <- 1
for(i in 1:length(lmes)){
  
  # Subset MSY trajectories
  lme <- lmes[i]
  print(paste(i, lme))
  msy_trajs_lme <- msy_trajs[,,which(stocks$lme_name %in% lme)]
  
  # Sum MSY trajectories
  msy_trajs_lme_sum <- apply(msy_trajs_lme, 1:2, sum)
  msy_trajs_sum_lme[,,i] <- msy_trajs_lme_sum
  
  # Calculate MSY quantiles
  msy_trajs_lme_q <- apply(msy_trajs_lme_sum, 1, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
  
  # Format MSY quantiles
  msy_trajs_lme_q_df <- format_data(msy_trajs_lme_q)
  
  # Add LME
  msy_trajs_lme_q_df1 <- msy_trajs_lme_q_df %>% 
    mutate(lme=lme) %>% 
    select(lme, year, msy, msy_lo, msy_hi)
  
  # Merge data
  if(i==1){msy_lme_df <- msy_trajs_lme_q_df1}else{msy_lme_df <- rbind(msy_lme_df, msy_trajs_lme_q_df1)}
  
}


# Summarize for total, pos, neg, and non-significant stocks
################################################################################

# Subset significant positve, negative, and non-sig
msy_trajs_pos <- msy_trajs[,,which(stocks$betaT_inf=="positive")]
msy_trajs_neg <- msy_trajs[,,which(stocks$betaT_inf=="negative")]
msy_trajs_non <- msy_trajs[,,which(stocks$betaT_inf=="none")]

# Sum total MSY trajectories
# Creates a matrix with rows=years and columns=sum of all stocks for each runs
msy_trajs_sum <- apply(msy_trajs, 1:2, sum)
msy_trajs_pos_sum <- apply(msy_trajs_pos, 1:2, sum)
msy_trajs_neg_sum <- apply(msy_trajs_neg, 1:2, sum)
msy_trajs_non_sum <- apply(msy_trajs_non, 1:2, sum)

# Calculate quantiles
msy_trajs_q <- apply(msy_trajs_sum, 1, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
msy_trajs_pos_q <- apply(msy_trajs_pos_sum, 1, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
msy_trajs_neg_q <- apply(msy_trajs_neg_sum, 1, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
msy_trajs_non_q <- apply(msy_trajs_non_sum, 1, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))

# Format quantiles
msy_trajs_df <- format_data(msy_trajs_q)
msy_trajs_pos_df <- format_data(msy_trajs_pos_q)
msy_trajs_neg_df <- format_data(msy_trajs_neg_q)
msy_trajs_non_df <- format_data(msy_trajs_non_q)


# Export data
################################################################################

# Export data
save(stocks, sst,
     msy_trajs_sum, msy_trajs_sum_lme,
     msy_avg, msy_avg_pos, msy_avg_neg, msy_avg_non,
     msy_trajs_df, msy_trajs_pos_df, msy_trajs_neg_df, msy_trajs_non_df, msy_lme_df,
     file=paste(datadir, "spsst_pella40%_cobe_lme_msy_hindcast_1000traj.Rdata", sep="/"))


