
# Clear workspace
rm(list = ls())

# Turn off scientific notation
options(scipen=999)

# Read data
################################################################################

# Data requirements:
# r/K parameters - estimated
# B0 - observed
# SST time series - observed
# Catch time series - observed
# Exploitation time series - observed

# Packages
library(plyr)
library(dplyr)
library(ggplot2)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"

# Read model output and input
load(paste(datadir, "ramldb_v3.8_spmodel_temp_dpdt_sst_yr_t.Rdata", sep="/"))
rm(sd, results.df, output, model, hess, input.data, problem.stocks, params, nstocks, stocks)


# Build data
################################################################################

# Identify B0 for stocks
b0 <- data %>% 
  group_by(stockid) %>% 
  summarize(tb_max=max(tb),
            b0=tb[year==min(year)],
            nyr=n())

# Build params
stocks <- results.wide %>% 
  # Important columns
  select(stockid, r, B0) %>% 
  rename(k_orig=B0) %>% 
  # Add B0, TBmax, and nyr
  left_join(b0, by="stockid") %>% 
  # Calculate K in 1000s MT
  mutate(k=k_orig*tb_max)

# Add exploitation rate
data <- data %>% 
  mutate(u_orig=catch/tb, 
         u=pmin(u_orig, 0.95))

# Inspect ER > 1
range(data$u_orig)
hist(data$u_orig, breaks=seq(0,3,0.1), main="", xlab="Exploitation rate", col="grey60", border=F)
er_large <- subset(data, u_orig>0.95)


# Simulation parameters
################################################################################

# Simulation parameters
betaT_mu_vals <- c(-0.1, 0, 0.1)
betaT_sd_vals <- c(0.2, 0.3, 0.4)
sigmaP_vals <- c(0.05, 0.1, 0.2)

# Build scenarios
scenarios <- expand.grid(betaT_mu=betaT_mu_vals, 
                         betaT_sd=betaT_sd_vals, 
                         sigmaP=sigmaP_vals, stringsAsFactors=F)
scenarios <- scenarios %>% 
  mutate(scenario=paste(betaT_mu, betaT_sd, sigmaP, sep="_")) %>% 
  select(scenario, everything())


# Simulate data
################################################################################

# Simulate time series
# For testing: l <- i <- 1; j <- 2
for(l in 1:nrow(scenarios)){
  
  # Get scenario data
  scenario <- scenarios$scenario[l]
  betaT_mu <- scenarios$betaT_mu[l]
  betaT_sd <- scenarios$betaT_sd[l]
  sigmaP <- scenarios$sigmaP[l]
  print(paste0("Scenario ", l, " of ", nrow(scenarios), " - ", scenario))
  
  # Loop through stocks
  for(i in 1:nrow(stocks)){
  
    # Get stock data
    stock <- stocks$stockid[i]
    r <- stocks$r[i]
    k <- stocks$k[i]
    b0_obs <- stocks$b0[i]
    nyr <- stocks$nyr[i]
    betaT <- rnorm(1, betaT_mu, betaT_sd)
    
    # Get stock time series
    yrs <- data$year[data$stockid==stock]
    C_t_obs <- data$catch[data$stockid==stock]
    B_t_obs <- data$tb[data$stockid==stock]
    U_t_obs <- data$u[data$stockid==stock]
    P_t_obs <- data$sp[data$stockid==stock]
    Temp_t_obs <- data$sst_t_sd[data$stockid==stock]
  
    # Setup empty time series vectors
    P_t <- B_t <- C_t <- rep(NA, nyr)
  
    # Initial conditions
    B_t[1] <- b0_obs
  
    # Simulate biomass
    for(j in 2:nyr){
      P_t[j-1] <- r*B_t[j-1]*(1-B_t[j-1]/k) + rnorm(1, mean=0, sd=sigmaP)
      C_t[j-1] <- B_t[j-1] * U_t_obs[j-1]
      B_t[j] <- B_t[j-1] + P_t[j-1] - C_t[j-1] 
    }
    
    # Merge results
    mat <- data.frame(scenario=scenario, betaT_mu=betaT_mu, betaT_sd=betaT_sd, sigmaP=sigmaP,
                      stockid=stock, year=yrs, 
                      betaT=betaT, sst_obs=Temp_t_obs, u_obs=U_t_obs,
                      ct_obs=C_t_obs, bt_obs=B_t_obs, pt_obs=P_t_obs,
                      ct_sim=C_t, bt_sim=B_t, pt_sim=P_t, stringsAsFactors=F)
    if(i==1){sim_s <- mat}else{sim_s <- rbind(sim_s, mat)}
    
  }
  
  # Merge results
  if(l==1){sim_all <- sim_s}else{sim_all <- rbind(sim_all, sim_s)}
  
}

# Add ID to simulations
sim_all <- sim_all %>% 
  mutate(id=paste(scenario, stockid, sep="_")) %>% 
  select(id, everything())


# Inspect simulations
################################################################################

# Build simulation key
sim_key_all <- sim_all %>%
  group_by(id, scenario, betaT_mu, betaT_sd, sigmaP, stockid) %>% 
  summarize(betaT=unique(betaT),
            bneg = sum(bt_sim<0)>0 | sum(is.na(bt_sim))>0)

# Build final simulation key
sim_key <- sim_key_all %>% 
  filter(bneg==F) %>% 
  select(-bneg) %>% 
  left_join(stocks, by=c("stockid"))

# Inspect sample size across scenarios
scenarios_n <- sim_key %>% 
  group_by(scenario, betaT_mu, betaT_sd, sigmaP) %>% 
  summarize(n=n())
hist(scenarios_n$n, xlab="# of stocks (201 max)", main="")

# Format simulation data
sim_data <- sim_all %>% 
  filter(id %in% sim_key$id)


# Export simulations
################################################################################

# Export data
save(sim_key, sim_data, file=paste(datadir, "simulations_self.Rdata", sep="/"))


# Plot sample simulations
################################################################################

# Sample scenario
ex_scenario <- scenarios$scenario[1]
ex_stocks <- subset(sim_key, scenario==ex_scenario)

# Setup figure
figname <- "AppendixE_simulations_self_example.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(8,6), oma=c(4,4,4,4), mar=c(0.5,1.3,0.5,1), mgp=c(3, 0.5, 0))

# Top of page index
top_i <- seq(1,nrow(ex_stocks),8)

# Loop through stocks
for(i in 1:nrow(ex_stocks)){
  
  # Subset data
  stock <- ex_stocks$stockid[i]
  sdata <- subset(sim_data, stockid==stock & scenario==ex_scenario)
  yr1 <- min(sdata$year)
  yr2 <- max(sdata$year)
  print(paste(i, stock))
  
  # Stock info
  r <- round(ex_stocks$r[i], 2)
  k <- round(ex_stocks$k[i], 0)
  b0 <- ex_stocks$b0[i]
  betaT <- round(ex_stocks$betaT[i], 3)
  betaT_mu <- ex_stocks$betaT_mu[i]
  betaT_sd <- ex_stocks$betaT_sd[i]
  sigmaP <- ex_stocks$sigmaP[i]
  
  # Derive ref points
  fmsy <- r/2
  bmsy <- k/2
  msy <- r*k/4

  # Plot scenario details
  plot(0:10, 0:10, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  if(i%in%top_i){title("Scenario", line=1, xpd=NA)}
  sim_text <- paste0(stock, "\n",
                     "r = ", r, "\n",
                     "K = ", k, "\n",
                     "B0 = ", b0, "\n",
                     "0i = ", betaT, "\n",
                     "0i ~ N(", betaT_mu, ", ", betaT_sd, ")\n",
                     "SigmaP = ", sigmaP)
  text(x=9, y=5, adj=1, labels=sim_text, cex=1, xpd=NA)
  
  # Plot SST
  tmin <- floor(min(sdata$sst_obs) / 0.1) * 0.1
  tmax <- ceiling(max(sdata$sst_obs) / 0.1) * 0.1
  plot(sst_obs ~ year, sdata, type="l", bty="n", xaxt="n", yaxt="n", ylim=c(tmin, tmax),
       xlab="", ylab="", col="red", lwd=0.6)
  axis(2, at=c(tmin,tmax), labels=T, las=1, lwd=0.5, lwd.ticks=-1)
  lines(x=c(yr1,yr2), y=c(0, 0), lty=3)
  if(i%in%top_i){title("SST anomaly (obs)", line=1, xpd=NA)}
  text(x=yr1-1, y=0, pos=2, labels="0°C", xpd=NA)

  # Plot exploitation rate
  umax <- ceiling(max(sdata$u_obs, na.rm=T) / 0.1) * 0.1
  plot(u_obs ~ year, sdata, type="l", bty="n", xaxt="n", yaxt="n", ylim=c(0, umax),
       xlab="", ylab="", col="orange", lwd=0.6)
  axis(2, at=c(0,umax), labels=T, las=1, lwd=0.5, lwd.ticks=-1)
  if(i%in%top_i){title("ER (obs)", line=1, xpd=NA)}

  # Plot catch
  cmax <- ceiling(max(sdata$ct_sim, na.rm=T) / 0.1) * 0.1
  plot(ct_sim ~ year, sdata, type="l", bty="n", xaxt="n", yaxt="n", ylim=c(0, cmax),
       xlab="", ylab="", col="blue", lwd=0.6)
  axis(2, at=c(0,cmax), labels=T, las=1, lwd=0.5, lwd.ticks=-1)
  if(i%in%top_i){title("Catch", line=1, xpd=NA)}
  lines(x=c(yr1,yr2), y=c(msy, msy), lty=3)
  
  # Plot biomass
  bmax <- ceiling(max(sdata$bt_sim) / 0.1) * 0.1
  plot(bt_sim ~ year, sdata, type="l", bty="n", xaxt="n", yaxt="n", ylim=c(0,bmax),
       xlab="", ylab="", col="darkgreen", lwd=0.6)
  axis(2, at=c(0,bmax), labels=T, las=1, lwd=0.5, lwd.ticks=-1)
  if(i%in%top_i){title("Biomass", line=1, xpd=NA)}
  lines(x=c(yr1,yr2), y=c(bmsy, bmsy), lty=3)
  
  # Plot production
  pmin <- floor(min(sdata$pt_sim, na.rm=T) / 0.1) * 0.1
  pmax <- ceiling(max(sdata$pt_sim, na.rm=T) / 0.1) * 0.1
  plot(pt_sim ~ year, sdata, type="l", bty="n", xaxt="n", yaxt="n", ylim=c(pmin, pmax),
       xlab="", ylab="", col="purple", lwd=0.6)
  axis(2, at=c(pmin,pmax), labels=T, las=1, lwd=0.5, lwd.ticks=-1)
  if(i%in%top_i){title("Production", line=1, xpd=NA)}

}

# Off
dev.off()
graphics.off()

