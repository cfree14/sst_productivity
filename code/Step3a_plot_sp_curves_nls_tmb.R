
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
inputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/input"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"

# Read NLS SP-model fits
sp_nls <- read.csv(paste(outputdir, "ramldb_v3.8_spmodel_nls.csv", sep="/"), as.is=T)

# Read TMB SP-model fits
load(paste(outputdir, "ramldb_v3.8_spmodel.Rdata", sep="/"))
sp_tmb <- results.wide
rm(data, hess, results.df, results.wide, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)

# Read TMB SP-SST-model fits
load(paste(outputdir, "ramldb_v3.8_spmodel_temp_dpdt_sst_yr_t.Rdata", sep="/"))
sp_sst_tmb <- results.wide
rm(hess, results.df, results.wide, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)

# Read data
load(paste(inputdir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))
data_orig <- data; data <- final; rm(final, spfits)


# Build data
################################################################################

# Build data
spfits <- sp_nls %>% 
  # Rename NLS columns
  rename(r_nls=r, k_nls=k,
         r_nls_sst=r2, k_nls_sst=k2, betaT_nls_sst=betaT) %>% 
  # Add TMB values 
  left_join(select(sp_tmb, stockid, r, B0), by="stockid") %>% 
  rename(r_tmb=r, k_tmb=B0) %>% 
  # Add TMB-SST values
  left_join(select(sp_sst_tmb, stockid, r, B0, BetaT), by="stockid") %>% 
  rename(r_tmb_sst=r, k_tmb_sst=B0, betaT_tmb_sst=BetaT)


# Plot data
################################################################################

# Plot surplus production curves
stocks <- sort(unique(spfits$stockid))

# Setup figure
figname <- "AppendixC_spmodel_curves_nls_tmb.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 4), mar=c(1, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,3,3,3), lwd=0.8)

# Loop through stocks: i <- 1
for(i in 1:length(stocks)){
  
  # Subset data
  stock <- as.character(stocks[i])
  sdata <- data[data$stockid==stock,]
  sp <- sdata$sp_sd
  tb <- sdata$tb_sd
  print(paste(i, stock))
  
  # Plot data
  ymax <- max(sp)
  plot(sp ~ tb, bty="n", las=1, pch=16, xlab="", ylab="", xlim=c(0,1))
  mtext(stock, side=3, adj=0.5, line=0.2, cex=0.7, font=2)
  
  # Plot NLS fit
  r <- spfits$r_nls[spfits$stockid==stock]
  k <- spfits$k_nls[spfits$stockid==stock]
  curve(r*x*(1-x/k), 
        from=0, to=1, add=T, lty=1, lwd=1.1, col="darkgreen")
  
  # Plot NLS-SST fit
  r <- spfits$r_nls_sst[spfits$stockid==stock]
  k <- spfits$k_nls_sst[spfits$stockid==stock]
  curve(r*x*(1-x/k), 
        from=0, to=1, add=T, lty=1, lwd=1.1, col="green")
  
  # Plot TMB fit
  r <- spfits$r_tmb[spfits$stockid==stock]
  k <- spfits$k_tmb[spfits$stockid==stock]
  curve(r*x*(1-x/k), 
        from=0, to=1, add=T, lty=1, lwd=1.1, col="red")
  
  # Plot TMB fit
  r <- spfits$r_tmb_sst[spfits$stockid==stock]
  k <- spfits$k_tmb_sst[spfits$stockid==stock]
  curve(r*x*(1-x/k), 
        from=0, to=1, add=T, lty=1, lwd=1.1, col="orange")
  
}

# Off
dev.off()
graphics.off()


