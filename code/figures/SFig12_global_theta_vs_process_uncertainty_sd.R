
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
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Read random effects results
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)


# Plot data
################################################################################

# Global SST influence values
theta_mu <- results.df$estimate[results.df$param=="mu_T"]
theta_sd <- results.df$estimate[results.df$param=="sd_T"]

# Setup figure
figname <- "SFig12_global_theta_uncertainty_vs_process_uncertainty.png"
png(paste(plotdir, figname, sep="/"), width=4, height=3, units="in", res=600)
par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), mgp=c(1.9,0.6,0))

# Plot data
hist(results.wide$sigmaP, las=1, cex.axis=0.8, cex.lab=0.9, breaks=seq(0,0.35,0.01), border=F,
     main="", xlab=expression("σ"["P,i"]), col="grey80", ylim=c(0,30))

# Add standard deviation of global SST influence
lines(x=rep(theta_sd,2), y=c(0,30), lty=3)
text(x=theta_sd, y=30, pos=2, labels=expression("σ"["SST"]), cex=0.7, offset=0.1)

# Off
dev.off()
graphics.off()



