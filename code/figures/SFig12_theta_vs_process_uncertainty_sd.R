
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

# Setup figure
figname <- "SFig12_theta_vs_process_uncertainty.png"
png(paste(plotdir, figname, sep="/"), width=5, height=5, units="in", res=600)
par(mfrow=c(2,2), mar=c(3.3,4,1,0.5), mgp=c(2.3,0.8,0))

# betaT ~ sigmaP
colors <- c("red", "grey60", "blue")[factor(results.wide$betaT_inf)]
plot(betaT ~ sigmaP, results.wide, bty="n", las=1, col=colors, cex.axis=0.9,
     xlab=expression("σ"["P,i"]*" estimate"), ylab=expression("θ"["i"]*" estimate"))
lines(x=c(0, 0.3), y=c(0,0), lty=3, col="grey10")

# betaT ~ sigmaP_se
plot(betaT ~ sigmaP_se, results.wide, bty="n", las=1, col=colors, cex.axis=0.9,
     xlab=expression("σ"["P,i"]*" standard error"), ylab=expression("θ"["i"]*" estimate"))
lines(x=c(0, 0.04), y=c(0,0), lty=3, col="grey10")

# betaT ~ sigmaP
plot(betaT_se ~ sigmaP, results.wide, bty="n", las=1, col=colors, cex.axis=0.9,
     xlab=expression("σ"["P,i"]*" estimate"), ylab=expression("θ"["i"]*" standard error"))

# betaT_se ~ sigmaP
plot(betaT_se ~ sigmaP_se, results.wide, bty="n", las=1, col=colors, cex.axis=0.9,
     xlab=expression("σ"["P,i"]*" standard error"), ylab=expression("θ"["i"]*" standard error"))

# Off
dev.off()
graphics.off()


