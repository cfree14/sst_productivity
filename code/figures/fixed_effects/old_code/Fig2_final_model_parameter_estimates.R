
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
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures/fixed_effects"

# Read fixed effects results
load(paste(datadir, "ramldb_v3.8_pella_0.20_cobe_fixed.Rdata", sep="/"))
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd, results.df)


# Plot data
################################################################################

# Setup figure
figname <- "Fig2_final_model_parameter_estimates.png"
png(paste(plotdir, figname, sep="/"), width=6, height=6, units="in", res=600)
par(mfrow=c(2,2), mar=c(3,0.5,0.5,0.5), mgp=c(2,0.8,0))

# Loop through params
params <- c("r", "k", "betaT", "sigmaP")
labels <- c( expression("Intrinsic rate of growth (r"["i"]*")"),
            expression("Carrying capacity (K"["i"]*")"),
            expression("SST influence (θ"["i"]*")"),
            expression("Process uncertainty (σ"["P,i"]*")"))
xlim <- matrix(c(0,2,0.5,
                 0,8,2,
                 -8,14,2,
                 0,0.4,0.1), byrow=T, ncol=3)


# Loop through parameters
for(i in 1:length(params)){

  # Subset data
  vals <- data.frame(est=results.wide[,params[i]],
                     est_lo=results.wide[,paste0(params[i], "_lo")],
                     est_hi=results.wide[,paste0(params[i], "_hi")])
  vals <- arrange(vals, desc(est))

  # Plot data
  xmin <- xlim[i,1]
  xmax <- xlim[i,2]
  xbin <- xlim[i,3]
  plot(1:10, 1:10, type="n", bty="n", xaxt="n", yaxt="n", cex.axis=0.75,
       xlim=c(xmin, xmax), ylim=c(1,nrow(vals)),
       xlab=labels[i], ylab="", cex.main=0.8)
  axis(1, at=seq(xmin, xmax, xbin), cex.axis=0.75)

  # Add estimates
  sapply(1:nrow(vals), function(x) lines(x=c(vals$est_lo[x], vals$est_hi[x]), y=c(x,x), col="grey60", lwd=0.6))
  points(vals$est, 1:nrow(vals), pch=16, cex=0.5, col="black")

  # Add vertical line
  if(params[i]=="k"){lines(x=c(1,1), y=c(1, nrow(vals)), lty=3)}
  if(params[i]=="betaT"){lines(x=c(0,0), y=c(1, nrow(vals)), lty=3)}

}

# Off
dev.off()
graphics.off()


