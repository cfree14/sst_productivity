
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/code"

# Read helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# Read results
load(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.Rdata", sep="/"))
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)

# Format results
results.wide <- format_results(results.wide)


# Plot data
################################################################################

# Setup figure
figname <- "SFig7_final_model_parameter_estimates.png"
png(paste(plotdir, figname, sep="/"), width=6, height=4, units="in", res=600)
par(mfrow=c(2,3), mar=c(3,0.5,0.5,0.5), mgp=c(2,0.7,0))

# Loop through params
params <- c("r", "k", "betaT", "sigmaP", "rho")
labels <- c(expression("Intrinsic rate of growth (r"["i"]*")"),  
            expression("Carrying capacity (K"["i"]*")"), 
            expression("SST influence (θ"["i"]*")"), 
            expression("Residual process variability (σ"["P,i"]*")"),
            expression("AR1 correlation (ρ"["i"]*")"))
xlim <- matrix(c(0,2,
                 0,6,
                 -1.5,1.5,
                 0,0.4,
                 -1,1.5), byrow=T, ncol=2)


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
  plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=0.85,
       xlim=c(xmin, xmax), ylim=c(1,nrow(vals)),
       xlab=labels[i], ylab="", cex.main=0.8)
  
  # Add estimates
  sapply(1:nrow(vals), function(x) lines(x=c(vals$est_lo[x], vals$est_hi[x]), y=c(x,x), col="grey60", lwd=0.6))
  points(vals$est, 1:nrow(vals), pch=16, cex=0.5, col="black")
  
  # Add vertical line
  if(params[i]=="k"){lines(x=c(1,1), y=c(1, nrow(vals)), lty=3)}
  if(params[i]%in%c("betaT", "rho")){lines(x=c(0,0), y=c(1, nrow(vals)), lty=3)}
  
}

# Off
dev.off()
graphics.off()


