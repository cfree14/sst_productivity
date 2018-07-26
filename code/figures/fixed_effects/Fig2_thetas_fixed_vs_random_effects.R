
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output/fixed_effects"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures/fixed_effects"

# Read fixed effects results
load(paste(datadir, "ramldb_v3.8_pella_0.20_cobe_fixed.Rdata", sep="/"))
data_f <- results.wide
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd, results.df, results.wide)

# Read random effects results
load(paste(datadir, "ramldb_v3.8_pella_0.20_cobe_random_normal.Rdata", sep="/"))
data_r <- results.wide
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd, results.df, results.wide)

# Build data
################################################################################

# Merge data
data <- data_f %>%
  select(stockid, betaT, betaT_lo, betaT_hi, betaT_inf) %>%
  rename(betaT_f=betaT, betaT_lo_f=betaT_lo, betaT_hi_f=betaT_hi, betaT_inf_f=betaT_inf) %>%
  left_join(select(data_r, stockid, betaT, betaT_lo, betaT_hi, betaT_inf), by="stockid") %>%
  rename(betaT_r=betaT, betaT_lo_r=betaT_lo, betaT_hi_r=betaT_hi, betaT_inf_r=betaT_inf)


# Plot data
################################################################################

# Setup figure
figname <- "Fig2_thetas_fixed_vs_random_effects.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=4, units="in", res=600)
layout(matrix(c(1,2,
                1,3), ncol=2, byrow=T), widths=c(0.6, 0.4))
par(mar=c(3.5,3.5,1,0.5), mgp=c(2.2,0.8,0))

# A. Scatterplot
#####################################

# Setup empty plot
plot(betaT_f ~ betaT_r, data, type="n", bty="n",
     xlim=c(-1.5, 1.5), ylim=c(-8,14), pch=16,
     xlab=expression("θ"["random"]), ylab=expression("θ"["fixed"]), yaxt="n")
axis(2, at=seq(-8,14,2), las=2)
lines(x=c(-1.5,1.5), y=c(0,0), lty=2)
lines(x=c(0,0), y=c(-8,10), lty=2)

# Add error bars
for(i in 1:nrow(data)){lines(x=c(data$betaT_lo_r[i], data$betaT_hi_r[i]),
                             y=c(data$betaT_f[i], data$betaT_f[i]), col="grey60", lwd=0.6)}
for(i in 1:nrow(data)){lines(x=c(data$betaT_r[i], data$betaT_r[i]),
                             y=c(data$betaT_lo_f[i], data$betaT_hi_f[i]), col="grey60", lwd=0.6)}

# Add points
points(data$betaT_r, data$betaT_f, pch=16)

# B. Random effects
#####################################

# Fixed effects
hist(data$betaT_r, breaks=seq(-4,8,0.5), col="grey60", border=F, las=1,
     xlim=c(-4,8), xaxt="n", xlab=expression("θ"["random"]), main="")
axis(1, at=seq(-4,8,2))


# A. Fixed effects
#####################################

# Fixed effects
hist(data$betaT_f, breaks=seq(-4,8,0.5), col="grey60", border=F, las=1,
     xlim=c(-4,8), xaxt="n", xlab=expression("θ"["fixed"]), main="")
axis(1, at=seq(-4,8,2))


# Off
dev.off()
graphics.off()


