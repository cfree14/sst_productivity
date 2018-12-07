

# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)
library(reshape2)

# Directories
outdir1 <- "output"
outdir2 <- "revisions/output"
plotdir <- "revisions/figures"

# Load results
results1 <- read.csv(file.path(outdir1, "ramldb_v3.8_spsst_pella_cobe_lme.csv"), as.is=T)
results2 <- read.csv(file.path(outdir2, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.csv"), as.is=T)

# Build and plot data
################################################################################

# Build data
data <- results1 %>% 
  select(stockid, betaT, betaT_lo, betaT_hi, betaT_inf) %>% 
  left_join(select(results2, stockid, betaT, betaT_lo, betaT_hi, betaT_inf), by="stockid")

# Export figure
figname <- "FigX_comparison_of_initial_revised_influences.png"
png(paste(plotdir, figname, sep="/"), width=4, height=4, units="in", res=600)
par(mar=c(4,4,0.5,0.5), mgp=c(2.5,0.8,0))

# Plot data
plot(betaT.y ~ betaT.x, data, bty="n", las=1, xlim=c(-1,1), ylim=c(-1,1), 
     xlab="SST influence (indpt error, original)", ylab="SST influence (AR1 error, revised)", 
     col="grey70")
abline(h=0, lty=3)
abline(v=0, lty=3)
abline(a=0, b=1)

# Off
dev.off()



