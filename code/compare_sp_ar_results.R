

# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"

# Read SP data
load(file.path(outputdir, "ramldb_v3.8_sp.Rdata"))
sp <- results.wide

# Read SP AR1 data
load(file.path(outputdir, "ramldb_v3.8_sp_ar1_rho0.Rdata"))
ar <- results.wide

# Merge data
data <- sp %>% 
  left_join(ar, by="stockid")

# Plot r
par(mfrow=c(1,1))
plot(r.y ~ r.x, data)
plot(B0.y ~ B0.x, data)
plot(sigmaP.y ~ sigmaP.x, data)



