
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"

# A. Monotonic SST dependence
############################################

# Parameters
ssts <- seq(-1,1,0.1)
thetas <- seq(-0.8, 0.8, 0.2)
colors <- colorpal(brewer.pal(11, "RdBu"), length(thetas))

# Calculate scalars
scalars <- sapply(thetas, function(x) exp(ssts*x))

# Plot monotonic influence
range(scalars)
plot(1:10, 1:10, type="n", bty="n",
     xlim=c(-1,1), ylim=c(0, 2.5), las=1,
     xlab="SST anamoly (°C)", ylab="Production scalar")
for(i in 1:ncol(scalars)){
  lines(x=ssts, y=scalars[,i], col=colors[i])
}


# B. Dome-shaped SST dependence
############################################

# Parameters
ssts <- seq(0,30,0.1)
thetas <- seq(0, 1, 0.1)
colors <- colorpal(brewer.pal(9, "Oranges")[3:9], length(thetas))

# Calculate scalars
z <- 15
scalars <- sapply(thetas, function(x) exp(-(ssts-z)^2*x))

# Plot monotonic influence
range(scalars)
plot(1:10, 1:10, type="n", bty="n",
     xlim=c(0,30), ylim=c(0, 1), las=1,
     xlab="SST anamoly (°C)", ylab="Production scalar")
for(i in 1:ncol(scalars)){
  lines(x=ssts, y=scalars[,i], col=colors[i])
}




