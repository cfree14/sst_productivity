
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"

# Read data
data_orig <- read.csv(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.csv", sep="/"), as.is=T)


# Format data
################################################################################

# Format data
data <- data_orig %>% 
  select(stockid, betaT, betaT_inf, migratory, 
         repro_mode, repro_guild1, repro_guild2, spawning, spawning_ground) %>% 
  mutate(repro_guild1=revalue(repro_guild1, c("nonguarders"="non-\nguarders")),
         repro_guild1=factor(repro_guild1, levels=c("non-\nguarders", "guarders", "bearers")),
         repro_guild2=revalue(repro_guild2, c("External brooders"="external\nbrooders",
                                              "internal live bearers"="internal\nlive bearers",
                                              "open water/substratum egg scatterers"="egg\nscatterers")),
         repro_guild2=factor(repro_guild2, levels=c("egg\nscatterers", "nesters", "external\nbrooders", "internal\nlive bearers")),
         spawning=revalue(spawning, c("no obvious seasonal peak"="no peaks",
                                      "once in a lifetime"="1 / lifetime",
                                      "one clear seasonal peak per year"="1 peak / year",
                                      "two seasonal peaks per year"="2 peaks / year",
                                      "Two seasonal peaks per year"="2 peaks / year",
                                      "Variable throughout range"="variable\nthroughout range")),
         spawning=factor(spawning, levels=c("no peaks", "2 peaks / year", "1 peak / year", "1 / lifetime", "variable\nthroughout range")),
         spawning_ground=factor(spawning_ground, levels=c("coastal", "shelf", "oceanic")),
         migratory=factor(migratory, levels=c("catadromous", "anadromous", "oceano-estuarine", "oceanadromous", "non-migratory")))

# Values
table(data$repro_guild1)
table(data$repro_guild2)
table(data$spawning)

# Plot data
################################################################################

# Variables
vars <- c("repro_mode", "repro_guild1", "repro_guild2", "migratory", "spawning", "spawning_ground")
titles <- c("Reproductive mode", "Reproductive guild 1°", "Reproductive guild 2°", 
            "Migratory behaviour", "Spawning frequency", "Spawning grounds")
labels <- c("", "Parental investment →", "Parental investment →", 
            "Freshwater to saltwater →", "Parental investment →", "Inshore to offshore →")

# Setup figure
figname <- "SFig13_thetas_by_behavior.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=5, units="in", res=600)
par(mfrow=c(2,3), mar=c(5,2,1,1), mgp=c(2.5,0.8,0), xpd=NA, oma=c(2,2,0,0))

# Loop through vars
for(i in 1:length(vars)){

  # Plot
  boxplot(data$betaT ~ data[,vars[i]], frame=F, las=2, col="grey80", border="grey30",
          ylim=c(-1, 1), lty=1, lwd=0.8, cex.axis=0.9,
          xlab="", ylab="", main=titles[i])
  text(x=0.5, y=-0.9, adj=0, labels=labels[i], cex=0.9)
  
  # Sample size and line
  n <- table(data[,vars[i]])
  text(x=1:length(n), y=rep(0.9, length(n)), labels=n, col="grey20", cex=0.9)  
  lines(x=c(0.5,0.5+length(n)), y=c(0,0), lty=3, lwd=1.2)
  
}

# Add axis label
mtext(expression("SST influence (θ"["i"]*")"), outer=T, side=2, adj=0.57, line=0.5, cex=0.8)

# Off
dev.off()









