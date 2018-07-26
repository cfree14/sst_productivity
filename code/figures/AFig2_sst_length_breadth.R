
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Load data
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
rm(hess, results.wide, results.df, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)

# Read data
stocks <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)


# Stats for appendix
################################################################################

# Number of stocks in each species
nspp <- stocks %>% 
  group_by(species) %>% 
  summarize(n=n()) %>% 
  arrange(desc(n))
sum(nspp$n==1)

# Plot data
################################################################################

# Build data
df <- data %>%
  group_by(stockid) %>% 
  summarize(nyr=n(),
            sst_breadth=max(cobe_sst_c)-min(cobe_sst_c),
            startyr=min(year),
            endyr=max(year))

# Setup figure
figname <- "AFig2_sst_length_breadth.png"
png(paste(plotdir, figname, sep="/"), width=6, height=6, units="in", res=600)
par(mfrow=c(2,2), mar=c(3.5,1.5,0.5,0.5), mgp=c(2,0.7,0), oma=c(0,2,0,0))

# Plot SST breadths
range(df$sst_breadth)
sst_avg <- mean(df$sst_breadth)
hist(df$sst_breadth, breaks=seq(0,4.25,0.25),
     main="", xlab="SST breadth (°C)", ylab="", col="grey70", border=F,
     xlim=c(0,4.5), ylim=c(0,60), las=1, cex.axis=0.75, cex.lab=0.9, xaxt="n")
axis(1, at=seq(0, 4.5, 0.5), cex.axis=0.75)
lines(x=c(sst_avg,sst_avg), y=c(0,60), lty=2)
text(x=sst_avg, y=60, pos=4, label=paste0(roundf(sst_avg,1), "°C"), cex=0.8)
# mtext("A", side=3, adj=0.05, line=-1.8, font=2, cex=1)

# Plot time series lengths
range(df$nyr)
nyr_avg <- mean(df$nyr)
hist(df$nyr, main="", breaks=seq(0,100,5),
     xlab="Time series length (yr)", ylab="", col="grey70", border=F,
     xlim=c(0, 100), ylim=c(0,40), las=1, cex.axis=0.7, xpd=NA, cex.lab=0.9)
lines(x=c(nyr_avg,nyr_avg), y=c(0,40), lty=2)
text(x=nyr_avg, y=40, pos=4, label=paste(roundf(nyr_avg,1), " years"), cex=0.8)
# mtext("B", side=3, adj=0.05, line=-1.8, font=2, cex=1)

# Plot start years
range(df$startyr)
yr <- median(df$startyr)
hist(df$startyr, breaks=seq(1900,2000, 5),
     main="", xlab="Start year", ylab="", col="grey70", border=F,
     ylim=c(0,50), las=1, cex.axis=0.75, cex.lab=0.9)
lines(x=c(yr,yr), y=c(0,50), lty=2)
text(x=yr, y=50, pos=4, label=yr, cex=0.8)
# mtext("C", side=3, adj=0.05, line=-1.8, font=2, cex=1)

# Plot end years
range(df$endyr)
yr <- median(df$endyr)
hist(df$endyr, breaks=seq(1990,2020,1),
     main="", xlab="End year", ylab="", col="grey70", border=F,
     ylim=c(0,50), las=1, cex.axis=0.75, cex.lab=0.9)
lines(x=c(yr,yr), y=c(0,50), lty=2)
text(x=yr, y=50, pos=4, label=yr, cex=0.8)
# mtext("D", side=3, adj=0.05, line=-1.8, font=2, cex=1)

# Add y-axis label
mtext("Number of stocks", outer=T, side=2, adj=0.55, line=0.5, cex=1)

# Off
dev.off()
graphics.off()
