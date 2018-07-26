
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_spsst_cobe.csv", sep="/"), as.is=T)


# Plot data
################################################################################

# Add colors
data <- data %>% 
  mutate(color=revalue(betaT_inf, c("positive"="blue",
                                    "negative"="red",
                                    "none"="grey60")))

  

# Setup figure
figname <- "Fig2_theta_overdispersion.png"
png(paste(plotdir, figname, sep="/"), width=4, height=4, units="in", res=600)
par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), mgp=c(2.5,0.8,0))

# Build data
qq <- data.frame(qqnorm(data$betaT, plot.it=F), col=data$color, stringsAsFactors=F)

# Plot data
xlabel <- expression(paste("Theoretical ", "θ"["i"], " quantiles"))
ylabel <- expression(paste("Observed ", "θ"["i"], " quantiles"))
qqnorm(data$betaT, col=data$color, cex=1.3, las=1, type="n",
       xlab=xlabel, ylab=ylabel, main="", ylim=c(-1,1), bty="n", xpd=NA)
for(i in c("grey60", "red", "blue")){
  qqs <- subset(qq, col==i)
  points(x=qqs$x, y=qqs$y, col=qqs$col, xpd=NA, cex=1.2)
}
qqline(data$betaT)

# Legend
legend("topleft", bty="n", pch=1, col=c("blue", "grey60", "red"), pt.cex=1.2,
       legend=c("Positive", "Non-significant", "Negative"), title=expression(bold("SST influence")))

# Off
dev.off()
graphics.off()

