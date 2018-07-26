
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"


# Read random effects results
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)

# Read fixed effects results
data_f <- read.csv(paste(datadir, "ramldb_v3.8_sp_nls.csv", sep="/"), as.is=T)


# Build data
################################################################################

# Build data
data <- results.df %>% 
  # Reduce and rename
  filter(param=="BetaT") %>% 
  select(-param) %>% 
  rename(est=estimate, 
         se=stderror) %>% 
  # 95% confidence intervals and colors
  mutate(est_lo=est-se*1.96,
         est_hi=est+se*1.96,
         lcolor="grey60",
         pcolor="black",
         lcolor=ifelse(est_hi<0, "red", lcolor),
         pcolor=ifelse(est_hi<0, "red", pcolor),
         lcolor=ifelse(est_lo>0, "blue", lcolor),
         pcolor=ifelse(est_lo>0, "blue", pcolor)) %>% 
  # Add fixed effects estimates
  left_join(select(data_f, stockid, betaT), by="stockid") %>% 
  rename(est_fixed=betaT) %>% 
  arrange(desc(est))


# Plot data
################################################################################

# Setup figure
figname <- "SFig10_fixed_vs_random_effects.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6.5, units="in", res=600)
layout(matrix(c(1,2,
                1,3,
                1,4), byrow=T, ncol=2), widths=c(0.7,0.3))
par(mar=c(3.3,0.5,1,0.5), mgp=c(2.2,0.8,0))

# A. Spline plot
######################################

# Setup empty plot
range(data$est_fixed, na.rm=T)
xlabel <- expression("SST influence (θ"["i"]*")")
plot(1:10, 1:10, type="n", bty="n", yaxt="n",
     xlim=c(-4, 4), ylim=c(1,nrow(data)),
     xlab=xlabel, ylab="", cex.lab=1.3, cex.axis=1.1)
mtext("A", side=3, adj=0.1, line=-2, font=2)
  
# Add theta estimates (points) and errors (lines)
sapply(1:nrow(data), function(x) lines(x=c(data$est_lo[x], data$est_hi[x]), y=c(x,x), col=data$lcolor[x], lwd=1))
# points(data$est, 1:nrow(data), pch=16, cex=0.6, col=data$pcolor)

# Add theta=0 line
lines(x=c(0,0), y=c(1, nrow(data)), lty=3, col="black", lwd=1.5)
  
# Add fixed effects points
points(data$est_fixed, 1:nrow(data), pch=17, cex=0.8)

# Add legend
legend(y=15, x=-3.8, bty="n", legend=c("Fixed effect", "Random effect"), cex=1.3, pch=c(17, NA), lty=c(NA, 1))


# B. Random effects histogram
######################################

# Random effects
par(xpd=NA)
range(data$est, na.rm=T)
hist(data$est, breaks=seq(-3.4,9.7,0.1), 
     border=NA, col="grey50", las=1, xlim=c(-4,10),
     main="", xlab=expression("θ"["random"]), cex.lab=1.3, cex.axis=1.1)
mtext("B", side=3, adj=0.05, line=-2, font=2)

# C. Fixed effects histogram
######################################

# Fixed effects
range(data$est_fixed, na.rm=T)
hist(data$est_fixed, breaks=seq(-3.4,9.7,0.1), 
     border=NA, col="grey50", las=1, xlim=c(-4,10),
     main="", xlab=expression("θ"["fixed"]), cex.lab=1.3, cex.axis=1.1)
mtext("C", side=3, adj=0.05, line=-2, font=2)

# D.Fixed vs. random effects
######################################

# Fixed vs. random effects
plot(est_fixed ~ est, data, bty="n", xlim=c(-1,1.5), ylim=c(-4,10), las=1, col="grey50",
     xlab=expression("θ"["random"]), ylab=expression("θ"["fixed"]), cex.lab=1.3, cex.axis=1.1)
lines(x=c(0,0), y=c(-4,10), lty=3, col="grey10")
lines(x=c(-1,1.5), y=c(0,0), lty=3, col="grey10")
mtext("D", side=3, adj=0.05, line=-2, font=2)


# Off
dev.off()
graphics.off()


