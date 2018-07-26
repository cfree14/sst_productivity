
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
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Load data
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme_n1.Rdata", sep="/"))
rm(hess, results.df, results.wide, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)


# Plot data
################################################################################

# Pick example
nyr <- data %>%
  group_by(assessid, stockid) %>% 
  summarize(n=n()) %>% 
  arrange(desc(n))

# Subset example stock
sort(unique(data$stockid))
sdata <- data %>% 
  filter(stockid=="BLACKROCKNPCOAST") %>% 
  select(assessid, stockid, year, cobe_sst_c, cobe_sst_c_n1, cobe_sst_c_n2, cobe_sst_c_n3) %>% 
  ungroup()

# Calculate stats
ssts <- select(sdata, cobe_sst_c, cobe_sst_c_n1, cobe_sst_c_n2, cobe_sst_c_n3)
mus <- roundf(apply(ssts, 2, mean),1)
sds <- roundf(apply(ssts, 2, sd),2)
trends <- trimws(roundf(sapply(1:4, function(x) coef(lm(as.matrix(ssts[,x]) ~ sdata$year))[2])*10,2))

# Params
cex.axis <- 1.2
cex.lab <- 1.4
colors <- brewer.pal(4, "Set1")
col2rgb(colors)
sst_mu <- mean(sdata$cobe_sst_c)

# Setup figure
figname <- "SFig3_sst_null_model_examples.png"
png(paste(plotdir, figname, sep="/"), width=6, height=6, units="in", res=600)
par(mfrow=c(4,1), mar=c(2,5,0.5,0.5), mgp=c(2.5,0.7,0))

# Plot observed
plot(cobe_sst_c ~ year, sdata, type="l", las=1, bty="n", 
     cex.lab=cex.lab, cex.axis=cex.axis, col=colors[1],
     xlim=c(1960, 2010), ylim=c(10.5,13), xlab="", ylab="")
# axis(2, at=10:13, las=1, cex.axis=cex.axis)
lines(x=c(1960, 2010), y=rep(sst_mu,2), lty=3, lwd=0.8, col="grey30")
text(x=1959, y=12.8, pos=4, labels="A. Observed", cex=1.2, font=2)
stat_text <- paste("μ=", mus[1], "°C | σ=", sds[1], "°C | ", trends[1], "°C / decade", sep="")
text(x=2010, y=10.6, pos=2, cex=1, labels=stat_text)

# Plot Null 1
plot(cobe_sst_c_n1 ~ year, sdata, type="l", las=1, bty="n",
     cex.lab=cex.lab, cex.axis=cex.axis, col=colors[4],
     xlim=c(1960, 2010), ylim=c(10.5,13), xlab="", ylab="")
# axis(2, at=12:15, las=1, cex.axis=cex.axis)
lines(x=c(1880, 2020), y=rep(sst_mu,2), lty=3, lwd=0.8, col="grey30")
text(x=1959, y=12.8, pos=4, labels="B. Null 1 - same μ/σ/AR/trend", cex=1.2, font=2)
stat_text <- paste("μ=", mus[2], "°C | σ=", sds[2], "°C | ", trends[2], "°C / decade", sep="")
text(x=2010, y=10.6, pos=2, cex=1, labels=stat_text)

# Plot Null 2
plot(cobe_sst_c_n2 ~ year, sdata, type="l", las=1, bty="n",
     cex.lab=cex.lab, cex.axis=cex.axis, col=colors[3],
     xlim=c(1960, 2010), ylim=c(10.5,13), xlab="", ylab="")
# axis(2, at=12:15, las=1, cex.axis=cex.axis)
lines(x=c(1880, 2020), y=rep(sst_mu,2), lty=3, lwd=0.8, col="grey30")
text(x=1959, y=12.8, pos=4, labels="C. Null 2 - same μ/σ/AR", cex=1.2, font=2)
stat_text <- paste("μ=", mus[3], "°C | σ=", sds[3], "°C | ", trends[3], "°C / decade", sep="")
text(x=2010, y=10.6, pos=2, cex=1, labels=stat_text)

# Plot Null 3
plot(cobe_sst_c_n3 ~ year, sdata, type="l", las=1, bty="n", 
     cex.lab=cex.lab, cex.axis=cex.axis, col=colors[2],
     xlim=c(1960, 2010), ylim=c(10.5,13), xlab="", ylab="")
# axis(2, at=12:15, las=1, cex.axis=cex.axis)
lines(x=c(1880, 2020), y=rep(sst_mu,2), lty=3, lwd=0.8, col="grey30")
text(x=1959, y=12.8, pos=4, labels="D. Null 3 - same μ/σ", cex=1.2, font=2)
stat_text <- paste("μ=", mus[4], "°C | σ=", sds[4], "°C | ", trends[4], "°C / decade", sep="")
text(x=2010, y=10.6, pos=2, cex=1, labels=stat_text)

# Add y-axis label
mtext("SST (°C)", outer=T, side=2, adj=0.52, line=-1.5, cex=1)

# Off
dev.off()
graphics.off()


