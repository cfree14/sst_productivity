
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Read model data
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
spdata <- data; stocks <- results.wide
rm(data, hess, input.data, model, nstocks, output, params, problem.stocks, sd, results.df, results.wide)

# Read data
data_orig <- read.csv(paste(datadir, "msy_hindcast_time_series_using_lme_model.csv", sep="/"), as.is=T)


# Format data
################################################################################

# Hindcast period
yr1 <- 1930
yr2 <- 2010

# Format data
data <- data_orig %>%
  filter(year>=yr1 & year<=yr2) %>%
  mutate(msy_avg=msy_avg/1000,
         msy_real=msy_real/1000,
         msy_real_lo=msy_real_lo/1000,
         msy_real_hi=msy_real_hi/1000)


# Plot MSY hindcasts
################################################################################

# Colors
display.brewer.pal(11, "RdBu")
colors <- brewer.pal(11, "RdBu")

# For y-axis label
top.i <- seq(1, nrow(stocks), 24)

# Setup figure
figname <- "AppendixG_msy_hindcasts.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 4), mar=c(1, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,6,3,3), lwd=0.8)

# Loop through stocks: i <- 1
for(i in 1:nrow(stocks)){

  # Subset data
  stock <- stocks$stockid[i]
  sdata <- subset(data, stockid==stock)
  betaT <- format(round(unique(stocks$betaT[i]),2), nsmall=2)
  betaT_sig <- stocks$betaT_inf[i]
  msy_avg <- unique(sdata$msy_avg)
  if(length(msy_avg)>1){
    print(paste(i, msy_avg))
    msy_avg <- unique(sdata$msy_avg)[1] # WHAT? WHY DOES UNIQUE GIVE TWO?
  }

  # Colors
  if(betaT_sig=="none"){
    scolor <- ifelse(betaT<0, colors[5], colors[7])
    lcolor <- ifelse(betaT<0, colors[3], colors[9])
  }else{
    scolor <- ifelse(betaT<0, colors[3], colors[9])
    lcolor <- ifelse(betaT<0, colors[1], colors[11])
  }

  # Plot MSY hindcast
  ymin <- min(sdata$msy_real_lo)
  ymax <- max(sdata$msy_real_hi)
  plot(msy_real ~ year, sdata, type="n", bty="n", las=1,
       col="grey50", lwd=0.4, ylim=c(ymin, ymax),
       xaxt="n", yaxt="n", xlab="", ylab="", main=stock)
  axis(1, at=seq(yr1, yr2, 10), label=F, cex.axis=1.2)
  axis(2, at=c(ymin, ymax), label=round(c(ymin, ymax), 1), cex.axis=1.2)

  # Plot CI shading
  polygon(x=c(sdata$year, rev(sdata$year)),
          y=c(sdata$msy_real_lo, rev(sdata$msy_real_hi)), col=scolor, border=F)

  # Add MSY line
  lines(sdata$year, sdata$msy_real, col=lcolor)

  # Add mean MST line
  lines(x=c(min(sdata$year), max(sdata$year)), y=rep(msy_avg,2), lty=3, lwd=0.7)

  # Add BetaT
  text(x=yr1, y=ymax-(ymax-ymin)*0.05, pos=4, labels=betaT, cex=1.2, col=lcolor)

  # Add line showing model coverge
  spdata1 <- subset(spdata, stockid==stock)
  model_yr1 <- min(spdata1$year)
  model_yr2 <- max(spdata1$year)
  lines(x=c(model_yr1,model_yr2), y=rep(ymin,2), lwd=1.8)

  # Add SST time series
  par(new = T)
  tmin <- min(sdata$sst_c)
  tmax <- max(sdata$sst_c)
  trange <- tmax - tmin
  ymax <- tmax + trange*0.2
  ymin <- tmin - trange*0.7
  tcolor <- rgb(t(col2rgb("grey30"))/255, alpha=0.4)
  plot(sst_c ~ year, sdata, type="l", col=tcolor, lty=1, lwd=0.5,
       bty="n", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(ymin, ymax))

  # Add y-axis label
  if(i%in%top.i){mtext("MSY (1000s mt)", outer=T, side=2, adj=0.5, line=0.3)}

}

# Off
dev.off()
graphics.off()

