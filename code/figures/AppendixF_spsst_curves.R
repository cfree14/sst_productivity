
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
inputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/input"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Read TMB SP-SST-model fits
load(paste(outputdir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
stocks <- results.wide
rm(hess, results.df, results.wide, input.data, model, nstocks, output, params, problem.stocks, sd)


# Plot data
################################################################################

# Shape parameter
p <- 0.2

# For y-axis label
top.i <- seq(1, nrow(stocks), 24)

# Setup figure
figname <- "AppendixF_spsst_curves.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 4), mar=c(2.5, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(4,6,3,3), lwd=0.8)

# Loop through problem stocks
for(i in 1:nrow(stocks)){

  # Subset data
  stock <- stocks$stockid[i]
  sdata <- subset(data, stockid==stock)
  r <- stocks$r[i]
  b0 <- stocks$k[i]
  betaT <- stocks$betaT[i]
  betaT_inf <- stocks$betaT_inf[i]
  r_format <- round(r,2)
  b0_format <- round(b0, 1)
  betaT_format <- format(round(betaT, 2), nsmall=2)
  
  # Add color bins to data
  sdata$sst_bin <- cut(sdata$cobe_sst_c_sd, breaks=seq(-1.5,1.5,0.1))
  pcolors <- tcolor(rev(colorpal(brewer.pal(11,"RdBu"),nlevels(sdata$sst_bin))), 0.7)
  sdata$sst_col <- pcolors[sdata$sst_bin]

  # Plot surplus production curve
  ymin <- floor1(min(sdata$sp_sd), 0.05)
  ymax <- ceiling1(max(sdata$sp_sd), 0.05)
  plot(sp_sd ~ tb_sd, sdata, bty="n", las=1, bg=sdata$sst_col, pch=21, cex=1.3,
       xlim=c(0,1), ylim=c(ymin, ymax), xlab="", ylab="", main=stock)
  curve(r/p*x*(1-(x/b0)^p),
        from=0, to=1, add=T, lty=1, lwd=1.4, col="grey20")

  # Plot surplus produciton curves with change temp
  ssts <- c(-1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1)
  colors <- rev(brewer.pal(length(ssts), "RdBu"))
  for(j in 1:length(ssts)){
    curve(r/p*x*(1-(x/b0)^p)*exp(betaT*ssts[j]),
          from=0, to=1, add=T, lty=1, lwd=1.2, col=colors[j])
  }

  # Print theta
  col <- ifelse(betaT<0, "red", "blue")
  col <- ifelse(betaT_inf=="none", tcolor(col, 0.6), col)
  mtext(betaT_format, side=3, adj=0.1, line=-1, cex=0.7, col=col, font=2)

  # Add y-axis label each new page
  if(i%in%top.i){mtext("Standardized biomass", outer=T, side=1, adj=0.5, line=0.6, cex=1)}
  if(i%in%top.i){mtext("Standardized production", outer=T, side=2, adj=0.5, line=0.8, cex=1)}

}


# Off
dev.off()
graphics.off()





