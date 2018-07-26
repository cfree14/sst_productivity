
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output/fixed_effects"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures/fixed_effects"

# Read fixed effects results
load(paste(datadir, "ramldb_v3.8_pella_0.20_cobe_fixed.Rdata", sep="/"))
rm(hess, input.data, model, nstocks, stocks, output, params, problem.stocks, sd, results.df)
stocks <- results.wide
rm(results.wide)


# Build data
################################################################################

# Parameters
p <- 0.2

# Big stocks
stocks1 <- stocks %>%
  filter(abs(betaT)>2) %>%
  arrange(desc(betaT))



# Plot data
################################################################################

# Setup figure
figname <- "Fig3_stocks_w_lg_thetas.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6.5, units="in", res=600)
par(mfrow=c(3,3), mar=c(3.5, 3.5, 0.5, 0.5), oma=c(1,2,1,0))

# Loop through stocks
for(i in 1:nrow(stocks1)){

  # Subset data
  stock <- stocks1$stockid[i]
  r <- stocks1$r[i]
  b0 <- stocks1$k[i]
  betaT <- stocks1$betaT[i]
  betaT_inf <- stocks1$betaT_inf[i]
  sdata <- subset(data, stockid==stock)
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
       xlim=c(0,1), ylim=c(ymin, ymax), xlab="", ylab="")
  title(main=stock, line=0.5, xpd=NA)
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

}

# Axis labels
mtext("Standardized biomass", outer=1, side=1, adj=0.55, line=-0.5)
mtext("Standardized production", outer=1, side=2, adj=0.54)


# Off
dev.off()
graphics.off()


