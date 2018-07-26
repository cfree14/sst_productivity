
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
data_f <- results.df
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd, results.wide, results.df)

# Read random effects results
load(paste(datadir, "ramldb_v3.8_pella_0.20_cobe_random_normal.Rdata", sep="/"))
data_r <- results.wide
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd, results.df, results.wide)


# Plot data
################################################################################

# Setup figure
figname <- "Fig1_thetas_fixed_effects.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=4, units="in", res=600)
layout(matrix(data=c(1,2), ncol=2), widths=c(0.65, 0.35))
par(mar=c(3.3,0.5,1,0.5), mgp=c(2.2,0.8,0))

# Fixed effects
###################################################

# Individual rheta estimates
vals <- data_f %>%
  filter(param=="BetaT") %>%
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
  arrange(desc(est))

# Setup empty plot
plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=1,
     xlim=c(-8, 14), ylim=c(1,nrow(vals)), xaxt="n",
     xlab=expression("SST influence (θ"["i"]*")"), ylab="", main="Fixed effects", cex.main=1)
axis(1, at=seq(-8,14,2))

# Add theta estimates (points) and errors (lines)
sapply(1:nrow(vals), function(x) lines(x=c(vals$est_lo[x], vals$est_hi[x]), y=c(x,x), col=vals$lcolor[x], lwd=0.6))
points(vals$est, 1:nrow(vals), pch=16, cex=0.6, col=vals$pcolor)

# Add theta=0 line
lines(x=c(0,0), y=c(1, nrow(vals)), lty=3, col="black", lwd=0.8)

# Add positive/negative influence labels
n_pos <- sum(vals$est_lo>0)
n_neg <- sum(vals$est_hi<0)
n_neutral <- nrow(vals)-n_pos-n_neg
text_pos <- paste0(n_pos, " stocks\n", "positive")
text_neg <- paste0(n_neg, " stocks\n", "negative")
text(x=-8, y=6, labels=text_pos, pos=4, adj=1, cex=0.9, col="blue")
text(x=12, y=nrow(vals)-8, labels=text_neg, pos=2, adj=0, cex=0.9, col="red")

# Random effects
###################################################

# Format data
data_r <- data_r %>% 
  arrange(desc(betaT)) %>% 
  mutate(lcolor=revalue(betaT_inf, c("negative"="red", "positive"="blue", "none"="grey60")),
         pcolor=revalue(betaT_inf, c("negative"="red", "positive"="blue", "none"="black")))

# Setup empty plot
plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=1,
     xlim=c(-2, 2), ylim=c(1,nrow(data_r)), main="Random effects", 
     xlab=expression("SST influence (θ"["i"]*")"), ylab="", cex.main=1)

# Add theta estimates (points) and errors (lines)
sapply(1:nrow(data_r), function(x) lines(x=c(data_r$betaT_lo[x], data_r$betaT_hi[x]), y=c(x,x), col=data_r$lcolor[x], lwd=0.6))
points(data_r$betaT, 1:nrow(data_r), pch=16, cex=0.6, col=data_r$pcolor)

# Add theta=0 line
lines(x=c(0,0), y=c(1, nrow(data_r)), lty=3, col="black", lwd=0.8)

# Add positive/negative influence labels
n_pos <- sum(data_r$betaT_lo>0)
n_neg <- sum(data_r$betaT_hi<0)
n_neutral <- nrow(data_r)-n_pos-n_neg
text_pos <- paste0(n_pos, " stocks\n", "positive")
text_neg <- paste0(n_neg, " stocks\n", "negative")
text(x=-2, y=6, labels=text_pos, pos=4, adj=1, cex=0.9, col="blue")
text(x=2, y=nrow(data_r)-8, labels=text_neg, pos=2, adj=0, cex=0.9, col="red")

# Off
dev.off()
graphics.off()


