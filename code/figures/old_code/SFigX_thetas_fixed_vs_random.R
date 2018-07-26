
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"

# Read fixed effects results
fixed <- read.csv(paste(datadir, "ramldb_v3.8_spmodel_nls.csv", sep="/"), as.is=T)

# Read random effects results
load(paste(datadir, "ramldb_v3.8_spmodel_temp_dpdt_sst_yr_t.Rdata", sep="/"))
rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)

# Build data
################################################################################

# Global theta estimate
mu <- results.df %>% 
  filter(param%in%c("mu_T", "sd_T")) %>% 
  rename(est=estimate, 
         se=stderror) %>% 
  mutate(est_lo=est-se*1.96,
         est_hi=est+se*1.96)

# Individual theta estimates
vals <- results.df %>% 
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
  arrange(desc(est)) %>% 
  # Add fixed effects estimates
  left_join(select(fixed, stockid, betaT), by="stockid") %>% 
  rename(est_fixed=betaT)
  
# Plot data
################################################################################

# Setup figure
figname <- "SFig3_thetas_fixed_vs_random.png"
png(paste(plotdir, figname, sep="/"), width=4, height=6, units="in", res=600)
par(mfrow=c(1,1), mar=c(3,0.5,0.5,0.5), mgp=c(2,0.8,0))

# Setup empty plot
xlabel <- expression("SST influence (θ"["i"]*")")
plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=0.75,
     xlim=c(-6, 4), ylim=c(1,nrow(vals)),
     xlab=xlabel, ylab="", cex.main=0.8)

# Add theta mean
mu_est <- mu$est[mu$param=="mu_T"]
mu_min <- mu$est_lo[mu$param=="mu_T"]
mu_max <- mu$est_hi[mu$param=="mu_T"]
rect(xleft=mu_min, xright=mu_max, ytop=nrow(vals), ybottom=1, border=NA, col="grey80")
text(x=0, y=nrow(vals)-2, labels=paste0("μ = ", round(mu_est,2)), pos=4, cex=0.7, col="grey30")

# Add theta estimates (points) and errors (lines)
sapply(1:nrow(vals), function(x) lines(x=c(vals$est_lo[x], vals$est_hi[x]), y=c(x,x), col=vals$lcolor[x], lwd=0.9))
# points(vals$est, 1:nrow(vals), pch=16, cex=0.6, col=vals$pcolor)

# Add theta=0 line
lines(x=c(0,0), y=c(1, nrow(vals)), lty=3, col="black", lwd=0.8)

# Add fixed estimates
points(vals$est_fixed, 1:nrow(vals), pch=17, cex=0.6, col="black")

# # Add positive/negative influence labels
# n_pos <- sum(vals$est_lo>0)
# n_neg <- sum(vals$est_hi<0)
# n_neutral <- nrow(vals)-n_pos-n_neg
# text_pos <- paste0(n_pos, " stocks\n", "positive")
# text_neg <- paste0(n_neg, " stocks\n", "negative")
# text(x=-1.65, y=8, labels=text_pos, pos=4, adj=1, cex=0.7, col="blue")
# text(x=1.65, y=nrow(vals)-8, labels=text_neg, pos=2, adj=0, cex=0.7, col="red")

# Legend
legend("bottomleft", bty="n", pch=c(17,NA), lwd=c(0,1), col=c("black", "grey60"),
       legend=c("Fixed effect estimate", "Random effect estimate"), cex=0.8)

# Off
dev.off()
graphics.off()


