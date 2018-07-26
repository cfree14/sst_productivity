

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/code"

# Read helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# Read SP-SST results using 3 SST datasets
# COBE (v2), ERSST (v4), HadISST (v1.1)

# Read HadISST version
load(paste(datadir, "ramldb_v3.8_spsst_had.Rdata", sep="/"))
hadsst <- results.wide
had_mu <- format_mu(results.df)
rm(data, hess, input.data, model, output, results.df, results.wide, sd, nstocks, params, problem.stocks, stocks)

# Read COBE version
load(paste(datadir, "ramldb_v3.8_spsst_cobe.Rdata", sep="/"))
cobe <- results.wide
cobe_mu <- format_mu(results.df)
rm(data, hess, input.data, model, output, results.df, results.wide, sd, nstocks, params, problem.stocks, stocks)

# Read ERSST version
load(paste(datadir, "ramldb_v3.8_spsst_er.Rdata", sep="/"))
ersst <- results.wide
er_mu <- format_mu(results.df)
rm(data, hess, input.data, model, output, results.df, results.wide, sd, nstocks, params, problem.stocks, stocks)


# Build data
################################################################################

# Build data
data <- cobe %>% 
  select(stockid, betaT) %>% 
  rename(betaT_cobe=betaT) %>% 
  # Add ERSST results
  left_join(select(ersst, stockid, betaT), by="stockid") %>% 
  rename(betaT_er=betaT) %>% 
  # Add HadISST results
  left_join(select(hadsst, stockid, betaT), by="stockid") %>% 
  rename(betaT_had=betaT)


# Plot data
################################################################################

# Setup figure
figname <- "SFig1_thetas_by_sst_dataset.png"
png(paste(plotdir, figname, sep="/"), width=6.6, height=5, units="in", res=600)
layout(matrix(1:6, ncol=3, byrow=T), heights=c(0.7, 0.3))
par(mar=c(4,1,1,0.5), mgp=c(2.5, 0.8, 0))

# Loop through models and plot thetas
results <- list(cobe, ersst, hadsst)
mus <- list(cobe_mu, er_mu, had_mu)
model_names <- c("COBE v2", "ERSST v4", "HadISST v1.1")
for(i in 1:length(results)){
  
  # Extract and format data
  vals <- results[[i]] 
  vals <- vals %>% 
    mutate(lcolor="grey60",
           pcolor="black",
           lcolor=ifelse(betaT_hi<0, "red", lcolor),
           pcolor=ifelse(betaT_hi<0, "red", pcolor),
           lcolor=ifelse(betaT_lo>0, "blue", lcolor),
           pcolor=ifelse(betaT_lo>0, "blue", pcolor)) %>% 
    arrange(desc(betaT))
  
  # Setup empty plot
  plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=1,
       xlim=c(-1.5, 1.5), ylim=c(1,nrow(vals)),
       xlab="", ylab="", main=model_names[i], cex.main=1.2)
  
  # Add theta mean
  mu <- mus[[i]]
  mu_est <- mu$est[mu$param=="mu_T"]
  mu_min <- mu$est_lo[mu$param=="mu_T"]
  mu_max <- mu$est_hi[mu$param=="mu_T"]
  rect(xleft=mu_min, xright=mu_max, ytop=nrow(vals), ybottom=1, border=NA, col="grey80")
  text(x=0, y=nrow(vals), labels=paste0("μ = ", format(round(mu_est,2), nsmall=2)), pos=4, cex=0.9, col="grey30")
  
  # Add theta estimates (points) and errors (lines)
  sapply(1:nrow(vals), function(x) lines(x=c(vals$betaT_lo[x], vals$betaT_hi[x]), y=c(x,x), col=vals$lcolor[x], lwd=0.6))
  points(vals$betaT, 1:nrow(vals), pch=16, cex=0.6, col=vals$pcolor)
  
  # Add theta=0 line
  lines(x=c(0,0), y=c(1, nrow(vals)), lty=3, col="black", lwd=0.8)
  
  # Add positive/negative influence labels
  n_pos <- sum(vals$betaT_lo>0)
  n_neg <- sum(vals$betaT_hi<0)
  n_neutral <- nrow(vals)-n_pos-n_neg
  text_pos <- paste0(n_pos, " stocks\n", "positive")
  text_neg <- paste0(n_neg, " stocks\n", "negative")
  text(x=-1.65, y=8, labels=text_pos, pos=4, adj=1, cex=0.9, col="blue")
  text(x=1.65, y=nrow(vals)-8, labels=text_neg, pos=2, adj=0, cex=0.9, col="red")
  
}

# Add x-axis label
mtext(expression("SST influence (θ"["i"]*")"), 
      outer=T, side=1, adj=0.53, line=-12.8, cex=0.8)

# New par
par(mar=c(4,3.8,0.2,0.5), mgp=c(2.5, 0.8, 0))

# Plot COBE vs. HadISST
plot(betaT_cobe ~ betaT_had, data, bty='n', las=1, col="grey50",
     xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), cex.axis=0.87, cex.lab=1.3,
     xlab=expression("θ"["HadISST"]),  ylab=expression("θ"["COBE"]))
lines(x=c(-1.5,1.5), y=c(-1.5, 1.5))
lines(x=c(-1.5,1.5), y=c(0, 0), lty=3, col="grey40")
lines(y=c(-1.5,1.5), x=c(0, 0), lty=3, col="grey40")
r2_format <- roundf(r2(lm(betaT_cobe ~ betaT_had, data)), 2)
text(x=1.6, y=-1.3, pos=2, labels=bquote("r"^"2"*" = "*.(r2_format)), cex=0.9)

# Plot COBE vs. ERSST
plot(betaT_cobe ~ betaT_er, data, bty='n', las=1, col="grey50",
     xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), cex.axis=0.87, cex.lab=1.3,
     xlab=expression("θ"["ERSST"]),  ylab=expression("θ"["COBE"]))
lines(x=c(-1.5,1.5), y=c(-1.5, 1.5))
lines(x=c(-1.5,1.5), y=c(0, 0), lty=3, col="grey40")
lines(y=c(-1.5,1.5), x=c(0, 0), lty=3, col="grey40")
r2_format <- roundf(r2(lm(betaT_cobe ~ betaT_er, data)), 2)
text(x=1.6, y=-1.3, pos=2, labels=bquote("r"^"2"*" = "*.(r2_format)), cex=0.9)

# Plot ERSST vs. HadISST
plot(betaT_er ~ betaT_had, data, bty='n', las=1, col="grey50",
     xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), cex.axis=0.87, cex.lab=1.3,
     xlab=expression("θ"["HadISST"]),  ylab=expression("θ"["ERSST"]))
lines(x=c(-1.5,1.5), y=c(-1.5, 1.5))
lines(x=c(-1.5,1.5), y=c(0, 0), lty=3, col="grey40")
lines(y=c(-1.5,1.5), x=c(0, 0), lty=3, col="grey40")
r2_format <- roundf(r2(lm(betaT_er ~ betaT_had, data)), 2)
text(x=1.6, y=-1.3, pos=2, labels=bquote("r"^"2"*" = "*.(r2_format)), cex=0.9)

# Off
dev.off()
graphics.off()



