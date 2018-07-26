

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(irr)
library(freeR)
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/code"

# Read helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# Read SP-SST results using 4 PT shape parameters

# Read Schaefer (p=1)
load(paste(datadir, "ramldb_v3.8_spsst_cobe.Rdata", sep="/"))
p50 <- results.wide
p50_mu <- format_mu(results.df)
rm(data, hess, input.data, model, output, results.df, results.wide, sd, nstocks, params, problem.stocks, stocks)

# Read MSY @ 45% (p=0.55)
load(paste(datadir, "ramldb_v3.8_spsst_pella_0.55_cobe.Rdata", sep="/"))
p45  <- results.wide
p45_mu <- format_mu(results.df)
rm(data, hess, input.data, model, output, results.df, results.wide, sd, nstocks, params, problem.stocks, stocks)

# Read MSY @ 40% (p=0.20)
load(paste(datadir, "ramldb_v3.8_spsst_pella_0.20_cobe.Rdata", sep="/"))
p40  <- results.wide
p40_mu <- format_mu(results.df)
rm(data, hess, input.data, model, output, results.df, results.wide, sd, nstocks, params, problem.stocks, stocks)

# Read MSY @ 37% (p=0.01)
load(paste(datadir, "ramldb_v3.8_spsst_pella_0.01_cobe.Rdata", sep="/"))
p37  <- results.wide
p37_mu <- format_mu(results.df)
rm(data, hess, input.data, model, output, results.df, results.wide, sd, nstocks, params, problem.stocks, stocks)


# Build data
################################################################################

# Build data
data <- p50 %>% 
  select(stockid, betaT, betaT_inf) %>% 
  rename(betaT_p50=betaT, betaT_inf_p50=betaT_inf) %>% 
  # Add 45% results
  left_join(select(p45, stockid, betaT, betaT_inf), by="stockid") %>% 
  rename(betaT_p45=betaT, betaT_inf_p45=betaT_inf) %>% 
  # Add 40% results
  left_join(select(p40, stockid, betaT, betaT_inf), by="stockid") %>% 
  rename(betaT_p40=betaT, betaT_inf_p40=betaT_inf) %>% 
  # Add 37% results
  left_join(select(p37, stockid, betaT, betaT_inf), by="stockid") %>% 
  rename(betaT_p37=betaT, betaT_inf_p37=betaT_inf)


# Plot data
################################################################################

# Setup figure
figname <- "SFig2_thetas_by_pt_shape.png"
png(paste(plotdir, figname, sep="/"), width=6.6, height=4.5, units="in", res=600)
layout(matrix(1:8, ncol=4, byrow=T), heights=c(0.7, 0.3))
par(mar=c(4,1,1,0.5), mgp=c(2.5, 0.8, 0))

# Loop through models and plot thetas
results <- list(p50, p45, p40, p37)
mus <- list(p50_mu, p45_mu, p40_mu, p37_mu)
model_names <- c("MSY @ 50% K", "MSY @ 45% K", "MSY @ 40% K", "MSY @ 37% K")
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
  plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=0.75,
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
  text(x=1.65, y=nrow(vals)-20, labels=text_neg, pos=2, adj=0, cex=0.9, col="red")
  
}

# Add x-axis label
mtext(expression("SST influence (θ"["i"]*")"), 
      outer=T, side=1, adj=0.53, line=-11.8, cex=0.8)

# New par
par(mar=c(3.5,3.2,0.0,0.5), mgp=c(2.2, 0.8, 0), xpd=NA)

# Empty plot
plot.new()

# Plot 50% vs. 45%
plot(betaT_p45 ~ betaT_p50, data, bty='n', las=1, col="grey50",
     xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), cex.axis=0.75, cex.lab=1.3,
     xlab=expression("θ"["50%"]),  ylab=expression("θ"["45%"]))
lines(x=c(-1.5,1.5), y=c(-1.5, 1.5))
lines(x=c(-1.5,1.5), y=c(0, 0), lty=3, col="grey40")
lines(y=c(-1.5,1.5), x=c(0, 0), lty=3, col="grey40")
perc <- agree(select(data, betaT_inf_p50, betaT_inf_p45))
r2_format <- roundf(r2(lm(betaT_p45 ~ betaT_p50, data)), 2)
text(x=1.68, y=-1.3, pos=2, labels=bquote("r"^"2"*" = "*.(r2_format)), cex=0.8)
text(x=1.68, y=-0.9, pos=2, labels=paste0(round(perc$value,1), "% match"), cex=0.8)

# Plot 50% vs. 40%
plot(betaT_p40 ~ betaT_p50, data, bty='n', las=1, col="grey50",
     xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), cex.axis=0.75, cex.lab=1.3,
     xlab=expression("θ"["50%"]),  ylab=expression("θ"["40%"]))
lines(x=c(-1.5,1.5), y=c(-1.5, 1.5))
lines(x=c(-1.5,1.5), y=c(0, 0), lty=3, col="grey40")
lines(y=c(-1.5,1.5), x=c(0, 0), lty=3, col="grey40")
perc <- agree(select(data, betaT_inf_p50, betaT_inf_p40))
r2_format <- roundf(r2(lm(betaT_p40 ~ betaT_p50, data)), 2)
text(x=1.68, y=-1.3, pos=2, labels=bquote("r"^"2"*" = "*.(r2_format)), cex=0.8)
text(x=1.68, y=-0.9, pos=2, labels=paste0(round(perc$value,1), "% match"), cex=0.8)

# Plot 50% vs. 37%
plot(betaT_p37 ~ betaT_p50, data, bty='n', las=1, col="grey50",
     xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), cex.axis=0.75, cex.lab=1.3,
     xlab=expression("θ"["50%"]),  ylab=expression("θ"["37%"]))
lines(x=c(-1.5,1.5), y=c(-1.5, 1.5))
lines(x=c(-1.5,1.5), y=c(0, 0), lty=3, col="grey40")
lines(y=c(-1.5,1.5), x=c(0, 0), lty=3, col="grey40")
perc <- agree(select(data, betaT_inf_p50, betaT_inf_p37))
r2_format <- roundf(r2(lm(betaT_p37 ~ betaT_p50, data)), 2)
text(x=1.68, y=-1.3, pos=2, labels=bquote("r"^"2"*" = "*.(r2_format)), cex=0.8)
text(x=1.68, y=-0.9, pos=2, labels=paste0(round(perc$value,1), "% match"), cex=0.8)

# Off
dev.off()
graphics.off()



