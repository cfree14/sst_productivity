
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


# Plot data
################################################################################

# Setup figure
figname <- "SFig6_thetas_group_models.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6, units="in", res=600)
par(mfcol=c(2,3), mar=c(2,0.5,1,0.5), mgp=c(2.8,0.8,0), oma=c(2,0,0,0))

# Loop through models
model_objs <- paste0("ramldb_v3.8_spsst_pella_cobe_", c("order", "family", "fao_area", "lme", "method", "method1"), ".Rdata")
model_names <- c("Order", "Family", "FAO area", "Ecoregion", "Generic assessment method", "Specific assessment method")
for(i in 1:length(model_objs)){

  # Load data
  model_obj <- model_objs[i]
  load(paste(datadir, model_obj, sep="/"))
  rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)
  
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
    arrange(desc(est))
  
  # Setup empty plot
  plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=1,
       xlim=c(-1.5, 1.5), ylim=c(1,nrow(vals)),
       xlab="", ylab="", main=model_names[i], cex.main=1.1)
  
  # Add theta mean
  mu_est <- mu$est[mu$param=="mu_T"]
  mu_min <- mu$est_lo[mu$param=="mu_T"]
  mu_max <- mu$est_hi[mu$param=="mu_T"]
  rect(xleft=mu_min, xright=mu_max, ytop=nrow(vals), ybottom=1, border=NA, col="grey80")
  text(x=0, y=nrow(vals), labels=paste0("μ = ", format(round(mu_est,2), nsmall=2)), pos=4, cex=0.9, col="grey30", xpd=NA)

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
  text(x=-1.65, y=8, labels=text_pos, pos=4, adj=1, cex=0.95, col="blue")
  text(x=1.65, y=nrow(vals)-18, labels=text_neg, pos=2, adj=0, cex=0.95, col="red")
  
}

# Add x-axis label
xlabel <- expression("SST influence (θ"["i"]*")")
mtext(xlabel, outer=T, side=1, adj=0.49, line=0.5, cex=0.9)

# Off
dev.off()
graphics.off()


