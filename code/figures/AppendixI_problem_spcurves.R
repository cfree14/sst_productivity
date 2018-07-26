
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)

# Directories
inputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Read data
load(paste(inputdir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))
data_orig <- data; rm(data, spfits)

# Read results (to get problem stocks)
load(paste(outputdir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
probs <- problem.stocks
rm(hess, results.df, results.wide, input.data, model, nstocks, output, params, problem.stocks, sd, data, stocks)

# Subset problem stocks
data <- data_orig %>% 
  filter(stockid %in% probs)

# Reduce prob stocks to only those present in data
probs1 <- sort(probs[probs%in%data$stockid])


# Helper functions
################################################################################

# Function to fit surplus production model
# sp <- subset(data, assessid==unique(data$assessid)[1])$sp
# tb <- subset(data, assessid==unique(data$assessid)[1])$tb
fit.sp.model <- function(sp, tb){
  r_start <- log(0.4)
  k_start <- log(max(tb) * 1.5)
  spfit <- try(nls(sp ~ exp(r)*tb*(1-tb/exp(k)),
                   start=list(r=r_start, k=k_start)))
  return(spfit)
}


# Plot problem stocks
################################################################################

# For y-axis label
top.i <- seq(1, length(probs), 24)

# Setup figure
figname <- "AppendixI_problem_spcurves.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 4), mar=c(2.5, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(4,6,3,3), lwd=0.8)

# Loop through problem stocks
for(i in 1:length(probs1)){
  
  # Subset data
  stock <- probs1[i]
  sdata <- subset(data, stockid==stock)
 
  # Plot surplus production curve
  ymin <- floor(min(sdata$sp_sd)/0.05) * 0.05
  ymax <- ceiling(max(sdata$sp_sd)/0.05) * 0.05
  plot(sp_sd ~ tb_sd, sdata, bty="n", xlim=c(0,1), ylim=c(ymin, ymax),
       las=1, xlab="", ylab="", main=stock)
  
  # Fit and plot SP curve if possible
  spfit <- fit.sp.model(sp=sdata$sp_sd, tb=sdata$tb_sd)
  if(!(inherits(spfit, "try-error"))){
    r <- exp(coef(spfit)["r"])
    k <- exp(coef(spfit)["k"])
    curve(r*x*(1-x/k), from=0, to=1, add=T, lty=1, lwd=0.9, col="black")
    rk_text <- paste0("k = ", roundf(k,1), "\nr = ",roundf(r,2))
    text(x=0, y=ymax-(ymax-ymin)*0.1, pos=4, labels=rk_text, cex=1, xpd=NA)
  }
  
  # Add y-axis label each new page
  if(i%in%top.i){mtext("Standardized biomass", outer=T, side=1, adj=0.5, line=0.6, cex=1)}
  if(i%in%top.i){mtext("Standardized production", outer=T, side=2, adj=0.5, line=0.8, cex=1)}
  
}

# Off
dev.off()
graphics.off()





