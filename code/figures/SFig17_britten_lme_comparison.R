
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"
tabledir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/tables"
brittendir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/references/Britten_etal_2016_supp"

# Read my data
# load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
# rm(hess, input.data, model, output, sd, nstocks, params, problem.stocks, data, results.wide)
data_orig <- read.csv(paste(tabledir, "STable6_lme_msy_hindcast.csv", sep="/"), as.is=T)

# Read Britten data
britten <- read.csv(paste(brittendir, "Britten_etal_2016_Fig2B_data.csv", sep="/"), as.is=T)


# Build data
################################################################################

# Format Britten data
britten1 <- britten %>%
  select(-number) %>% 
  mutate(lme_name=revalue(lme_name, c("Iberian Coast"="Iberian Coastal",
                                      "Newfoundland-Labrador"="Labrador - Newfoundland",
                                      "Northeast U.S. Shelf"="Northeast U.S. Continental Shelf",
                                      "Southeast U.S. Shelf"="Southeast U.S. Continental Shelf",
                                      "Southwest Australian Shelf"="South West Australian Shelf",
                                      "Celtic Sea"="Celtic-Biscay Shelf",
                                      "Iceland Shelf"="Iceland Shelf and Sea")))
  
# Merge Britten data with ours
data <- data_orig %>% 
  # Subset our data and add Britten's data
  select(lme, n, sst_c_slope_10yr, betaT, ts_slope, ts_slope_sd, pdiff) %>% 
  mutate(betaTmult=sst_c_slope_10yr*betaT) %>% 
  left_join(britten1, by=c("lme"="lme_name")) %>% 
  filter(!is.na(rmax)) %>% 
  # Add slope colors
  mutate(slope_col="grey60",
         slope_col=ifelse(ts_slope_sd<0 & rmax<0, "red", slope_col),
         slope_col=ifelse(ts_slope_sd>0 & rmax>0, "blue", slope_col), 
         pdiff_col="grey60",
         pdiff_col=ifelse(pdiff<0 & rmax<0, "red", pdiff_col),
         pdiff_col=ifelse(pdiff>0 & rmax>0, "blue", pdiff_col),
         betaT_col="grey60",
         betaT_col=ifelse(betaTmult<0 & rmax<0, "red", betaT_col),
         betaT_col=ifelse(betaTmult>0 & rmax>0, "blue", betaT_col))

# Are any of theirs not in our data?
britten1$lme_name[!britten1$lme_name%in%data_orig$lme]

# Are any of ours not in their data?
data_orig$lme[!data_orig$lme%in%britten1$lme_name]

# Plot data
################################################################################

# Setup figure
figname <- "SFig17_britten_lme_comparison.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=2.75, units="in", res=600)
par(mfrow=c(1,3), mar=c(4.7, 3.5, 0.5, 0.5), mgp=c(2.1,0.8,0))

# Parameters
metrics <- c("ts_slope_sd", "pdiff", "betaTmult")
metric_cols <- c("slope_col", "pdiff_col", "betaT_col")
metric_labels <- c(expression("Δ MSY/MSY"["max"]*" per decade"), 
                   " % difference in MSY\n(1930-39 to 2000-10)", 
                   expression("SST influence (θ"["i"]*") x SST trend (°C/10 yr)"))
ylim <- matrix(c(-0.04, 0.04, 0.02,
                 -80, 40, 20,
                 -0.03, 0.02, 0.01), ncol=3, byrow=T)

# Loop through metrics
for(i in 1:length(metrics)){
  
  # Plot percent difference against Rmax
  plot(x=data$rmax, y=data[,metrics[i]], bty="n", pch=16, yaxt="n", cex=1.3,
       xlim=c(-8,4), ylim=c(ylim[i,1], ylim[i,2]),
       xlab="", ylab=metric_labels[i], col=data[,metric_cols[i]], xpd=NA)
  axis(side=2, at=seq(ylim[i,1], ylim[i,2], ylim[i,3]))
  lines(x=c(-8,4), y=c(0,0), lty=3)
  lines(x=c(0,0), y=c(ylim[i,1],ylim[i,2]), lty=3)
  mtext(LETTERS[i], side=3, adj=0.05, line=-2, font=2, cex=0.8)
  
  # Add trend line
  lmfit <- lm(data[,metrics[i]] ~ data$rmax)
  curve(coef(lmfit)[1]+coef(lmfit)[2]*x, from=-8, to=4, add=T)
  
  # Calculate and plot r2 and agreement
  r2 <- cor(data$rmax, data[,metrics[i]])
  pagree <- sum(data[,metric_cols[i]]!="grey60") / nrow(data) * 100
  stext <- paste0("r2 = ", roundf(r2,2), "\n", roundf(pagree,1), "% agree")
  text(x=4, y=ylim[i,1], label=stext, cex=0.75, xpd=NA, pos=2, offset=0.01)
}

# X-axis labels
mtext(expression("Trend in recruitment potential (R"["MAX"]*")"), outer=T, side=1, adj=0.55, line=-2.4, cex=0.7)
mtext("(Britten et al. 2016)", outer=T, side=1, adj=0.55, line=-1.4, cex=0.7, font=3)

# Off
dev.off()
graphics.off()

