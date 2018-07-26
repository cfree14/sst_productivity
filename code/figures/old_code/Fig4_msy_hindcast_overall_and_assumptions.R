
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"

# Read model data
load(paste(datadir, "ramldb_v3.8_spsst_cobe.Rdata", sep="/"))
spdata <- data
rm(data, hess, input.data, model, nstocks, output, params, problem.stocks, sd, stocks, results.df, results.wide)

# Read results
data <- read.csv(paste(datadir, "ramldb_v3.8_spsst_cobe.csv", sep="/"), as.is=T)

# SST time series data
sst <- read.csv(paste(sstdir, "ramldb_sst_yearly_cobe.csv", sep="/"), as.is=T)


# Examine SST data
################################################################################

# Min/max sst experience
sst_exp <- spdata %>% 
  group_by(assessid) %>% 
  summarize(sst_in_min=min(cobe_sst_c),
            sst_in_max=max(cobe_sst_c))

# Merge full SST data with used SST data
assessids <- sort(unique(data$assessid))
sst_type <- sst %>%
  # Reducee to stocks used in analysis
  filter(assessid%in%assessids & year<2016) %>% 
  # Identify years used in "model" or are "outside" model
  left_join(select(spdata, assessid, year, use), by=c("assessid", "year")) %>% 
  rename(model=use) %>% 
  mutate(model=ifelse(is.na(model), "outside", "model")) %>% 
  # Identify years "outside" model that are "inside" model SST years or are "cooler" or "warmer"
  group_by(assessid) %>% 
  mutate(model_sst_min=min(sst_c[model=="model"]),
         model_sst_max=max(sst_c[model=="model"])) %>% 
  ungroup() %>% 
  mutate(sst_type=model,
         sst_type=ifelse(sst_c<model_sst_min, "cooler", sst_type), 
         sst_type=ifelse(sst_c>model_sst_max, "warmer", sst_type),
         sst_type=ifelse(sst_c>=model_sst_min & sst_c<=model_sst_max & model!="model", "inside", sst_type),
         sst_type_num=as.numeric(revalue(sst_type, c("model"=1, 
                                                     "inside"=2,
                                                     "cooler"=3,
                                                     "warmer"=4))),
         color=revalue(sst_type, c("model"="black", 
                                   "inside"="grey60",
                                   "cooler"="blue",
                                   "warmer"="red")))
str(sst_type)

# Convert to wide format and add some stuff
sst_type_wide <- dcast(sst_type, assessid ~ year, value.var="sst_type_num")
sst_type_wide1 <- sst_type_wide %>% 
  left_join(select(data, assessid, yr1), by="assessid") %>% 
  select(assessid, yr1, everything()) %>% 
  arrange(yr1)

# Reshape SST type data for plotting
sst_type_mat <- as.matrix(sst_type_wide1[,3:ncol(sst_type_wide1)])
sst_type_mat_t <- t(sst_type_mat)

# Plot data
image(x=1850:2015,
      y=1:ncol(sst_type_mat_t), 
      z=sst_type_mat_t, col=c("black", "grey80", "blue", "red"), 
      xlim=c(1840, 2020), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
axis(1, at=seq(1840,2020,20), las=2)

# Extrapolation percent
hind <- sst_type_mat[,as.character(1940:2015)]
n_model <- sum(hind==1)
n_outside <-  sum(hind!=1)
n_warm <- sum(hind==4)
n_cool <- sum(hind==3)
n_in <- sum(hind==2)
n_cool/n_outside
n_warm/n_outside
(n_cool+n_warm)/n_outside


# Build MSY hindcasts
################################################################################

# Add K and MSY to results 
data <- data %>%
  rename(k_orig = k) %>% 
  mutate(k = k_orig * tb_max,
         msy_avg = (r * k) / 4)
msy_avg <- sum(data$msy_avg)

# Build realized MSY data
stockids <- sort(unique(data$stockid))
assessids <- sort(unique(data$assessid))
msys <- sst %>% 
  # Reduce to stocks in analysis
  filter(assessid %in% assessids & year<2016) %>% 
  # Add stockid/r/k/theta/SSTavg from data
  left_join(select(data, assessid, stockid, r, k, betaT, sst_c_avg), by="assessid") %>% 
  # Scale SST time series and calculate realized MSY
  mutate(sst_c_sd=sst_c-sst_c_avg, 
         msy_real=r*k/4*exp(sst_c_sd*betaT)) %>% 
  # Rearrange columns
  select(assessid, stockid, r, k, betaT, sst_c_avg, year, sst_c, sst_c_sd, msy_real)

# Inspect completeness
apply(msys, 2, function(x) sum(is.na(x)))

# Summarize realized MSY over years
msy_ts <- msys %>% 
  group_by(year) %>% 
  summarize(sst=mean(sst_c),
            msy=sum(msy_real))

# Plot quickly
plot(sst ~ year, msy_ts, type="l")
plot(msy/1000 ~ year, msy_ts, type="b")
abline(h=msy_avg/1000)

# Export MSY time series
write.csv(msys, paste(datadir, "msy_hindcast_time_series.csv", sep="/"), row.names=F)


# Plot MSY hindcasts
################################################################################

# Setup figure
figname <- "AppendixE_msy_hindcasts.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 4), mar=c(1, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,3,3,3), lwd=0.8)

# Loop through stocks: i <- 1
for(i in 1:length(stockids)){
  
  # Subset data
  stock <- stockids[i]
  sdata <- subset(msys, stockid==stock)
  betaT <- round(data$betaT[data$stockid==stock],3)
  
  # Plot MSY hindcast
  ymin <- min(sdata$msy_real/1000)
  ymax <- max(sdata$msy_real/1000)
  plot(msy_real/1000 ~ year, sdata, type="l", bty="n", las=1, 
       col="grey50", lwd=0.4, ylim=c(ymin, ymax),
       xaxt="n", yaxt="n", xlab="", ylab="", main=stock)
  axis(2, at=c(ymin, ymax), label=round(c(ymin, ymax), 1))
  
  # Add BetaT
  text(x=1870, y=ymax-(ymax-ymin)*0.05, pos=4, labels=betaT, font=2, cex=1.2)
  
  # Add SST time series
  par(new = T)
  plot(sst_c ~ year, sdata, type="l", col="red", lty=3,
       bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  
}

# Off
dev.off()
graphics.off()


# Plot manuscript figure
################################################################################

# Setup figure
figname <- "Fig4_msy_hindcast.png"
png(paste(plotdir, figname, sep="/"), width=6, height=5, units="in", res=600)
layout(matrix(1:4, ncol=2, byrow=F), heights=c(0.7,0.3), widths=c(0.6, 0.4))
par(mar=c(1,4,0.5,0.5), mgp=c(2.8,0.7,0), oma=c(2,0,0,0))

# A. MSY over time
######################################

# Plot realized MSY over time
range(msy_ts$msy)/1E6
plot(msy/1E6 ~ year, msy_ts, type="l", bty="n", las=1, lwd=0.9, col="grey20",
     xaxt="n", xlim=c(1840, 2020), xlab="", ylab="Global MSY (millions of mt)")
axis(1, at=seq(1840,2020,20), las=2, labels=F)
mtext("A", side=3, adj=0.05, line=-1.5, cex=1, font=2)

# Label MSY at average temperature
lines(x=c(1840,2020), y=rep(msy_avg/1E6,2), lty=3, col="grey50")
text(x=1846, y=37.6, label="MSY @ SST average", pos=4, cex=0.75, col="grey50")

# Label regimes
# yr_shift <- 1988
# regime1avg <- mean(msy_ts$msy[msy_ts$year <= yr_shift]) / 1E6
# regime2avg <- mean(msy_ts$msy[msy_ts$year > yr_shift]) / 1E6
# lines(x=c(1870,yr_shift), y=rep(regime1avg,2), lwd=5, col=rgb(t(col2rgb("blue"))/255, alpha=0.5))
# lines(x=c(yr_shift,2015), y=rep(regime2avg,2), lwd=5, col=rgb(t(col2rgb("red"))/255, alpha=0.5))

# B. SST over time
######################################

# Plot SST over time
plot(sst ~ year, msy_ts, bty="n", type="l", col="red", las=1,
     xaxt="n", xlim=c(1840, 2020), ylim=c(12.5,14.5), xlab="", ylab="SST (°C)")
axis(1, at=seq(1840,2020,20), las=2)
mtext("B", side=3, adj=0.05, line=-1.5, cex=1, font=2)


# C. Full SST vs. model SST
######################################

# Plot data
par(mar=c(3,0.5,0.5,0.8)) # increase bottom padding, decrease left padding, right padding a bit
image(x=1870:2015,
      y=1:ncol(sst_type_mat_t), 
      z=sst_type_mat_t, col=c("black", "grey95", "blue", "red"), 
      xlim=c(1860, 2020), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
axis(1, at=seq(1860,2020,20), las=2, cex.axis=0.8)
mtext("C", side=3, adj=-0.05, line=-1.5, cex=1, font=2)

# Add vertical line
yr <- 1940
lines(x=c(yr, yr), y=c(1, ncol(sst_type_mat_t)), col="grey40", lty=1, lwd=2)
# lines(x=c(1935, 1935), y=c(1, ncol(sst_type_mat_t)), col="black", lty=1, lwd=0.8)

# D. MSY correlation
######################################

# Axis labels
xlabels <- c("10", 
             expression("10"^2), 
             expression("10"^3),
             expression("10"^4),
             expression("10"^5),
             expression("10"^6),
             expression("10"^7))

# Plot MSY correlation
range(data$msy_avg)
range(data$msy_true, na.rm=T)
par(mgp=c(2,0.7,0), mar=c(1,4,0.5,0.5)) # return padding but change axis spacing
plot(msy_avg ~ msy_true, data, log="xy", bty="n",
     xlim=c(10, 1E7), ylim=c(10, 1E7), xaxt="n", yaxt="n", xpd=NA,
     xlab=expression("MSY"["RAM"]*" (mt)"), ylab=expression("MSY"["est"]*" (mt)"), pch=16, col="grey60")
axis(1, at=c(1E1,1E2,1E3,1E4,1E5,1E6,1E7), labels=xlabels, cex.axis=0.9)
axis(2, at=c(1E1,1E2,1E3,1E4,1E5,1E6,1E7), labels=xlabels, las=2, cex.axis=0.9)
mtext("D", side=3, adj=0.05, line=-1.5, cex=1, font=2)
# Add one-to-one line
lines(x=c(1,1E7), y=c(1,1E7))
# Fit linear regression
lmfit <- lm(msy_avg ~ msy_true, data)
r2 <- format(round(summary(lmfit)$r.squared, 2), nsmall=2)
pvalue <- format(round(anova(lmfit)$'Pr(>F)'[1], 3), nsmall=3)
rmse <- sqrt(mean(lmfit$residuals^2))
# Add sample size text
n <- sum(!is.na(data$msy_true))
n_text <- paste0("n=", n)
r2_text <- bquote("r"^2*"="*.(r2))
mtext(n_text, side=1, adj=0.95, line=-2.1, cex=0.7)
mtext(r2_text, side=1, adj=0.95, line=-1.2, cex=0.7)

# Off
dev.off()
graphics.off()

