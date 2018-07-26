
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(forecast)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Real SST simulation fits
load(paste(datadir, "ramldb_v3.8_sst_null_model_fits.Rdata", sep="/"))

# How many didn't match perfectly?
sum(arima_matched==F)

# Plot data
################################################################################

# For y-axis label
top.i <- seq(1, length(stocks), 6)

# Setup figure
figname <- "AppendixE_sst_obs_simulated.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 4), mar=c(1, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,6,3,3), lwd=0.8)

# Loop through stocks
# for(i in 1:12){
for(i in 1:length(stocks)){

  # Subset data
  stock <- stocks[i]
  sdata <- subset(data, stockid==stock)
  print(paste(i, stock))

  # Analyze observed time series
  lmfit <- lm(cobe_sst_c ~ year, sdata)
  arfit <- arfit_list[[i]]
  arfit_n1 <- arfit_n1_list[[i]]

  # Extract statistics
  mu_obs <- mean(sdata$cobe_sst_c)
  sd_obs <- sd(sdata$cobe_sst_c)
  slope_obs <- coef(lmfit)[2]
  ar_pdq_obs <- arimaorder(arfit)

  # Format statistics
  ########################################

  # Format statistics
  mu_obs_f <- format(round(mu_obs,1),nsmall=1)
  sd_obs_f <- format(round(sd_obs,2),nsmall=2)
  slope_obs_f <- format(round(slope_obs*10,2),nsmall=2)
  ar_pdq_obs_f <- paste0("ARIMA(", paste(ar_pdq_obs, collapse=", "), ")")

  # Format N1 statistics
  mu_n1 <- mean(sdata$cobe_sst_c_n1)
  sd_n1 <- sd(sdata$cobe_sst_c_n1)
  lmfit_n1 <- lm(cobe_sst_c_n1 ~ year, sdata)
  slope_n1 <- coef(lmfit_n1)[2]
  mu_n1_f <- format(round(mu_n1,1),nsmall=1)
  sd_n1_f <- format(round(sd_n1,2),nsmall=2)
  slope_n1_f <- format(round(slope_n1*10,2),nsmall=2)
  ar_pdq_n1 <- arimaorder(arfit_n1)
  ar_pdq_n1_f <- paste0("ARIMA(", paste(ar_pdq_n1, collapse=", "), ")")

  # Format N2 statistics
  mu_n2 <- mean(sdata$cobe_sst_c_n2)
  sd_n2 <- sd(sdata$cobe_sst_c_n2)
  lmfit_n2 <- lm(cobe_sst_c_n2 ~ year, sdata)
  slope_n2 <- coef(lmfit_n2)[2]
  mu_n2_f <- format(round(mu_n2,1),nsmall=1)
  sd_n2_f <- format(round(sd_n2,2),nsmall=2)
  slope_n2_f <- format(round(slope_n2*10,2),nsmall=2)
  arfit_n2 <- auto.arima(sdata$cobe_sst_c_n2)
  ar_pdq_n2 <- arimaorder(arfit_n2)
  ar_pdq_n2_f <- paste0("ARIMA(", paste(ar_pdq_n2, collapse=", "), ")")

  # Format N3 statistics
  mu_n3 <- mean(sdata$cobe_sst_c_n3)
  sd_n3 <- sd(sdata$cobe_sst_c_n3)
  mu_n3_f <- format(round(mu_n3,1),nsmall=1)
  sd_n3_f <- format(round(sd_n3,2),nsmall=2)
  lmfit_n3 <- lm(cobe_sst_c_n3 ~ year, sdata)
  slope_n3 <- coef(lmfit_n3)[2]
  slope_n3_f <- format(round(slope_n3*10,2),nsmall=2)
  arfit_n3 <- auto.arima(sdata$cobe_sst_c_n3)
  ar_pdq_n3 <- arimaorder(arfit_n3)
  ar_pdq_n3_f <- paste0("ARIMA(", paste(ar_pdq_n3, collapse=", "), ")")

  # A. Plot SST observed
  ########################################

  # Parameters
  titles <- rep("", 4)
  if(i %in% top.i){titles <- c("Observed", "Null 1", "Null 2", "Null 3")}
  xmin <- min(sdata$year)
  xmax <- max(sdata$year)
  tmin <- min(sdata[,c("cobe_sst_c", "cobe_sst_c_n1", "cobe_sst_c_n2", "cobe_sst_c_n3")])
  tmax <- max(sdata[,c("cobe_sst_c", "cobe_sst_c_n1", "cobe_sst_c_n2", "cobe_sst_c_n3")])
  tdiff <- tmax-tmin
  ymin <- floor((tmin-tdiff*0.3) / 0.1) * 0.1
  ymax <- ceiling(tmax / 0.1) * 0.1

  # Plot data
  plot(cobe_sst_c ~ year, sdata, type="l", bty="n", col="grey50", #las=1,
       xaxt="n", yaxt="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=titles[1])
  axis(1, at=c(xmin, xmax))
  axis(2, at=c(ymin, ymax), las=2)
  lines(x=c(xmin, xmax), y=c(mu_obs, mu_obs), lty=3)
  curve(coef(lmfit)[1]+coef(lmfit)[2]*x, from=xmin, to=xmax, add=T)

  # Add stat text
  stat_text <- paste0(ar_pdq_obs_f, "\n", slope_obs_f, "°C / decade\n", "mu=", mu_obs_f, "; sd=", sd_obs_f)
  text(x=xmax, y=ymin+tdiff*0.15, labels=stat_text, pos=2, cex=0.9)

  # Add stockid
  text(x=xmin, y=ymax, pos=4, labels=stock, font=2, cex=0.9)


  # B. Plot Null 1 simulated
  ########################################

  # Plot data
  plot(cobe_sst_c_n1 ~ year, sdata, type="l", bty="n", col="grey50", las=1,
       xaxt="n", yaxt="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=titles[2])
  axis(1, at=c(xmin, xmax))
  axis(2, at=c(ymin, ymax), las=2)
  lines(x=c(xmin, xmax), y=c(mu_n1, mu_n1), lty=3)
  curve(coef(lmfit_n1)[1]+coef(lmfit_n1)[2]*x, from=xmin, to=xmax, add=T)

  # Add stat text
  stat_text <- paste0(ar_pdq_n1_f, "\n", slope_n1_f, "°C / decade\n", "mu=", mu_n1_f, "; sd=", sd_n1_f)
  text(x=xmax, y=ymin+tdiff*0.15, labels=stat_text, pos=2, cex=0.9)

  # Mark not perfectly matched ones
  matched <- arima_matched[i]
  if(matched==F){text(x=xmin, y=ymax-tdiff*0.15, pos=4, labels="*", font=2, col="blue", cex=2)}


  # C. Plot Null 2 simulated
  ########################################

  # Plot data
  plot(cobe_sst_c_n2 ~ year, sdata, type="l", bty="n", col="grey50", las=1,
       xaxt="n", yaxt="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=titles[3])
  axis(1, at=c(xmin, xmax))
  axis(2, at=c(ymin, ymax), las=2)
  lines(x=c(xmin, xmax), y=c(mu_n2, mu_n2), lty=3)
  curve(coef(lmfit_n2)[1]+coef(lmfit_n2)[2]*x, from=xmin, to=xmax, add=T)

  # Add stat text
  stat_text <- paste0(ar_pdq_n2_f, "\n", slope_n2_f, "°C / decade\n", "mu=", mu_n2_f, "; sd=", sd_n2_f)
  text(x=xmax, y=ymin+tdiff*0.15, labels=stat_text, pos=2, cex=0.9)


  # D. Plot Null 3 simulated
  ########################################

  # Plot data
  plot(cobe_sst_c_n3 ~ year, sdata, type="l", bty="n", col="grey50", las=1,
       xaxt="n", yaxt="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=titles[4])
  axis(1, at=c(xmin, xmax))
  axis(2, at=c(ymin, ymax), las=2)
  lines(x=c(xmin, xmax), y=c(mu_n2, mu_n2), lty=3)
  curve(coef(lmfit_n3)[1]+coef(lmfit_n3)[2]*x, from=xmin, to=xmax, add=T)

  # Add stat text
  stat_text <- paste0(ar_pdq_n3_f, "\n", slope_n3_f, "°C / decade\n", "mu=", mu_n3_f, "; sd=", sd_n3_f)
  text(x=xmax, y=ymin+tdiff*0.15, labels=stat_text, pos=2, cex=0.9)


  # Final touches
  ########################################

  # Add y-axis label
  if(i%in%top.i){mtext("SST (°C)", outer=T, side=2, adj=0.5, line=0.3)}

}

# Off
dev.off()
graphics.off()

