
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

