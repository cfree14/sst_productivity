
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"

# Read model data
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
spdata <- data
rm(data, hess, input.data, model, nstocks, output, params, problem.stocks, sd, stocks, results.df, results.wide)

# Read results
data <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)

# SST time series data
sst <- read.csv(paste(sstdir, "ramldb_sst_yearly_cobe.csv", sep="/"), as.is=T)


# Catch data
################################################################################

# Total catch by year
# 2000=catch in all years=27,996,497.951 mt in 2000
catch <- spdata %>% 
  group_by(year) %>% 
  summarise(n=n(),
            tot=sum(catch))


# Build data
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
  filter(assessid%in%assessids & year<=2015) %>% 
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

# Plot data quickly
image(x=1850:2015,
      y=1:ncol(sst_type_mat_t), 
      z=sst_type_mat_t, col=c("black", "grey80", "blue", "red"), 
      xlim=c(1840, 2020), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
axis(1, at=seq(1840,2020,20), las=2)

# Calculate % extrapolated by year
# model, inside, cooler, warmer
ext_perc <- sst_type %>%
  group_by(year) %>% 
  summarize(n=n(), # total
            nout=sum(sst_type!="model"), # n outside model
            n_c=sum(sst_type=="cooler"),
            n_w=sum(sst_type=="warmer"),
            n_cw=sum(sst_type%in%c("cooler", "warmer")), # n cooler/warmer
            p_c_all=n_c/n*100,
            p_w_all=n_w/n*100,
            p_cw_all=n_cw/n*100,
            p_cw_out=n_cw/nout*100)

# Plot % extrapolated by year
plot(p_cw_all ~ year, ext_perc, type="l", bty="n", las=1,
     xaxt="n", xlim=c(1840,2020), ylab="% of years")
axis(1, at=seq(1840,2020,20), las=2)
abline(h=15, lty=2)


# Stats for manuscript
################################################################################

# Get hindcast window
hind_data <- sst_type_mat[, as.character(1930:2010)]
n_outside <- sum(hind_data!=1)
n_warm <- sum(hind_data==4)
n_cool <- sum(hind_data==3)
p_cool <- n_cool / n_outside
p_warm <- n_warm / n_outside
p_out <- (n_cool+n_warm) / n_outside
(n_cool+n_warm) / (ncol(hind_data) * nrow(hind_data)) * 100

# Plot data
################################################################################

# Years: 1940-2010

# Setup figure
figname <- "SFig4_hindcast_extrapolation.png"
png(paste(plotdir, figname, sep="/"), width=6, height=6, units="in", res=600)
layout(matrix(c(1,2,4,3), ncol=2, byrow=F), heights=c(0.6,0.4))
par(mar=c(3, 3,0.5,0.5), mgp=c(2,0.7,0))

# A. Classified SST time series
########################################

# Plot data
image(x=1850:2015,
      y=1:ncol(sst_type_mat_t), 
      z=sst_type_mat_t, col=c("black", "grey95", "blue", "red"), 
      xlim=c(1840, 2020), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
axis(1, at=seq(1840,2020,20), las=2, cex.axis=1)
mtext("A", side=3, adj=-0.05, line=-1.5, cex=1, font=2)

# Add cutoff line
# yr <- 1940
# lines(x=c(yr, yr), y=c(1, ncol(sst_type_mat_t)), col="grey40", lty=1, lwd=2)

# B. % extrapolated by year
########################################

# Setup empty plot
plot(p_cw_all ~ year, ext_perc, type="n", bty="n", las=1, cex.lab=1, 
     xaxt="n", xlim=c(1840,2020), xlab="", ylab="% extrapolated", cex.axis=1)
axis(1, at=seq(1840,2020,20), las=2, cex.axis=1)
mtext("B", side=3, adj=0.05, line=-1.5, cex=1, font=2)

# Plot % cooler polygon
polygon(x=c(ext_perc$year, rev(ext_perc$year)), 
        y=c(rep(0, nrow(ext_perc)), rev(ext_perc$p_c_all)), 
        col=rgb(t(col2rgb("blue"))/255, alpha=0.8), border=F)
polygon(x=c(ext_perc$year, rev(ext_perc$year)), 
        y=c(ext_perc$p_c_all, rev(ext_perc$p_cw_all)), 
        col=rgb(t(col2rgb("red"))/255, alpha=0.8), border=F)

# Add cutoff threshold
lines(x=c(1840, 2020), y=c(15,15), lty=3)

# Add legend
legend(x=1930, y=65,  bty="n", cex=0.8, 
       fill=c("blue", "red", "grey60", "black"), 
       legend=c("Cooler", "Warmer", "Inside model", "Model years"))

# C. MSY comparison
########################################

# Axis labels
xlabels <- c("10", 
             expression("10"^2), 
             expression("10"^3),
             expression("10"^4),
             expression("10"^5),
             expression("10"^6))

# Add K and MSY to results
p <- 0.2
div <- (p+1)^((p+1)/p)
data <- data %>%
  rename(k_orig = k) %>% 
  mutate(k = k_orig * tb_max,
         msy_avg = (r * k) / div)

# Plot MSY correlation
range(data$msy_avg)
range(data$msy_true, na.rm=T)
plot(msy_avg ~ msy_true, data, log="xy", bty="n",
     xlim=c(10, 1E6), ylim=c(10, 1E6), xaxt="n", yaxt="n", xpd=NA,
     xlab=expression("MSY"["RAM"]*" (mt)"), ylab=expression("MSY"["est"]*" (mt)"), pch=16, col="grey60")
axis(1, at=c(1E1,1E2,1E3,1E4,1E5,1E6), labels=xlabels, cex.axis=1)
axis(2, at=c(1E1,1E2,1E3,1E4,1E5,1E6), labels=xlabels, las=2, cex.axis=1)
mtext("C", side=3, adj=0.05, line=-1.5, cex=1, font=2)
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

########################################

# Off
dev.off()
graphics.off()
