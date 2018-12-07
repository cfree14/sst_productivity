
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(freeR)
library(mblm) # For Thiel-Sen slope
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"

# Read data
data <- read.csv(paste(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.csv", sep="/"), as.is=T)


# Plot data
################################################################################

# Sample size
nspp <- data %>% 
  group_by(species, comm_name) %>% 
  summarize(n=n()) %>% 
  arrange(desc(n)) %>% 
  filter(n>=5) %>% 
  ungroup() %>% 
  mutate(a=NA, b=NA, c=NA,
         comm_name=revalue(comm_name, c("Pollock"="Saithe",
                                        "Herring"="Atlantic herring",
                                        "common European sole"="Common sole",
                                        "Walleye pollock"="Alaska pollock")))

# Setup figure
figname <- "SFig14_thetas_by_thermal_niche_position.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3.5, units="in", res=600)
par(mfrow=c(2,4), mar=c(2, 3, 1.5, 0.5), oma=c(1.5,1,0,0), mgp=c(2.8,0.8,0))

# Loop through species
for(i in 1:nrow(nspp)){
  
  # Subset data
  spp <- nspp$species[i]
  comm_name1 <- nspp$comm_name[i]
  sdata <- subset(data, species==spp)
  
  # Plot data
  plot(betaT ~ sst_c_avg, sdata, bty="n", las=1, xlab="", ylab="", 
       xlim=c(0, 20), ylim=c(-1, 1), 
       cex.axis=0.8, pch=16, col="grey30", main="", cex.main=0.8)
  title(main=comm_name1, line=0.9, cex.main=0.8)
  title(main=bquote(italic(.(spp))), line=0.4, cex.main=0.8)
  
  # Fit and plot Thiel-Sen slope
  tsfit <- mblm(betaT ~ sst_c_avg, sdata, repeated=F)
  pvalue <- roundf(summary(tsfit)$coefficients[2,4],2)
  # lty <- ifelse(pvalue<0.1, 1, 2)
  lty <- 1
  curve(coef(tsfit)[1] + coef(tsfit)[2]*x, from=0, to=20, n=100, add=T, lty=lty)
  text(x=20, y=0.95, pos=2, labels=paste0("n=", nspp$n[i]), cex=0.7)
  # text(x=20, y=0.8, pos=2, labels=bquote("r"^2*"="*.(r2)), cex=0.7)
  # text(x=20, y=0.8, pos=2, labels=paste0("p=", pvalue), cex=0.7)
  
  
  # Fit and plot lm
  # lmfit <- lm(betaT ~ sst_c_avg, sdata)
  # r2 <- format(round(summary(lmfit)$r.squared, 2), nsmall=2)
  # pvalue <- format(round(anova(lmfit)$'Pr(>F)'[1],3), nsmall=3)
  # lty <- ifelse(pvalue<0.1, 1, 2)
  # curve(coef(lmfit)[1] + coef(lmfit)[2]*x, from=0, to=20, n=100, add=T, lty=lty)
  # text(x=20, y=0.95, pos=2, labels=paste0("n=", nspp$n[i]), cex=0.7)
  # text(x=20, y=0.8, pos=2, labels=bquote("r"^2*"="*.(r2)), cex=0.7)
  # text(x=20, y=0.63, pos=2, labels=paste0("p=", pvalue), cex=0.7)
  
  # Add 0 influence line
  lines(x=c(0, 20), y=c(0,0), lty=3, lwd=0.6)
  
}

# Add x-axis label
ylabel <- expression("SST influence (θ"["i"]*")")
mtext("Mean SST (°C)", outer=T, side=1, adj=0.53, line=0, cex=0.8)
mtext(ylabel, outer=T, side=2, adj=0.52, line=-0.5, cex=0.8)

# Off
dev.off()
graphics.off()


