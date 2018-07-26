
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Read data
data_orig <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)


# Plot data
################################################################################

# Format data
table(data_orig$habitat)
pel.habitats <- c("pelagic-oceanic", "pelagic-neritic")
data <- data_orig %>% 
  select(stockid, betaT, betaT_inf, habitat) %>% 
  mutate(habitat1=ifelse(habitat%in%pel.habitats, "pelagic", "demersal"),
         habitat=factor(habitat, 
                           levels=c("benthic", "demersal", "bathydemersal",
                                    "bathypelagic", "benthopelagic", "reef-associated", 
                                    "pelagic-neritic", "pelagic-oceanic")))

# Sample size
n_hab <- data %>% 
  group_by(habitat) %>% 
  summarize(n=n(),
            n_pos=sum(betaT_inf=="positive"),
            n_neg=sum(betaT_inf=="negative")) %>% 
  filter(!is.na(habitat))
n_hab1 <- data %>% 
  group_by(habitat1) %>% 
  summarize(n=n(),
            n_pos=sum(betaT_inf=="positive"),
            n_neg=sum(betaT_inf=="negative")) %>% 
  filter(!is.na(habitat1))

# Setup figure
figname <- "SFig15_thetas_by_habitat.png"
png(paste(plotdir, figname, sep="/"), width=6, height=4, units="in", res=600)
layout(matrix(1:2, ncol=2), widths=c(0.8,0.2))
par(mar=c(6.5, 0.5, 0.8, 0.5), mgp=c(2.5,0.8,0), oma=c(0,4,0,0), xpd=NA)

# Boxplot by habitat
boxplot(betaT ~ habitat, data, las=2, frame=F, ylim=c(-1.5, 1.5),
        ylab=expression("SST influence (θ"["i"]*")"), cex.axis=0.9,
        col=tcolor(c(rep("brown4", 6), rep("navy", 2)), 0.7), lty=1)
lines(x=c(0.5,8.5), y=c(0, 0), lty=3)
text(x=1:8, y=-1.4, labels=n_hab$n, cex=0.8)
text(x=1:8, y=1.7, labels=n_hab$n_pos, cex=0.8, col=ifelse(n_hab$n_pos>0, "blue", "white"))
text(x=1:8, y=1.5, labels=n_hab$n_neg, cex=0.8, col=ifelse(n_hab$n_neg>0, "red", "white"))
text(x=-0.80, y=1.63, pos=2, labels="A", font=2, cex=1.2, xpd=NA)

# Boxplot by habitat
boxplot(betaT ~ habitat1, data, las=2, frame=F, ylim=c(-1.5, 1.5), lty=1,
        ylab="", yaxt="n", cex.axis=0.9, col=tcolor(c("brown4", "navy"), 0.7))
lines(x=c(0.5,2.5), y=c(0, 0), lty=3)
text(x=1:2, y=-1.4, labels=n_hab1$n, cex=0.8)
text(x=1:2, y=1.75, labels=n_hab1$n_pos, cex=0.8, col=ifelse(n_hab1$n_pos>0, "blue", "white"))
text(x=1:2, y=1.55, labels=n_hab1$n_neg, cex=0.8, col=ifelse(n_hab1$n_neg>0, "red", "white"))
text(x=0.5, y=1.63, pos=2, labels="B", font=2, cex=1.2, xpd=NA)

# Off
dev.off()


