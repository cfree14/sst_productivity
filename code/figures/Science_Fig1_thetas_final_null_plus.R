
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"

# Plot data
################################################################################

# Setup figure
figname <- "Science_Fig1_thetas_final_null_plus.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=4, units="in", res=600)
layout(matrix(c(1,2,3,
                1,2,4,
                1,2,5), ncol=3, byrow=T), widths=c(0.70/2, 0.70/2, 0.3))
par(mar=c(1,0.5,1,0.5), mgp=c(2.8,0.8,0), oma=c(3,0,0,0))

# SST influence splines
#######################################

# Loop through models
models <- c("ramldb_v3.8_spsst_pella_cobe_lme.Rdata",
            "ramldb_v3.8_spsst_pella_cobe_lme_n1.Rdata")
model_names <- c("Final model", "Null model")
for(i in 1:length(models)){

  # Load data
  load(paste(datadir, models[i], sep="/"))
  rm(hess, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)
  
  # Global theta estimate
  mu <- results.df %>% 
    filter(param%in%c("mu_T", "sd_T")) %>% 
    rename(est=estimate, 
           se=stderror) %>% 
    mutate(est_lo=est-se*1.96,
           est_hi=est+se*1.96)
  
  # Individual rheta estimates
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
       xlab="", ylab="", main=model_names[i], cex.main=1.2)
  title(LETTERS[i], cex.main=1.2, adj=0)
  # mtext(LETTERS[i], side=3, line=-0.1, adj=0, font=2, cex=0.9)
  
  # Add theta mean
  mu_est <- mu$est[mu$param=="mu_T"]
  mu_min <- mu$est_lo[mu$param=="mu_T"]
  mu_max <- mu$est_hi[mu$param=="mu_T"]
  rect(xleft=mu_min, xright=mu_max, ytop=nrow(vals), ybottom=1, border=NA, col="grey80")
  text(x=0, y=nrow(vals), labels=paste0("μ = ", format(round(mu_est,2), nsmall=2)), pos=4, cex=0.9, col="grey30")

  # Add theta estimates (points) and errors (lines)
  sapply(1:nrow(vals), function(x) lines(x=c(vals$est_lo[x], vals$est_hi[x]), y=c(x,x), col=vals$lcolor[x], lwd=0.6))
  points(vals$est, 1:nrow(vals), pch=16, cex=0.6, col=vals$pcolor)
  
  # Add theta=0 line
  lines(x=c(0,0), y=c(1, nrow(vals)), lty=3, col="black", lwd=0.8)
  
  # Add positive/negative influence labels
  n_pos <- sum(vals$est_lo>0)
  n_neg <- sum(vals$est_hi<0)
  n_neutral <- nrow(vals)-n_pos-n_neg
  if(i==1){
    text_pos <- paste0(n_pos, " stocks\n", "positively\n", "influenced\n", "by warming")
    text_neg <- paste0(n_neg, " stocks\n", "negatively\n", "influenced\n", "by warming")
  }else{
    text_pos <- paste0("\n \n \n", n_pos, " stocks")
    text_neg <- paste0(n_neg, " stocks\n \n \n")
  }
  text(x=-1.65, y=12, labels=text_pos, pos=4, adj=1, cex=0.9, col="blue")
  text(x=1.65, y=nrow(vals)-18, labels=text_neg, pos=2, adj=0, cex=0.9, col="red")
  
}

# Add x-axis label
# xlabel <- expression("SST influence (θ"["i"]*")")
xlabel <- "Temperature influence"
mtext(xlabel, outer=T, side=1, adj=0.32, line=1.2, cex=0.8)

# SST influence examples
#######################################

# Read data
stocks <- read.csv(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.csv", sep="/"), as.is=T)

# Example stocks
ex_stocks <- c("BSBASSMATLC", "CODIS", "HERRNWATLC")
ex_text <- c("Warming ↑\nproductivity", "Warming ↓\nproductivity", "Warming has\nno effect")

# Parameters
p <- 0.2
ymin <- c(0, -4, -400)
ymax <- c(7, 20, 1200)
ybin <- c(1, 4, 400)
xmin <- c(0, 0, 0)
xmax <- c(16, 30, 2000)
xbin <- c(4, 5, 500)


# Plotting par
par(mar=c(1,3.5,0.5,0.5), mgp=c(2,0.7,0))

# Loop through example stocks
for(i in 1:length(ex_stocks)){

  # Stock info
  stock <- filter(stocks, stockid==ex_stocks[i])
  r <- stock$r
  k_orig <- stock$k
  k <- k_orig * stock$tb_max / 1000
  betaT <- stock$betaT

  # Subset data
  sdata <- data %>%
    filter(stockid==ex_stocks[i]) %>%
    mutate(sp=sp/1000,
           tb=tb/1000)

  # Add color bins to data
  sdata$sst_bin <- cut(sdata$cobe_sst_c_sd, breaks=seq(-2,2,0.1))
  pcolors <- tcolor(rev(colorpal(brewer.pal(11,"RdBu"),nlevels(sdata$sst_bin))), 0.8)
  sdata$sst_col <- pcolors[sdata$sst_bin]

  # Plot data
  plot(sp ~ tb, sdata, bty="n", las=1, cex=1.5, pch=21, bg=sdata$sst_col, #bg="grey80",
       xaxt="n", yaxt="n", xlim=c(xmin[i], xmax[i]), ylim=c(ymin[i], ymax[i]), cex.lab=1.2,
       xlab="", ylab="", xpd=NA)
  axis(1, at=seq(xmin[i], xmax[i], xbin[i]), las=1, cex.axis=0.9)
  axis(2, at=seq(ymin[i], ymax[i], ybin[i]), las=1, cex.axis=0.9)

  # Add SP curve
  curve(r/p*x*(1-(x/k)^p),
        from=xmin[i], to=xmax[i], n=101, add=T, lty=1, lwd=1.8, col="black")

  # Add SP-SST curves
  ssts <- c(-1, -0.5, 0.5, 1)
  colors <- rev(brewer.pal(10, "RdBu"))[c(1,3,8,10)]
  for(j in 1:length(ssts)){
    curve(r/p*x*(1-(x/k)^p)*exp(betaT*ssts[j]),
          from=xmin[i], to=xmax[i], n=101, add=T, lty=1, lwd=1.3, col=colors[j])
  }
  
  # Add degree labels
  if(i==1){
    text(x=4.5, y=4, pos=2, labels="+1°C", cex=0.9, col=colors[4])
    text(x=4, y=1.3, pos=4, labels="-1°C", cex=0.9, col=colors[1])
  }
  if(i==2){
    text(x=10, y=11.5, pos=2, labels="-1°C", cex=0.9, col=colors[1])
    text(x=20, y=1, pos=4, labels="+1°C", cex=0.9, col=colors[4])
  }
  
  # Add warming effect label
  x <- 0-xmax[i]*0.05
  y <- ymax[i]-(ymax[i]-ymin[i])*0.12
  text(x=x, y=y, labels=ex_text[i], pos=4,  cex=0.9, col="grey50")
  
  # Add letter
  text(x=xmax, y=ymax, labels=LETTERS[i+2], pos=1, font=2, cex=1.1, offset=0.2)
  

}

# Axis labels
mtext("Biomass (1000s mt)", outer=T, side=1, adj=0.97, line=1.2, cex=0.8)
mtext("Production (1000s mt)", outer=T, side=2, adj=0.52, line=-35.3, cex=0.8)

# Off
dev.off()
graphics.off()


