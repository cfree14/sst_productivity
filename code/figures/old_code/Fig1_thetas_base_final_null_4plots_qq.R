
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)
library(RColorBrewer)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/figures"


# Plot data
################################################################################

# Setup figure
figname <- "Fig1_thetas_final_null1.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=4, units="in", res=600)
layout(matrix(1:8, ncol=4), heights=c(0.7, 0.3))
par(oma=c(2,2,0,0))

# Loop through models
models <- c("", "_lme", "_family", "_n1")
model_names <- c("Base model", "LME model", "Family model", "Null model")
for(i in 1:length(models)){

  # Load data
  model_obj <- paste0("ramldb_v3.8_spsst_cobe", models[i], ".Rdata")
  load(paste(datadir, model_obj, sep="/"))
  rm(hess, data, input.data, model, nstocks, stocks, output, params, problem.stocks, sd)
  
  # Format data
  ##########################################
  
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
  
  # A. Spline plot
  ##########################################
  
  # Par
  par(mar=c(3.3,0.5,1,0.5), mgp=c(2.8,0.5,0))
  
  # Setup empty plot
  plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=0.75,
       xlim=c(-1.5, 1.5), ylim=c(1,nrow(vals)),
       xlab="", ylab="", main=model_names[i], cex.main=0.8)
  
  # Add theta mean
  mu_est <- mu$est[mu$param=="mu_T"]
  mu_min <- mu$est_lo[mu$param=="mu_T"]
  mu_max <- mu$est_hi[mu$param=="mu_T"]
  rect(xleft=mu_min, xright=mu_max, ytop=nrow(vals), ybottom=1, border=NA, col="grey80")
  text(x=0, y=nrow(vals), labels=paste0("μ = ", format(round(mu_est,2), nsmall=2)), pos=4, cex=0.7, col="grey30")

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
  text(x=-1.65, y=8, labels=text_pos, pos=4, adj=1, cex=0.7, col="blue")
  text(x=1.65, y=nrow(vals)-8, labels=text_neg, pos=2, adj=0, cex=0.7, col="red")
  
  # B. QQ plot
  ##########################################
  
  # Par
  par(mar=c(1,2,0.5,1), mgp=c(2.8,0.8,0))
  
  # Build data
  qq <- data.frame(qqnorm(vals$est, plot.it=F), col=vals$lcolor, stringsAsFactors=F)
  
  # Plot data
  if(i%in%1:3){
    ymin <- -1.5
    ymax <- 1.5
    ybin <- 0.5
  }else{
    ymin <- -0.25
    ymax <- 0.5
    ybin <- 0.25
  }
  qqnorm(vals$est, col=vals$pcolor, cex=1, las=1, type="n", yaxt="n", cex.axis=0.75,
         xlab="", ylab="", ylim=c(ymin, ymax), main="",  bty="n", xpd=NA)
  axis(2, at=seq(ymin, ymax, ybin), las=1, cex.axis=0.75)
  for(i in c("grey60", "red", "blue")){
    qqs <- subset(qq, col==i)
    points(x=qqs$x, y=qqs$y, col=qqs$col, xpd=NA, cex=1)
  }
  qqline(vals$est)
  
  
}

# Add axis labels
xlabel1 <- expression("SST influence (θ"["i"]*")")
xlabel2 <- "Theoretical quantiles"
ylabel2 <- "Observed quantiles"
mtext(xlabel1, outer=T, side=1, adj=0.52, line=-9.5, cex=0.75)
mtext(xlabel2, outer=T, side=1, adj=0.52, line=0.95, cex=0.75)
mtext(ylabel2, outer=T, side=2, adj=0.01, line=0.3, cex=0.75)

# Off
dev.off()
graphics.off()


