

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/figures"

# Read base model
load(paste(datadir, "ramldb_v3.8_spsst_pella_0.20_cobe.Rdata", sep="/"))
rm(data, hess, input.data, model, output, results.df, sd, nstocks, params, problem.stocks, stocks)
base <- results.wide

# Read order model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_order.Rdata", sep="/"))
rm(data, hess, input.data, model, output, results.df, sd, nstocks, params, problem.stocks, stocks)
order <- results.wide

# Read family model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_family.Rdata", sep="/"))
rm(data, hess, input.data, model, output, results.df, sd, nstocks, params, problem.stocks, stocks)
family <- results.wide

# Read LME model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_lme.Rdata", sep="/"))
rm(data, hess, input.data, model, output, results.df, sd, nstocks, params, problem.stocks, stocks)
lme <- results.wide

# Read FAO area model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_fao_area.Rdata", sep="/"))
rm(data, hess, input.data, model, output, results.df, sd, nstocks, params, problem.stocks, stocks)
fao_area <- results.wide

# Read method #1 model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_method.Rdata", sep="/"))
rm(data, hess, input.data, model, output, results.df, sd, nstocks, params, problem.stocks, stocks)
method1 <- results.wide

# Read method #2 model
load(paste(datadir, "ramldb_v3.8_spsst_pella_cobe_method1.Rdata", sep="/"))
rm(data, hess, input.data, model, output, results.df, sd, nstocks, params, problem.stocks, stocks)
method2 <- results.wide


# Build data
################################################################################

# Build data
data <- base %>% 
  select(stockid, betaT) %>% 
  rename(betaT_base=betaT) %>% 
  # Add ORDER results
  left_join(select(order, stockid, betaT), by="stockid") %>% 
  rename(betaT_order=betaT) %>% 
  # Add FAMILY results
  left_join(select(family, stockid, betaT), by="stockid") %>% 
  rename(betaT_family=betaT) %>% 
  # Add FAO AREA results
  left_join(select(fao_area, stockid, betaT), by="stockid") %>% 
  rename(betaT_fao=betaT) %>% 
  # Add LME results
  left_join(select(lme, stockid, betaT), by="stockid") %>% 
  rename(betaT_lme=betaT) %>% 
  # Add METHOD #1 results
  left_join(select(method1, stockid, betaT), by="stockid") %>% 
  rename(betaT_sa1=betaT) %>% 
  # Add METHOD #2 results
  left_join(select(method2, stockid, betaT), by="stockid") %>% 
  rename(betaT_sa2=betaT) 



# Plot data
################################################################################

# Setup figure
figname <- "SFig7_group_model_correlation.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6, units="in", res=600)
par(mfcol=c(7,7), mar=c(1,2,1,0.5), mgp=c(2.2, 0.8, 0), oma=c(4,4,0,0))

# Loop through models and plot thetas
cols <- c("base", "order", "family", "fao", "lme", "sa1", "sa2")
model_names <- c("Base", "Order", "Family", "FAO area", "Ecoregion", "Method #1", "Method #2")

# Loop through results
for(i in 1:length(cols)){
  for(j in 1:length(cols)){
    # If i==j, don't plot
    if(i==j | i>j){
      plot.new()
    }else{
      # Get data
      x_vals <- data[,paste0("betaT_", cols[i])]
      y_vals <- data[,paste0("betaT_", cols[j])]
      # Plot data
      xlab <- ifelse(j==7, model_names[i], "")
      ylab <- ifelse(i==1, model_names[j], "")
      plot(y_vals ~ x_vals, bty='n', las=1, col="grey50",
           xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), cex.axis=0.8, cex.lab=1,
           xlab=xlab,  ylab=ylab, xpd=NA)
      lines(x=c(-1.5,1.5), y=c(-1.5, 1.5))
      lines(x=c(-1.5,1.5), y=c(0, 0), lty=3, col="grey40")
      lines(y=c(-1.5,1.5), x=c(0, 0), lty=3, col="grey40")
      # Add r2 text
      r2_format <- roundf(r2(lm(y_vals ~ x_vals)), 2)
      text(x=1.7, y=-1.3, pos=2, labels=bquote("r"^"2"*"="*.(r2_format)), cex=0.7)
      
    }
  }
}

# Add more labels
mtext(expression("SST influence (θ"["i"]*")"), outer=1, side=1, adj=0.42, line=2.8, cex=0.9)
mtext(expression("SST influence (θ"["i"]*")"), outer=1, side=2, adj=0.42, line=2, cex=0.9)

# Off
dev.off()
graphics.off()

