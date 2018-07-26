
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"
faodir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/fao_stat_areas"


# Build data
################################################################################

# Groups
groups <- c("order", "family", "lme", "fao_area")
group_sds <- rep(NA, length(groups))
  
# Loop through groups
for(i in 1:length(groups)){
  
  # Read group model results
  group1 <- groups[i]
  model_name <- paste0("ramldb_v3.8_spsst_pella_cobe_", group1, ".Rdata")
  load(paste(datadir, model_name, sep="/"))
  rm(hess, input.data, model, output, sd, nstocks, params, problem.stocks)
  
  # Sample size by group
  colnames(data)[which(colnames(data)==group1)] <- "group"
  ngroup <- data %>% 
    group_by(group) %>% 
    summarize(n=n_distinct(stockid))
  
  # Subset group results
  vals <- results.df %>% 
    # Reduce to group means
    filter(param=="mu_group") %>% 
    rename(group=stockid,
           est=estimate,
           est_se=stderror) %>% 
    select(-param) %>% 
    # Add sample size
    left_join(ngroup, by="group") %>% 
    # Format data
    mutate(type=group1, 
           label=paste0(group, " (", n, ")"),
           est_lo=est-est_se*1.96, 
           est_hi=est+est_se*1.96, 
           est_inf="none",
           est_inf=ifelse(est_hi < 0, "negative", est_inf), 
           est_inf=ifelse(est_lo > 0, "positive", est_inf),
           pcolor=revalue(est_inf, c("positive"="blue",
                                     "negative"="red",
                                     "none"="black")),
           lcolor=revalue(est_inf, c("positive"="blue",
                                     "negative"="red",
                                     "none"="grey60")),
           tcolor=revalue(est_inf, c("positive"="blue",
                                     "negative"="red",
                                     "none"="grey30"))) %>% 
    # Rearrange and sort
    select(type, group, label, n, everything()) %>% 
    arrange(desc(est))
  
  # Group SD
  group_sds[i] <- results.df$estimate[results.df$param=="sd_group"]
  
  # Merge together
  if(i==1){results <- vals}else{results <- rbind(results, vals)}
  
}

# Format data
################################################################################

# Format data
results1 <- results %>% 
  mutate(group=revalue(group, c("California Current"="CA Current",
                                "East Bering Sea"="E Bering Sea",
                                "East China Sea"="E China Sea",
                                "Iceland Shelf and Sea"="Iceland Shelf & Sea",
                                "Labrador - Newfoundland"="Labrador-Newfoundland",
                                "North Atlantic Ocean"="N Atlantic Ocean",
                                "North Pacific Ocean"="N Pacific Ocean",
                                "Northeast U.S. Continental Shelf"="NE US Shelf",
                                "South Atlantic Ocean"="S Atlantic Ocean",
                                "South Pacific Ocean"="S Pacific Ocean",
                                "South West Australian Shelf"="SW Australian Shelf",
                                "Southeast Australian Shelf"="SE Australian Shelf",
                                "Southeast U.S. Continental Shelf"="SE US Shelf",
                                "West Bering Sea"="W Bering Sea",
                                "Atlantic, Northeast"="Atlantic, NE",
                                "Atlantic, Southeast"="Atlantic, SE",
                                "Atlantic, Southwest"="Atlantic, SW",
                                "Atlantic, Northwest"="Atlantic, NW",
                                "Atlantic, Eastern Central"="Atlantic, E Central",
                                "Atlantic, Western Central"="Atlantic, W Central",
                                "Indian Ocean, Western"="Indian, W",
                                "Indian Ocean, Eastern"="Indian, E",
                                "Mediterranean and Black Sea"="Mediterranean & Black Sea", 
                                "Pacific, Northwest"="Pacific, NW",
                                "Pacific, Eastern Central"="Pacific, E Central",
                                "Pacific, Southwest"="Pacific, SW",
                                "Pacific, Southeast"="Pacific, SE",
                                "Pacific, Northeast"="Pacific, NE")),
         label=paste0(group, " (", n, ")"))


# Plot data
################################################################################

# Groups
groups1 <- c("lme", "family", "fao_area", "order")
groups1_labels <- c("Marine\necoregion", "Family", "FAO major\nfishing area", "Order")
cexs <- c(0.7, 0.65, 0.7, 0.7)

# Setup figure
figname <- "Fig2_thetas_taxa_region_models.png"
png(paste(plotdir, figname, sep="/"), width=6, height=6, units="in", res=600)
layout(matrix(c(1,2,3,4), ncol=2, byrow=T), heights=c(0.7,0.3))
par(mar=c(0.7,4,0.1,0.5), oma=c(2.3,0,1,0), mgp=c(2, 0.5, 0))

# Loop through groups
for(i in 1:length(groups1)){
  
  # Subset data
  vals1 <- subset(results1, type==groups1[i])
  
  # Setup empty plot
  plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=0.75,
       xaxt="n", xlim=c(-1.5, 1.5), ylim=c(1,nrow(vals1)),
       xlab="", ylab="", cex.main=0.8)
  if(i<3){
    axis(1, at=seq(-1.5, 1.5, 0.5), labels=F)
  }else{    
    axis(1, at=seq(-1.5, 1.5, 0.5), labels=T, cex.axis=0.85)
  }
  
  # Add theta estimates (points) and errors (lines)
  sapply(1:nrow(vals1), function(x) lines(x=c(vals1$est_lo[x], vals1$est_hi[x]), y=c(x,x), col=vals1$lcolor[x], lwd=0.6))
  points(vals1$est, 1:nrow(vals1), pch=16, cex=0.6, col=vals1$pcolor)
  
  # Add theta=0 line
  lines(x=c(0,0), y=c(1, nrow(vals1)), lty=3, col="black", lwd=0.8)
  
  # Add labels
  text(x=-2.5, y=1:nrow(vals1), pos=4, labels=vals1$label, col=vals1$tcolor, cex=cexs[i], xpd=NA)
  
  # Add model
  text(x=1.5, y=nrow(vals1)-nrow(vals1)*0, labels=groups1_labels[i], adj=c(1,1), font=2, cex=0.85)
  
}

# Add x-axis label
mtext(expression("SST influence (θ"["i"*")"]), outer=T, side=1, adj=0.55, line=1.2, cex=0.9)

# Add titles
mtext("Geography", outer=T, side=3, adj=0.26, line=-0.2, cex=0.9, font=2)
mtext("Taxonomy", outer=T, side=3, adj=0.84, line=-0.2, cex=0.9, font=2)


# Off
dev.off()
graphics.off()

