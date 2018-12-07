
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/figures"
faodir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/fao_stat_areas"

# Read stock meta-data
mdata <- read.csv(file.path(datadir, "ramldb_v3.8_sp_ar1_pella0.55_cobe_lme.csv"), as.is=T)

# Build data
################################################################################

# Groups
groups <- c("method1", "method")
group_sds <- rep(NA, length(groups))
  
# Loop through groups
for(i in 1:length(groups)){
  
  # Read group model results
  group1 <- groups[i]
  model_name <- paste0("ramldb_v3.8_sp_ar1_pella0.55_cobe_", group1, ".Rdata")
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
results$group[results$group=="VPA" & results$type=="method"] <- "Virtual population analysis"
results1 <- results %>%
  mutate(group=revalue(group, c("Statistical catch at age model"="Statistical catch-at-age model",
                                "Statistical catch at length model"="Statistical catch-at-length model",
                                "Integrated Analysis"="Integrated analysis",
                                "unknown"="Unknown")),
         label=paste0(group, " (", n, ")"))



# Investigate significant influences
################################################################################

# Sig assessment methods
sig_sas <- filter(results1, est_inf%in%c("positive", "negative"))$group
mdata1 <- filter(mdata, method1 %in% sig_sas)

# Plot data
################################################################################

# Groups
groups1 <- c("method", "method1")
groups1_labels <- c("Generic\nassessment\nmethod", "Specific\nassessment\nmethod")

# Setup figure
figname <- "SFig4_thetas_assessment_models.png"
png(paste(plotdir, figname, sep="/"), width=6, height=6, units="in", res=600)
layout(matrix(c(1,2,3,2), ncol=2, byrow=T), heights=c(0.4,0.6))
par(mar=c(3,4,0.5,0.5), mgp=c(2, 0.5, 0))

# Loop through groups
for(i in 1:length(groups1)){
  
  # Subset data
  vals1 <- subset(results1, type==groups1[i])
  
  # Setup empty plot
  plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=0.75,
       xlim=c(-1.5, 1.5), ylim=c(1,nrow(vals1)),
       xlab=expression("SST influence (θ"["i"]*")"), ylab="", cex.main=0.8)
  
  # Add theta estimates (points) and errors (lines)
  sapply(1:nrow(vals1), function(x) lines(x=c(vals1$est_lo[x], vals1$est_hi[x]), y=c(x,x), col=vals1$lcolor[x], lwd=0.6))
  points(vals1$est, 1:nrow(vals1), pch=16, cex=0.6, col=vals1$pcolor)
  
  # Add theta=0 line
  lines(x=c(0,0), y=c(1, nrow(vals1)), lty=3, col="black", lwd=0.8)
  
  # Add labels
  text(x=-2.5, y=1:nrow(vals1), pos=4, labels=vals1$label, col=vals1$tcolor, cex=0.7, xpd=NA)
  
  # Add model
  text(x=1.5, y=nrow(vals1)-nrow(vals1)*0, labels=groups1_labels[i], adj=c(1,1), font=2, cex=0.85)
  
}


# Off
dev.off()
graphics.off()


# Plot data
################################################################################

# Setup figure
figname <- "SFig4_thetas_assessment_models_general_only.png"
png(paste(plotdir, figname, sep="/"), width=4, height=3, units="in", res=600)
par(mar=c(3,5,0.5,0.5), mgp=c(2, 0.5, 0))

# Subset data
vals1 <- subset(results1, type=="method")

# Setup empty plot
plot(1:10, 1:10, type="n", bty="n", yaxt="n", cex.axis=0.75,
     xlim=c(-1, 1), ylim=c(1,nrow(vals1)),
     xlab=expression("SST influence (θ"["i"]*")"), ylab="", cex.main=0.8)

# Add theta estimates (points) and errors (lines)
sapply(1:nrow(vals1), function(x) lines(x=c(vals1$est_lo[x], vals1$est_hi[x]), y=c(x,x), col=vals1$lcolor[x], lwd=0.6))
points(vals1$est, 1:nrow(vals1), pch=16, cex=0.6, col=vals1$pcolor)

# Add theta=0 line
lines(x=c(0,0), y=c(1, nrow(vals1)), lty=3, col="black", lwd=0.8)

# Add labels
text(x=-1.8, y=1:nrow(vals1), pos=4, labels=vals1$label, col=vals1$tcolor, cex=0.7, xpd=NA)


# Off
dev.off()
graphics.off()




