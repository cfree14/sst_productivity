
# Clear workspace
rm(list = ls())

# SETUP
################################################################################

# Packages
library(TMB)
library(plyr)
library(dplyr)
library(devtools)
library(reshape2)
library(freeR)
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Directories - original submission
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/code"
outputdir_orig <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"

# Directories - revisions
tmbdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/code/tmb_code"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisionsfigures"
setwd(tmbdir)

# Read data
load(paste(datadir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))
data_orig <- data; rm(data)

# Read starting values
starts <- read.csv(file.path(datadir, "sp_ar1_starting_values.csv"), as.is=T)

# Read FAO data
faodir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/fao_stat_areas"
fao <- read.csv(paste(faodir, "fao_stat_area_key.csv", sep="/"), as.is=T)
fao <- rename(fao, fao_area=name)

# Read results from other
key_orig <- read.csv(paste(outputdir_orig, "ramldb_v3.8_spsst_cobe.csv", sep="/"), as.is=T)
key <- key_orig %>% 
  rename(lme=lme_name, fao_code=fao_area) %>% 
  left_join(fao, by=c("fao_code"="code")) %>% 
  arrange(stockid)

# Add groups (order, family, LME, FAO area) to data
# Make sure the dataframe is sorted by STOCKID and YEAR
data <- data_orig %>%
  left_join(select(key, stockid, order, family, fao_area, lme, method, method1, method2), by="stockid") %>% 
  select(assessid, stockid, order, family, fao_area, lme, method, method1, method2, everything()) %>% 
  arrange(stockid, year)

# Source helper functions
source(paste(codedir, "helper_functions.R", sep="/"))


# FORMAT DATA
################################################################################

# Remove problem stocks
stocks.nls.no <- spfits$stockid[is.na(spfits$r)] 
stocks.w.big.k <- c("AMPL3LNO", "ARGHAKESARG", "BRMSOJ", "MORWONGSE", "RBRMECS", 
                    "REXSOLEGA", "RSNAPGM", "WPOLLNSO")
problem.stocks <- c(stocks.nls.no, stocks.w.big.k)
data <- subset(data, !(stockid%in%problem.stocks))

# Completeness - MUST HAVE FAMILY, ORDER, LME_NAME, 
complete(data)


# FIT PRODUCTION MODEL
################################################################################

# Which group? And shape parameter?
# Options: "order", "family", "lme", "fao_area", "method", "method1"
group <- "lme"
group_data <- key[,group]
p <- 0.55

# Which SST dataset?
# Options: "obs", "n1", "n2", "n3"
sst <- "n3"
if(sst=="obs"){
  sst.data <- data$cobe_sst_c_sd
  outputfile <- paste0("ramldb_v3.8_sp_ar1_pella0.55_cobe_", group, ".Rdata")
}else{
  sst.data <- as.matrix(data[,paste0("cobe_sst_c_", sst, "_sd")])
  outputfile <- paste0("ramldb_v3.8_sp_ar1_pella0.55_cobe_", group, "_", sst,  ".Rdata")
}

# Parameters
stocks <- unique(data$stockid)
nstocks <- length(stocks)
groups <- unique(group_data)
ngroups <- length(groups)

# Compile TMB code
# Only run once to compile code
if(FALSE){
  dyn.unload(paste(tmbdir, dynlib("spmodel_ar1_pella_sst_group"), sep="/"))
  file.remove(paste(tmbdir, c("spmodel_ar1_pella_sst_group.o", "spmodel_ar1_pella_sst_group.dll"), sep="/"))
  compile("spmodel_ar1_pella_sst_group.cpp")
}

# Load TMB code
dyn.load(dynlib("spmodel_ar1_pella_sst_group"))

# Input data and parameter starting values
# params <- list(ln_B0=log(rep(1.5, nstocks)),
#                ln_r=log(rep(0.3, nstocks)),
#                BetaT=rep(0.0, nstocks),
#                ln_sigmaP=log(rep(0.15, nstocks)), 
#                rho=starts$rho,
#                # rho=rep(0, nstocks),
#                mu_T=0.0,
#                ln_sd_T=-1.25,
#                mu_group=rep(0, ngroups),
#                ln_sd_group=-1.25)
params <- list(ln_B0=log(pmin(starts$k, 1.5)),
               ln_r=log(pmin(starts$r, 0.8)),
               BetaT=rep(0.0, nstocks),
               ln_sigmaP=log(starts$sigmaP),
               rho=starts$rho,
               mu_T=0.0,
               ln_sd_T=-1.25,
               mu_group=rep(0, ngroups),
               ln_sd_group=-1.25)
input.data <- list(Nstocks=nstocks,
                   Nobs=nrow(data),
                   Ngroups=ngroups,
                   p=p,
                   StockID=as.factor(data$stockid),
                   Group=as.factor(group_data),
                   B_t=data$tb_sd,
                   P_t=data$sp_sd,
                   Temp_t=sst.data)

# Initialization
model <- MakeADFun(data=input.data, parameters=params, random=c("BetaT", "mu_group"), DLL="spmodel_ar1_pella_sst_group")
model$control <- list(trace=1, parscale=rep(1,13), REPORT=1, reltol=1e-12, maxit=100)
model$hessian <- F
newtonOption(model, smartsearch=TRUE)

# Run model
# output <- TMBhelper::Optimize(obj=model, lower=-Inf, upper=Inf, loopnum=3, newtonsteps=3, bias.correct=FALSE, getsd=FALSE)
output <- nlminb(start=model$par, objective=model$fn, gradient=model$gr, lower=-Inf, upper=Inf, control=list(eval.max=1e4, iter.max=1e4))


# MODEL DIAGNOSTICS
################################################################################

# Use hessian to diagnose fixed effects that might cause a problem
hess <- optimHess(par=output$par, fn=model$fn, gr=model$gr)
problem.vals <- which(eigen(hess)$values<0)
if(length(problem.vals)>0 ){
  display <- eigen(hess)$vectors[,problem.vals]
  names(display) = (output$diagnostics$Param)
  cbind(1:length(output$par), output$par, display)
}

# Calculate SD
sd <- try(sdreport(model, hessian.fixed=hess))

# AIC of model
TMBhelper::TMBAIC(output)


# FORMAT RESULTS
################################################################################

# Build parameter estimate dataframe
############################################

# Format parameter estimates
results.mat <- summary.sdreport(sd)
results.df <- data.frame(param=rownames(results.mat),
                         estimate=results.mat[,1],
                         stderror=results.mat[,2])

# Add stockid to parameter estimates
n_params <- table(results.df$param)
global_params <- names(n_params[n_params==1])
group_params <- names(n_params[n_params==ngroups])
stock_params <- names(n_params[n_params==length(stocks)])
results.df$stockid[results.df$param%in%global_params] <- "global"
results.df$stockid[results.df$param%in%group_params] <- sort(groups)
for(i in 1:length(stock_params)){
  results.df$stockid[results.df$param==stock_params[i]] <- stocks
}
results.df <- subset(results.df, select=c(stockid, param, estimate, stderror))

# Reshape dataframe: long-to-wide
# results.wide <- dcast(results.df, stockid ~ param, value.var="estimate", subset=.(stockid%in%stocks))
results.wide <- format_results(results.df)

# Inspect for bad fits
############################################

# Histogram parameter values
params <- unique(results.df$param)
par(mfrow=c(2,5))
for(i in 1:length(params)){
  vals <- results.df$estimate[results.df$param==params[i]]
  hist(vals, main="", xlab=params[i], ylab="")
}

# Plot results
plot_ests(results.df)
plot_thetas(results.df)

# Plot group mus
group_mus <- results.df %>% 
  filter(param=="mu_group") %>% 
  arrange(desc(estimate))
hist(group_mus$estimate)
abline(v=results.df$estimate[results.df$param=="mu_T"])

# Plot spline
plot(x=group_mus$estimate, y=1:nrow(group_mus), bty="n", yaxt="n", ylab="")
axis(2, at=1:nrow(group_mus), labels=group_mus$stockid, las=2, cex.axis=0.6)
abline(v=0, lty=3)

# Inspect residuals over time
################################################################################

# # Create residuals dataframe
# resids <- model$report()
# resids_df <- data %>% 
#   ungroup() %>% 
#   select(stockid, year) %>% 
#   mutate(pt_pred=resids[[1]])
# 
# # Setup plot
# pdf(file.path(plotdir, "sp_pella_group_model_resid_check.pdf"), width=8.5, height=11)
# par(mfrow=c(6,4), oma=c(2,2,2,2))
# 
# # Loop through and plot
# for(i in 1:nstocks){
#   stock <- stocks[i]
#   sdata <- filter(resids_df, stockid==stock)
#   plot(pt_pred ~ year, sdata, type="l", bty="n", xlab="", ylab="SP residuals", las=1, main=stock)
#   abline(h=0, lty=3)
# }
# 
# dev.off()


# SAVE RESULTS
################################################################################

# Export model objects
save(data, stocks, problem.stocks, nstocks, 
     input.data, params, 
     model, output, sd, hess,
     results.df, results.wide,
     file=paste(outputdir, outputfile, sep="/"))
