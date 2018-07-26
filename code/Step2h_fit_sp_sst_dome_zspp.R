
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
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Directories
tmbdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/code/tmb_code"
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/input"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/code"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/output"
setwd(tmbdir)

# Read data
load(paste(datadir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))
data_orig <- data; data <- final; rm(final)

# Source helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# FORMAT DATA
################################################################################

# Inspect data completeness
apply(data, 2, function(x) sum(is.na(x)))

# Remove problem stocks
problem.stocks <- spfits$assessid[is.na(spfits$r)] # 26 stocks
data <- subset(data, !(assessid%in%problem.stocks))


# FIT PRODUCTION MODEL
################################################################################

# Parameters
spp <- unique(data$species)
nspp <- length(spp)
stocks <- unique(data$stockid)
nstocks <- length(stocks)

# Which SST?
# yr_t = present year
# shuf = SST times series shuffled among stocks
# ar = SST time series simulated with same mean/autocorrelation but no trend
# lin_ar = SST time series simulated with same mean/autocorrelation/trend
# lin = SST time series simulated with same trend but no mean/variance
sst <- "yr_t"
if(sst=="yr_t"){
  sst.data <- data$sst_t_sd
  outputfile <- "ramldb_v3.8_spmodel_temp_dpdt_sst_yr_t.Rdata"
}
if(sst=="shuf"){
  sst.data <- data$sst_shuf_sd
  outputfile <- "ramldb_v3.8_spmodel_temp_dpdt_sst_shuf.Rdata"
}
if(sst=="ar"){
  sst.data <- data$sst_ar_sd
  outputfile <- "ramldb_v3.8_spmodel_temp_dpdt_sst_ar.Rdata"
}
if(sst=="lin"){
  sst.data <- data$sst_lin_sd
  outputfile <- "ramldb_v3.8_spmodel_temp_dpdt_sst_lin.Rdata"
}
if(sst=="lin_ar"){
  sst.data <- data$sst_lin_ar_sd
  outputfile <- "ramldb_v3.8_spmodel_temp_dpdt_sst_lin_ar.Rdata"
}


# Compile TMB code
# Only run once to compile code
if(FALSE){
  dyn.unload(paste(tmbdir, dynlib("spmodel_sst_dome_zspp"), sep="/"))
  file.remove(paste(tmbdir, c("spmodel_sst_dome_zspp.o", "spmodel_sst_dome_zspp.dll"), sep="/"))
  compile("spmodel_sst_dome_zspp.cpp")
}

# Load TMB code
dyn.load(dynlib("spmodel_sst_dome_zspp"))

# Input data and parameter starting values
b0_starts <- log(sapply(stocks, function(x) max(data$tb_sd[data$stockid==x])) * 1.5)
z_starts <- sapply(spp, function(x) mean(data$sst_t_sd[data$species==x])); hist(z_starts)
params <- list(ln_B0=b0_starts,
               ln_r=rep(log(0.4), nstocks),
               ln_BetaT=rep(log(0.01), nstocks),
               Z=z_starts,
               ln_sigmaP=rep(-2.5, nstocks), # 0 before, -2.5 based on model fits
               mu_T=0.0,
               ln_sd_T=-1.25) # -3 before, -1.25 based on model fits
input.data <- list(Nstocks=nstocks,
                   Nobs=nrow(data),
                   StockID=as.factor(data$stockid),
                   SpeciesID=as.factor(data$species),
                   B_t=data$tb_sd,
                   P_t=data$sp_sd,
                   Temp_t=sst.data)

# Initialization
model <- MakeADFun(data=input.data, parameters=params, random="ln_BetaT", DLL="spmodel_sst_dome")
# model$control <- list(trace=1, parscale=rep(1,13), REPORT=1, reltol=1e-12, maxit=100)
# model$hessian <- F
# newtonOption(model, smartsearch=TRUE)

# Run model
output <- TMBhelper::Optimize(obj=model, lower=-Inf, upper=Inf, loopnum=3, newtonsteps=3, bias.correct=FALSE, getsd=FALSE)


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

# Extract MLE of fixed effects and Empirical Bayes prediction of random effects
ParHat <- model$env$parList()

# Compare sample statistics vs. population statistics for random effects
mean(ParHat$BetaT); output$par['mu_T']
sd(ParHat$BetaT); exp(output$par['ln_sd_T'])

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
stock_params <- names(n_params[n_params==length(stocks)])
results.df$stockid[results.df$param%in%global_params] <- "global"
for(i in 1:length(stock_params)){
  results.df$stockid[results.df$param==stock_params[i]] <- stocks
}
results.df <- subset(results.df, select=c(stockid, param, estimate, stderror))

# Reshape dataframe: long-to-wide
results.wide <- dcast(results.df, stockid ~ param, value.var="estimate", subset=.(stockid!="global"))

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


# SAVE RESULTS
################################################################################

# Export model objects
save(data, stocks, problem.stocks, nstocks, 
     input.data, params, 
     model, output, sd, hess,
     results.df, results.wide,
     file=paste(outputdir, outputfile, sep="/"))
