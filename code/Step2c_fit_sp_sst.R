
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(TMB)
library(plyr)
library(dplyr)
library(devtools)
library(reshape2)
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Directories
tmbdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/code/tmb_code"
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/input"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/code"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb1/output"
setwd(tmbdir)

# Read data
load(paste(datadir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))

# Source helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# Which SST and model?
# SST choices: cobe, er, had
# Model choices: final, n1, n2, n3
sst <- "had"
model <- "final"


# Format data
################################################################################

# Inspect data completeness
apply(data, 2, function(x) sum(is.na(x)))

# Remove problem stocks
stocks.nls.no <- spfits$stockid[is.na(spfits$r)]
stocks.w.big.k <- c("AMPL3LNO", "ARGHAKESARG", "BRMSOJ", "MORWONGSE", "RBRMECS", 
                    "REXSOLEGA", "RSNAPGM", "WPOLLNSO")
problem.stocks <- c(stocks.nls.no, stocks.w.big.k)
data <- subset(data, !(stockid%in%problem.stocks))


# Select data and model
################################################################################

# Parameters
stocks <- unique(data$stockid)
nstocks <- length(stocks)

# Specify SST data to use and output name
if(model=="final"){
  sst.data <- as.matrix(data[,paste0(sst, "_sst_c_sd")])
  outputfile <- paste0("ramldb_v3.8_spsst_", sst, ".Rdata")
}else{
  sst.data <- as.matrix(data[,paste0(sst, "_sst_c_", model, "_sd")])
  outputfile <- paste0("ramldb_v3.8_spsst_", sst, "_", model, ".Rdata")
}


# Fit production model
################################################################################

# Compile TMB code
# Only run once to compile code
if(FALSE){
  dyn.unload(paste(tmbdir, dynlib("spmodel_sst_link"), sep="/"))
  file.remove(paste(tmbdir, c("spmodel_sst_link.o", "spmodel_sst_link.dll"), sep="/"))
  compile("spmodel_sst_link.cpp")
}

# Load TMB code
dyn.load(dynlib("spmodel_sst_link"))

# Input data and parameter starting values
b0_starts <- log(sapply(stocks, function(x) max(data$tb_sd[data$stockid==x])) * 1.5)
params <- list(ln_B0=b0_starts,
               ln_r=rep(log(0.4), nstocks),
               BetaT=rep(0.0, nstocks),
               ln_sigmaP=rep(-2.5, nstocks), # 0 before, -2.5 based on model fits
               mu_T=0.0,
               ln_sd_T=-1.25) # -3 before, -1.25 based on model fits
input.data <- list(Nstocks=nstocks,
                   Nobs=nrow(data),
                   StockID=as.factor(data$stockid),
                   B_t=data$tb_sd,
                   P_t=data$sp_sd,
                   Temp_t=sst.data)

# Initialization
model <- MakeADFun(data=input.data, parameters=params, random="BetaT", DLL="spmodel_sst_link")
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


# SAVE RESULTS
################################################################################

# Export model objects
save(data, stocks, problem.stocks, nstocks,
     input.data, params,
     model, output, sd, hess,
     results.df, results.wide,
     file=paste(outputdir, outputfile, sep="/"))
