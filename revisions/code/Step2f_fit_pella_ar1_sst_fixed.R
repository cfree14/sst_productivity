
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
library(freeR)
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Directories - original submission
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/code"

# Directories - revisions
tmbdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/code/tmb_code"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisions/output"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/revisionsfigures"
setwd(tmbdir)

# Read data
load(paste(datadir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))

# Read starting values
starts <- read.csv(file.path(datadir, "sp_ar1_starting_values.csv"), as.is=T)

# Source helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# Which SST, model, and shape parameter?
# SST choices: cobe, er, had
# Model choices: final, n1, n2, n3
# Shape parameters: 0.55, 0.20, 0.01
sst <- "cobe"
model <- "final"
p <- 0.55

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
  outputfile <- paste0("ramldb_v3.8_sp_ar1_pella", roundf(p,2), "_", sst, "_fixed.Rdata")
}else{
  sst.data <- as.matrix(data[,paste0(sst, "_sst_c_", model, "_sd")])
  outputfile <- paste0("ramldb_v3.8_sp_ar1_pella", roundf(p,2), "_", sst, "_", model, "_fixed.Rdata")
}


# Fit production model
################################################################################

# Compile TMB code
# Only run once to compile code
if(FALSE){
  dyn.unload(paste(tmbdir, dynlib("spmodel_ar1_pella_sst_fixed"), sep="/"))
  file.remove(paste(tmbdir, c("spmodel_ar1_pella_sst_fixed.o", "spmodel_ar1_pella_sst_fixed.dll"), sep="/"))
  compile("spmodel_ar1_pella_sst_fixed.cpp")
}

# Load TMB code
dyn.load(dynlib("spmodel_ar1_pella_sst_fixed"))

# Input data and parameter starting values
params <- list(ln_B0=log(pmin(starts$k, 1.5)),
               ln_r=log(pmin(starts$r, 0.8)),
               BetaT=rep(0.0, nstocks),
               ln_sigmaP=log(starts$sigmaP),
               rho=starts$rho)
input.data <- list(Nstocks=nstocks,
                   Nobs=nrow(data),
                   p=p,
                   StockID=as.factor(data$stockid),
                   B_t=data$tb_sd,
                   P_t=data$sp_sd,
                   Temp_t=sst.data)

# Initialization
model <- MakeADFun(data=input.data, parameters=params, random="BetaT", DLL="spmodel_ar1_pella_sst_fixed")
# model$control <- list(trace=1, parscale=rep(1,13), REPORT=1, reltol=1e-12, maxit=100)
# model$hessian <- F
# newtonOption(model, smartsearch=TRUE)

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
