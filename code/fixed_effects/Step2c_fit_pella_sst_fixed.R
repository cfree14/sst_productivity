
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

# Directories
tmbdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/code/tmb_code"
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/code"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
setwd(tmbdir)

# Read data
load(paste(datadir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))

# Read NLS fits
spfits <- read.csv(paste(outputdir, "ramldb_v3.8_pella_0.2_cobe_nls.csv", sep="/"), as.is=T)

# Source helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# Which SST, model, and shape parameter?
# Shape parameters: 0.55 (45%), 0.20 (40%), 0.01 (37%)
# SST choices: cobe, er, had
# Model choices: final, n1, n2, n3
p <- 0.20
sst <- "cobe"
model <- "final"

# Format data
################################################################################

# Inspect data completeness
apply(data, 2, function(x) sum(is.na(x)))

# Remove problem stocks
stocks.nls.no <- spfits$stockid[is.na(spfits$betaT)]
stocks.w.big.k <- c("AMPL3LNO", "ARGHAKESARG", "BRMSOJ", "CHTRACCH", "CMACKPJPN", 
                    "COD3Pn4RS", "CODBA2532", "CODIS", "CODVIa", "CTRACSA", 
                    "FMEG8c9a", "GEMFISHSE", "GHALNEAR", "HERRNIRS", "JMACKPJPN",
                    "JMACKTSST", "MEG8c9a", "MORWONGSE", "PATGRENADIERCH", 
                    "PATGRENADIERSARG", "RBRMECS", "REXSOLEGA", "RKCRABBB", "RSNAPGM",
                    "SNAPSAUSNSG", "SOLEVIII", "WPOLLNSO", "WPOLLWBS")
stocks.w.big.r <- c("ACADRED2J3K", "ACADREDUT3", "COD3Ps", "LISQUIDATLC", 
                    "SEELNSSA1", "SKJEATL", "SKJWATL", "HERRTOG", "SEELNSSA2", "SEELNSSA2")
problem.stocks <- c(stocks.nls.no, stocks.w.big.k, stocks.w.big.r)
data <- subset(data, !(stockid%in%problem.stocks))


# Select data and model
################################################################################

# Parameters
stocks <- unique(data$stockid)
nstocks <- length(stocks)

# Specify SST data to use and output name
if(model=="final"){
  sst.data <- as.matrix(data[,paste0(sst, "_sst_c_sd")])
  outputfile <- paste0("ramldb_v3.8_pella_", roundf(p,2), "_", sst, "_fixed.Rdata")
}else{
  sst.data <- as.matrix(data[,paste0(sst, "_sst_c_", model, "_sd")])
  outputfile <- paste0("ramldb_v3.8_pella_", roundf(p,2), "_", sst, "_", model, "_fixed.Rdata")
}


# Fit production model
################################################################################

# Compile TMB code
# Only run once to compile code
if(FALSE){
  dyn.unload(paste(tmbdir, dynlib("pella_sst_fixed"), sep="/"))
  file.remove(paste(tmbdir, c("pella_sst_fixed.o", "pella_sst_fixed.dll"), sep="/"))
  compile("pella_sst_fixed.cpp")
}

# Load TMB code
dyn.load(dynlib("pella_sst_fixed"))

# Input data and parameter starting values
b0_starts <- log(sapply(stocks, function(x) max(data$tb_sd[data$stockid==x])) * 1.5)
params <- list(ln_B0=b0_starts,
               ln_r=rep(log(0.4), nstocks),
               BetaT=rep(0.0, nstocks),
               ln_sigmaP=rep(-2.5, nstocks)) # -3 before, -1.25 based on model fits
input.data <- list(Nstocks=nstocks,
                   Nobs=nrow(data),
                   p=p,
                   StockID=as.factor(data$stockid),
                   B_t=data$tb_sd,
                   P_t=data$sp_sd,
                   Temp_t=sst.data)

# Initialization
model <- MakeADFun(data=input.data, parameters=params, DLL="pella_sst_fixed")
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
