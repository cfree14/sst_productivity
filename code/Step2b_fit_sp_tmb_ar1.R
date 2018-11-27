
# Clear workspace
rm(list = ls())

# Read data
################################################################################

# Packages
library(TMB)
library(devtools)
library(reshape2)
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Directories
tmbdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/code/tmb_code"
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
codedir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/code"
outputdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output"
setwd(tmbdir)

# Source helper functions
source(paste(codedir, "helper_functions.R", sep="/"))

# Read data
load(paste(datadir, "ramldb_v3.8_production_data_final.Rdata", sep="/"))
data_orig <- data


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


# Fit production model
################################################################################

# Parameters
stocks <- unique(data$stockid)
nstocks <- length(stocks)

# Compile TMB code
# Only run once to compile code
if(FALSE){
  dyn.unload(paste(tmbdir, dynlib("spmodel_ar1"), sep="/"))
  file.remove(paste(tmbdir, c("spmodel_ar1.o", "spmodel_ar1.dll"), sep="/"))
  compile("spmodel_ar1.cpp")
}

# Load TMB code
dyn.load(dynlib("spmodel_ar1"))

# Input data and parameter starting values
phi_start <- 0.2
b0_starts <- log(sapply(stocks, function(x) max(data$tb_sd[data$stockid==x])) * 1.5) 
params <- list(ln_B0=b0_starts,
               ln_r=rep(-1.0, nstocks),
               ln_sigmaP=rep(0.0, nstocks),
               # logit_phi=rep(log(phi_start/(1-phi_start)), nstocks) # Kiva
               logit_phi=rep(0, nstocks) # Jim
               )
input.data <- list(Nstocks=nstocks,
                   Nobs=nrow(data),
                   StockID=as.factor(data$stockid),
                   B_t=data$tb_sd,
                   P_t=data$sp_sd)

# Initialization
model <- MakeADFun(data=input.data, parameters=params, DLL="spmodel_ar1")

# Run model
# A model with good fit should have a gradient near zero
# output <- nlminb(start=model$par, objective=model$fn, gradient=model$gr, lower=-Inf, upper=Inf, control=list(eval.max=1e4, iter.max=1e4))
output <- TMBhelper::Optimize(obj=model, lower=-Inf, upper=Inf, loopnum=3, newtonsteps=3, bias.correct=FALSE, getsd=FALSE)


# Model diagnostics
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
results.wide <- dcast(results.df, stockid ~ param, value.var="estimate")

# Inspect for bad fits
############################################

# Histogram parameter values
params <- unique(results.df$param)
par(mfrow=c(2,5))
for(i in 1:length(params)){
  vals <- results.df$estimate[results.df$param==params[i]]
  hist(vals, main="", xlab=params[i], ylab="")
}


# Save results
################################################################################

# Export model objects
outputfile <-"ramldb_v3.8_sp.Rdata"
save(data, stocks, problem.stocks, nstocks, 
     input.data, params, 
     model, output, sd, hess,
     results.df, results.wide,
     file=paste(outputdir, outputfile, sep="/"))
