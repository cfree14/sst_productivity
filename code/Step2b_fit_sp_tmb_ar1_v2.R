
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

# Read starting values
starts <- read.csv(file.path(datadir, "sp_ar1_starting_values.csv"), as.is=T)


# Format data
################################################################################

# Inspect data completeness
apply(data, 2, function(x) sum(is.na(x)))

# Remove problem stocks
stocks.nls.no <- spfits$stockid[is.na(spfits$r)]
stocks.w.big.k <- c("AMPL3LNO", "ARGHAKESARG", "BRMSOJ", "MORWONGSE", "RBRMECS", 
                    "REXSOLEGA", "RSNAPGM", "WPOLLNSO")
# stocks.w.flat.resids <- c("ACADRED2J3K", "AMPL23K", "COD3M", "COD3Pn4RS", "FLSOLEGA", "GHAL4RST",
#                           "KELPGREENLINGORECOAST", "PILCHPJPN", "PILCHTSST", "REDDEEP2J3K−3LNO",
#                           "SAILWATL", "SARDPCOAST", "SKJCIO", "SKJEATL", "SKJWATL", "SPANMACKSATLC", "WPOLLAI",
#                           "YELLSNEMATL")
# stocks.w.smooth.resids <- c("ALPLAICBSAI", "ARFLOUNDBSAI", "ARFLOUNDGA", "BIGHTREDSE", "BLACKROCKSPCOAST",
#                             "BLUEROCKCAL", "BSQLOBSTERCH", "CABEZNCAL", "CABEZORECOAST", "CABEZSCAL", "CHAKESA", "CTRACSA",
#                             "DEEPCHAKESA", "DEEPFLATHEADSE", "DKROCKPCOAST", "DSOLEPCOAST", "DUSROCKGA", "LSTHORNHPCOAST",
#                             "MONKGOMNGB", "MONKSGBMATL", "NROCKGA", "PTOOTHFISHMI", "RPORGYSATLC", "SABLEFPCAN", "SURFMATLC", "YSOLEBSAI")
problem.stocks <- c(stocks.nls.no, stocks.w.big.k)
data <- subset(data, !(stockid%in%problem.stocks))

# Reduce to 10 stocks
# ten <- unique(data$stockid)[1:50]
# data <- subset(data, stockid%in%ten)


# Fit production model
################################################################################

# Parameters
stocks <- unique(data$stockid)
nstocks <- length(stocks)

# Compile TMB code
# Only run once to compile code
if(FALSE){
  dyn.unload(paste(tmbdir, dynlib("spmodel_ar1_v2"), sep="/"))
  file.remove(paste(tmbdir, c("spmodel_ar1_v2.o", "spmodel_ar1_v2.dll"), sep="/"))
  compile("spmodel_ar1_v2.cpp")
}

# Load TMB code
dyn.load(dynlib("spmodel_ar1_v2"))

# Input data and parameter starting values
# params <- list(ln_B0=rep(log(1.5), nstocks),
#                ln_r=rep(-1.0, nstocks),
#                ln_sigmaP=rep(0.0, nstocks),
#                # rho=0.5) # rho for all stocks
#                rho=rep(0.1, nstocks)) # uncomment when fixing rho to 0
#                # logit_rho=rep(logit_rho_start, nstocks)) # comment out when fixing rho to 0
params <- list(ln_B0=log(pmin(starts$k, 1.5)),
               ln_r=log(pmin(starts$r, 0.8)),
               ln_sigmaP=log(starts$sigmaP),
               rho=starts$rho)
input.data <- list(Nstocks=nstocks,
                   Nobs=nrow(data),
                   StockID=as.factor(data$stockid),
                   B_t=data$tb_sd,
                   P_t=data$sp_sd)

# Set rho to 0
# map <- list()# uncomment when fixing rho to 0
# map[["rho"]] <- rep(NA, nstocks)# uncomment when fixing rho to 0
# map[["rho"]] <- factor(map[["rho"]])# uncomment when fixing rho to 0

# Initialization
model <- MakeADFun(data=input.data, parameters=params, DLL="spmodel_ar1_v2")
# model <- MakeADFun(data=input.data, parameters=params, DLL="spmodel_ar1_v2", map=map) # uncomment when fixing rho to 0
model$control <- list(trace=1, parscale=rep(1,13), REPORT=1, reltol=1e-12, maxit=100)
model$hessian <- F
newtonOption(model, smartsearch=TRUE)

# Run model
# A model with good fit should have a gradient near zero
output <- nlminb(start=model$par, objective=model$fn, gradient=model$gr, lower=-Inf, upper=Inf, control=list(eval.max=1e4, iter.max=1e4))
# output <- TMBhelper::Optimize(obj=model, lower=c(-Inf,-Inf,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf), loopnum=3, newtonsteps=3, bias.correct=FALSE, getsd=FALSE)


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
outputfile <-"ramldb_v3.8_sp_ar1_rho0.Rdata"
save(data, stocks, problem.stocks, nstocks, 
     input.data, params, 
     model, output, sd, hess,
     results.df, results.wide,
     file=paste(outputdir, outputfile, sep="/"))
