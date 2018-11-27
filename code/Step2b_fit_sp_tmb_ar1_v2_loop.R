
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


# Setup model
################################################################################

# Compile TMB code
# Only run once to compile code
if(FALSE){
  dyn.unload(paste(tmbdir, dynlib("spmodel_ar1_v2"), sep="/"))
  file.remove(paste(tmbdir, c("spmodel_ar1_v2.o", "spmodel_ar1_v2.dll"), sep="/"))
  compile("spmodel_ar1_v2.cpp")
}

# Load TMB code
dyn.load(dynlib("spmodel_ar1_v2"))


# Fit production model
################################################################################

# Parameters
stocks <- unique(data$stockid)
nstocks <- length(stocks)

# Setup container
df <- data.frame(stockid=stocks, r=NA, k=NA, sigmaP=NA, rho=NA, stringsAsFactors=F)

# Loop through stocks
for(i in 1:nrow(df)){
  
  # Subset data
  stock <- df$stockid[i]
  sdata <- filter(data, stockid==stock)
  print(paste(i, stock))

  # Input data and parameter starting values
  params <- list(ln_B0=log(1.5),
                 ln_r=-1,
                 ln_sigmaP=0,
                 rho=0.1)
  input.data <- list(Nstocks=1,
                     Nobs=nrow(sdata),
                     StockID=as.factor(sdata$stockid),
                     B_t=sdata$tb_sd,
                     P_t=sdata$sp_sd)

  # Initialize model
  model <- MakeADFun(data=input.data, parameters=params, DLL="spmodel_ar1_v2")
  model$control <- list(trace=1, parscale=rep(1,13), REPORT=1, reltol=1e-12, maxit=100)
  model$hessian <- F
  newtonOption(model, smartsearch=TRUE)
  
  # Run model
  output <- try(TMBhelper::Optimize(obj=model, lower=c(-Inf,-Inf,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf), loopnum=3, newtonsteps=3, bias.correct=FALSE, getsd=FALSE))
  
  # Record results
  if(!inherits(output, "try-error")){
    hess <- optimHess(par=output$par, fn=model$fn, gr=model$gr)
    sd <- try(sdreport(model, hessian.fixed=hess))
    mat <- summary.sdreport(sd)
    df$r[i] <- mat["r","Estimate"]
    df$k[i] <- mat["B0","Estimate"]
    df$sigmaP[i] <- mat["sigmaP","Estimate"]
    df$rho[i] <- mat["rho","Estimate"]
  }

}

# Histograms
hist(df$rho)
hist(df$r)
hist(df$k)
hist(df$sigmaP)

# Export values
write.csv(df, file.path(datadir, "sp_ar1_starting_values.csv"), row.names=F)
