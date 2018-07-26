
# Function to fit SST-linked SP model
# For testing: id_t <- sdata$stockid; p_t <- sdata$pt_sim; b_t <- sdata$bt_sim; sst_t <- sdata$sst_obs
fit_sp_sst <- function(id_t, p_t, b_t, sst_t){
  
  # Merge and standardize data
  data <- data.frame(stockid=id_t, sp=p_t, b=b_t, sst=sst_t, stringsAsFactors=F) %>% 
    filter(!is.na(sp)) %>% 
    group_by(stockid) %>% 
    mutate(sp_sd=sp/max(b),
           b_sd=b/max(b),
           sst_sd=sst-mean(sst)) 
  
  # Directories
  tmbdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/spmodel_tb/code/tmb_code"
  setwd(tmbdir)
  
  # Load TMB code
  dyn.load(dynlib("spmodel_sst_link"))
  
  # Model parameters
  stocks <- unique(data$stockid)
  nstocks <- length(stocks)
  
  # Input data and parameter starting values
  params <- list(ln_B0=rep(log(1.5), nstocks),
                 ln_r=rep(log(0.4), nstocks),
                 BetaT=rep(0.0, nstocks),
                 ln_sigmaP=rep(-2.5, nstocks), # 0 before, -2.5 based on model fits
                 mu_T=0.0,
                 ln_sd_T=-1.25) # -3 before, -1.25 based on model fits
  input.data <- list(Nstocks=nstocks,
                     Nobs=nrow(data),
                     StockID=as.factor(data$stockid),
                     B_t=data$b_sd,
                     P_t=data$sp_sd,
                     Temp_t=data$sst_sd)
  
  # Initialize model
  model <- MakeADFun(data=input.data, parameters=params, random="BetaT", DLL="spmodel_sst_link")
  
  # Fit model
  output <- TMBhelper::Optimize(obj=model, lower=-Inf, upper=Inf, loopnum=3, newtonsteps=3, bias.correct=FALSE, getsd=FALSE)
  
  # Calculate SD
  sd <- try(sdreport(model, hessian.fixed=hess))
  
  
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
  
  # Return
  return(results.wide)

}
