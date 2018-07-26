
# Clear workspace
rm(list = ls())

# Turn off scientific notation
options(scipen=999)

# Read data
################################################################################

# Packages
library(freeR)
library(plyr)
library(dplyr)
library(forecast)

# Directories
datadir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/input"
plotdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/figures"
sstdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/sst/data/averages"
ramdir <- "/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/data/ramldb/ramldb_v3.8"
  
# RAMLDB stock and time series data
ram.stocks <- read.csv(paste(datadir, "all_stocks_in_ramldb.csv", sep="/"), as.is=T)
ram.ts.vals <- read.csv(paste(ramdir, "ramldb_v38_timeseries_values_views.csv", sep="/"), as.is=T)
ram.ts.units <- read.csv(paste(ramdir, "ramldb_v38_timeseries_units_views.csv", sep="/"), as.is=T)
colnames(ram.ts.vals) <- tolower(colnames(ram.ts.vals))

# SST time series data
cobe <- read.csv(paste(sstdir, "ramldb_sst_yearly_cobe.csv", sep="/"), as.is=T)
ersst <- read.csv(paste(sstdir, "ramldb_sst_yearly_ersst.csv", sep="/"), as.is=T)
hadsst <- read.csv(paste(sstdir, "ramldb_sst_yearly_hadsst.csv", sep="/"), as.is=T)


# Format data
################################################################################

# Reduce to candidate stocks
ram_use <- ram.stocks %>% 
  filter(sp_model=="yes") %>% 
  mutate(method_abbrev=revalue(method, c("Biomass dynamics model"="BDM", 
                                         "Integrated Analysis"="IA",
                                         "Statistical catch at age model"="SCAA",
                                         "Statistical catch at length model"="SCAL",
                                         "Survey index"="Survey")))
sort(unique(ram_use$method))
sort(unique(ram_use$method_abbrev))

# Merge SST data
sst <- cobe %>% 
  rename(cobe_sst_c=sst_c) %>% 
  # Add ERSST
  left_join(ersst, by=c("assessid", "year")) %>% 
  rename(er_sst_c=sst_c) %>% 
  # Add HadISST
  left_join(hadsst, by=c("assessid", "year")) %>% 
  rename(had_sst_c=sst_c)


# Helper functions
################################################################################

# Function to calculate surplus production: 
# SP(t) = TB(t+1) - TB(t) + C(t)
# tb <- subset(bdata, assessid==unique(bdata$assessid)[1])$tb
# catch <- subset(bdata, assessid==unique(bdata$assessid)[1])$catch
calc_sp <- function(tb, catch){
  sp <- c(tb[2:length(tb)] - tb[1:(length(tb)-1)] + catch[1:(length(catch)-1)], NA)
  return(sp)
}
  
# Function to fit surplus production model
# sp <- subset(data, assessid==unique(data$assessid)[1])$sp
# tb <- subset(data, assessid==unique(data$assessid)[1])$tb
fit.sp.model <- function(sp, tb){
  r_start <- log(0.4)
  k_start <- log(max(tb, na.rm=T) * 1.5)
  spfit <- try(nls(sp ~ exp(r)*tb*(1-tb/exp(k)),
                   start=list(r=r_start, k=k_start)))
  return(spfit)
}


# Build data
################################################################################

# Candidate stocks
# Identified in building RAM v3.8 stock key
stocks <- ram.stocks$stockid[ram.stocks$sp_model=="yes" & !is.na(ram.stocks$sp_model)]

# Assemble time series data
bdata <- ram.ts.vals %>%
  # Subset data
  filter(stockid%in%stocks) %>%
  select(assessid, stockid, year, tb, tc, tl, r, b.bmsytouse, u.umsytouse) %>%
  rename(bbmsy=b.bmsytouse, ffmsy=u.umsytouse) %>% 
  # Add catch type
  left_join(select(ram.stocks, assessid, catch_type), by="assessid") %>%
  mutate(catch=ifelse(catch_type=="tc", tc, tl)) %>%
  # Calculate surplus production
  group_by(assessid) %>%
    mutate(sp=calc_sp(tb, catch)) %>%
  # Add SST data
  left_join(sst, by=c("assessid", "year")) %>% 
  # Rearrange columns
  select(assessid, stockid, year, bbmsy, ffmsy,
         catch_type, catch, tb, sp, r, cobe_sst_c, er_sst_c, had_sst_c)

# Confirm that years are 100% sequential (no year is ever skipped)
# If a year got skipped, the calculations of surplus production would be incorrect
year.check <- bdata %>% 
  group_by(stockid) %>% 
  summarize(check=mean(year[2:length(year)]-year[1:(length(year)-1)]))
sum(year.check$check!=1)

# Check for completeness
apply(bdata, 2, function(x) sum(is.na(x)))


# Build trimming key
################################################################################

# Stocks with too strong SP or SR relationships
sp.stocks <- c("ATHAL5YZ", "BGROCKPCOAST", "BIGEYEATL",  "BLACKOREOWECR", 
               "BLKMARLINIO", "COWCODSCAL", "CRLOBSTERSA12", "CRLOBSTERSA34", 
               "CRLOBSTERSA56", "CRLOBSTERSA7", "CRLOBSTERSA8",
               "GRSPROCKNCAL", "GRSPROCKSCAL", "LNOSESKAPCOAST", "OROUGHYCASCADE", 
               "OROUGHYSE", "SAABALONESA", "SBWHITARGS", "SPSDOGPCOAST",
               "SWORDNATL", "SWORDSATL", "YEYEROCKPCOAST", "YFINATL")

# Inspect assessment method for strong assumption stocks
ram.stocks$method[ram.stocks$stockid%in%sp.stocks]

# Identify trim years
# Columns: stock id, biomass year, recruitment year, catch year, 
trim.years <- matrix(data=c("ARFLOUNDPCOAST", 1939, 1965, NA,
                            "ATHAL5YZ", 1898, NA, NA,
                            "BKINGCRABPI", 1982, NA, NA,
                            "BLACKROCKNPCOAST", 1943, 1968, 1965,
                            "BLACKROCKSPCOAST", 1943, 1968, NA,
                            "BLUEROCKCAL", 1940, 1960, 1940,
                            "BOCACCSPCOAST", NA, 1955, NA,
                            "CABEZORECOAST", NA, 1979, NA,
                            "CABEZSCAL", NA, 1969, NA,
                            "CHAKESA", NA, 1984, NA,
                            "CHILISPCOAST", NA, 1963, NA,
                            "COWCODSCAL", NA, NA, 1918,
                            "CRLOBSTERSA12", NA, 1980, NA,
                            "CRLOBSTERSA34", NA, 1980, NA,
                            "CRLOBSTERSA56", NA, 1980, NA,
                            "CRLOBSTERSA7", NA, 1980, NA,
                            "CRLOBSTERSA8", NA, 1980, NA,
                            "CROCKPCOAST", 1943, 1961, NA,
                            "CTRACSA", NA, 1984, NA,
                            "DEEPCHAKESA", 1955, 1983, 1955,
                            "DKROCKPCOAST", NA, 1975, NA,
                            "DSOLEPCOAST", NA, 1960, NA,
                            "ESOLEPCOAST", 1918, 1938, 1918,
                            "GEMFISHNZ", NA, 1979, NA,
                            "GEMFISHSE", 1968, NA, NA,
                            "GOPHERSPCOAST", NA, 1980, NA,
                            "GRNSTROCKPCOAST", 1942, 1969, NA,
                            "GRSPROCKNCAL", 1942, NA, NA,
                            "LINGCODNPCOAST", 1940, 1960, NA,
                            "LINGCODSPCOAST", NA, 1971, NA,
                            "LNOSESKAPCOAST", 1950, NA, 1950,
                            "MORWONGSE", NA, 1942, NA,
                            "NZLINGESE", 1977, 1983, 1977,
                            "NZLINGWSE", 1982, 1982, NA,
                            "NZSNAPNZ8", NA, 1970, NA,
                            "OROUGHYSE", 1987, NA, 1987,
                            "PGEODWA", 1995, NA, NA,
                            "PSOLEPCOAST", NA, 1939, NA,
                            "RSNAPGM", NA, 1971, NA,
                            "RSNAPSATLC", NA, 1973, NA,
                            "SABLEFPCOAST", NA, 1965, NA,
                            "SBELLYROCKPCOAST", NA, 1960, NA,
                            "SNROCKPCOAST", 1950, 1960, NA,
                            "SPANMACKGM", NA, 1984, NA,
                            "SPSDOGPCOAST", 1938, 1942, 1938,
                            "SSTHORNHPCOAST", NA, 1984, NA,
                            "STFLOUNNPCOAST", NA, 1984, NA,
                            "STFLOUNSPCOAST", NA, 1978, NA,
                            "SWHITSE", NA, 1981, NA,
                            "VSNAPSATLC", NA, 1975, 1960,
                            "WMARLINATL", NA, 1977, NA,
                            "WROCKPCOAST", 1941, 1946, 1941), ncol=4, byrow=T)
trim.years.df <- as.data.frame(trim.years, stringsAsFactors=F)
trim.years <- trim.years.df %>% 
  rename(stockid=V1, year_tb=V2, year_r=V3, year_catch=V4) %>% 
  mutate(year_tb=as.numeric(year_tb),
         year_catch=as.numeric(year_catch),
         year_trim=pmax(year_tb, year_r, year_catch, na.rm=T))
trim.years$stockid[!(trim.years$stockid%in%ram.stocks$stockid)] # must report 0 = confirms stock IDs spelled correctly
trim.stocks <- ram.stocks$assessid[ram.stocks$stockid%in%trim.years$stockid]


# Mark usable years
################################################################################

# Mark usable years
# If the stockid is in "trim.years", mark usable years
# If the stockid isn't in "trim.years", all years are usable
bdata$use <- NA
for(i in 1:length(stocks)){
  stock <- stocks[i]
  sdata <- subset(bdata, stockid==stock)
  stockid <- unique(sdata$stockid)
  if(stockid%in%trim.years$stockid){
    trim.year <- trim.years$year_trim[trim.years$stockid==stockid]
    usable.years <- trim.year:max(sdata$year)
    bdata$use[bdata$stockid==stock] <- ifelse(bdata$year[bdata$stockid==stock] %in% usable.years, "yes", "no")
  }else{
    bdata$use[bdata$stockid==stock] <- "yes"
  }
}
bdata$use <- as.factor(bdata$use)

# BIGHTREDSE is a special case (bad years are at end)
bdata$use[bdata$stockid=="BIGHTREDSE" & bdata$year>1990] <- "no"


# Plot trimming
################################################################################

# Data frame to record SP fits
spfits <- bdata %>%
  group_by(assessid, stockid) %>%
  summarize(tb_max=NA, r=NA, k=NA)

# For headers
top.i <- seq(1, length(stocks), 6)

# Plot data and trimming decisions
figname <- "AppendixC_ramldb_data_and_trimming.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(6, 5), mar=c(1.5, 1.0, 2.5, 0.5), mgp=c(2.5,0.5,0), oma=c(3,3,3,3), lwd=0.8)

# Loop through stocks
for(i in 1:length(stocks)){
  
  # Subset data
  stock <- stocks[i]
  sdata <- subset(bdata, stockid==stock & (!is.na(tb) | !is.na(catch)))
  print(paste(i, stock))
  
  # Format data
  sdata <- sdata %>% 
    mutate(tb=tb/1000,
           catch=catch/1000,
           sp=sp/1000,
           r=r/1E6)
  method <- ram_use$method_abbrev[ram_use$stockid==stock]
  
  # Year limits
  tmin <- floor(min(sdata$year)/10)*10
  tmax <- ceiling(max(sdata$year)/10)*10 
  
  # Biology limits
  bmin <- floor(min(sdata$tb, na.rm=T))
  bmax <- ceiling(max(sdata$tb, na.rm=T))
  smin <- floor(min(sdata$sp, na.rm=T))
  smax <- ceiling(max(sdata$sp, na.rm=T))
  cmin <- floor(min(sdata$catch, na.rm=T))
  cmax <- ceiling(max(sdata$catch, na.rm=T))
  rmin <- floor(min(sdata$r, na.rm=T))
  rmax <- ceiling(max(sdata$r, na.rm=T))
  if(rmax==1){rmax <- ceiling1(max(sdata$r, na.rm=T), 0.1)}
  if(rmax==0.1){rmax <- ceiling1(max(sdata$r, na.rm=T), 0.01)}
  
  # A. Plot biomass
  #######################################
  
  # Plot data
  plot(tb ~ year, sdata, type="l", bty="n", 
       xlim=c(tmin, tmax), ylim=c(bmin, bmax), yaxt="n",
       xlab="", ylab="", las=2, col="blue")
  axis(2, at=c(bmin, bmax))
  title(stock, line=0.1, xpd=NA, col.main=ifelse(stock%in%sp.stocks, "red", "black"))
  if(i %in% top.i){title("Biomass (1000s mt)\ntime series", col.main="blue", line=2, xpd=NA)}
  
  # Add trim year
  trim.yr <- trim.years$year_tb[trim.years$stockid==stock]
  if(length(trim.yr)>0){
    lines(x=c(trim.yr,trim.yr), y=c(bmin,bmax), col="black", lty=3, lwd=1.5)
    text(x=trim.yr, y=bmin+(bmax-bmin)*0.10, pos=2, labels=trim.yr, col="black", cex=1.2, xpd=NA)
  }

  # B. Plot recruitment
  #######################################
  
  # Plot data
  if(sum(sdata$r, na.rm=T)>0){
    plot(r ~ year, sdata, type="l", bty="n", 
         xlim=c(tmin, tmax), ylim=c(rmin, rmax), yaxt="n",
         xlab="", ylab="", las=2, col="purple")
    axis(2, at=c(rmin, rmax))
  }else{
    plot.new()
  }
  if(i %in% top.i){title("Recruitment (millions)\ntime series", col.main="purple", line=2, xpd=NA)}
  
  # Add trim year
  trim.yr <- trim.years$year_r[trim.years$stockid==stock]
  if(length(trim.yr)>0){
    lines(x=c(trim.yr,trim.yr), y=c(rmin,rmax), col="black", lty=3, lwd=1.5)
    text(x=trim.yr, y=rmin+(rmax-rmin)*0.10, pos=2, labels=trim.yr, col="black", cex=1.2, xpd=NA)
  }
  
  # C. Plot catch
  #######################################
  
  # Plot data
  plot(catch ~ year, sdata, type="l", bty="n", 
       xlim=c(tmin, tmax), ylim=c(cmin, cmax), yaxt="n",
       xlab="", ylab="", las=2, col="darkgreen")
  axis(2, at=c(cmin, cmax))
  if(i %in% top.i){title("Catch (1000s mt)\ntime series", col.main="darkgreen", line=2, xpd=NA)}
  
  # Add trim year
  trim.yr <- trim.years$year_catch[trim.years$stockid==stock]
  if(length(trim.yr)>0){
    lines(x=c(trim.yr,trim.yr), y=c(cmin,cmax), col="black", lty=3, lwd=1.5)
    text(x=trim.yr, y=cmin+(cmax-cmin)*0.10, pos=2, labels=trim.yr, col="black", cex=1.2, xpd=NA)
  }
  
  # D. Plot surplus production
  #######################################

  # Plot data
  plot(sp ~ tb, sdata, type="p", bty="n", 
       xlim=c(bmin, bmax), ylim=c(smin, smax), xaxt="n", yaxt="n",
       xlab="", ylab="", col=c("red", "grey60")[as.factor(sdata$use)])
  axis(1, at=c(bmin, bmax))
  axis(2, at=c(smin, smax))
  if(i %in% top.i){title("Surplus production\n(1000s mt) curve", col.main="black", line=2, xpd=NA)}
  if(stock%in%sp.stocks){text(x=bmin, y=smax, pos=4, labels="*", cex=3, col="red", xpd=NA)}
  text(x=bmax, y=smax, pos=2, labels=method, cex=1.1, col="black", xpd=NA)
  
  # Fit and plot SP model
  sdata1 <- subset(sdata, use=="yes")
  spfit <- fit.sp.model(sp=sdata1$sp, tb=sdata1$tb)
  if(!(inherits(spfit, "try-error"))){
    r <- exp(coef(spfit)["r"])
    k <- exp(coef(spfit)["k"])
    spfits$r[spfits$stockid==stock] <- r
    spfits$k[spfits$stockid==stock] <- k
    spfits$tb_max[spfits$stockid==stock] <- max(sdata1$tb, na.rm=T)
    curve(r*x*(1-x/k), from=bmin, to=bmax, add=T, lty=1, lwd=0.9, col="black")
  }

  # E. Plot stock-recruit relationship
  #######################################
  
  # Plot stock-recruit relationship
  if(sum(sdata$r, na.rm=T)>0){
    plot(r ~ tb, sdata, type="p", bty="n", 
         xlim=c(bmin, bmax), ylim=c(rmin, rmax), xaxt="n", yaxt="n", col=c("red", "grey60")[as.factor(sdata$use)])
    axis(1, at=c(bmin, bmax))
    axis(2, at=c(rmin, rmax))
  }else{
    plot.new()
  }
  if(i %in% top.i){title("Stock-recruit\nrelationship", col.main="black", line=2, xpd=NA)}
  
}

# Off
dev.off()
graphics.off()


# Inspect SP fits
################################################################################

# Format SP fits
spfits <- spfits %>% 
  mutate(k_scalar=k/tb_max)

# Histograms
summary(spfits$r)
summary(spfits$k_scalar)
hist(spfits$r, breaks=seq(0,2.2,0.1))
hist(spfits$k_scalar, breaks=seq(0,30,0.5))

# Number not fit
sum(is.na(spfits$r))


# Filter data
################################################################################

# Remove:
# (1) Stocks with strong SP/SR relationships
# (2) Stocks with fewer than 20 years of observations after trimming
# (3) Stocks without SST data

# 1. Remove stocks with strong SP/SR relationships
#########################################################

# Remove stocks with strong SP/SR relationship
length(stocks)
length(sp.stocks)
data <- bdata %>% 
  filter(!stockid%in%sp.stocks)
length(unique(data$stockid))

# 2. Remove stocks with <20 years after trimming
#########################################################

# Remove unusable years
data1 <- data %>% 
  filter(use=="yes" & !is.na(sp) & !is.na(tb))

# Inspect sample size
nyears <- data1 %>%
  group_by(stockid) %>%
  summarize(nyr=n()) %>% 
  arrange(nyr)
stocks.wout.20yr <- subset(nyears, nyr<20)$stockid
length(stocks.wout.20yr)

# Remove stocks with <20 years
data2 <- data1 %>% 
  filter(!stockid%in%stocks.wout.20yr)
length(unique(data2$stockid))


# 3. Remove stocks missing SST data
#########################################################

# Which stocks are missing SSTs?
stocks.wout.sst <- unique(subset(data2, is.na(cobe_sst_c))$stockid)
stocks.wout.sst.df <- subset(data2, stockid%in%stocks.wout.sst)
length(stocks.wout.sst)

# Remove stocks missing SST data
data3 <- data2 %>%
  filter(!stockid%in%stocks.wout.sst)
length(unique(data3$stockid))


# Standardize data
################################################################################

# Standardize TB, catch, surplus production data
# Standarize SST data and add SST null models
final <- data3 %>%
  group_by(stockid) %>% 
  mutate(tb_sd=tb/max(tb, na.rm=T),
         catch_sd=catch/max(tb, na.rm=T), 
         sp_sd=sp/max(tb, na.rm=T),
         cobe_sst_c_sd=cobe_sst_c-mean(cobe_sst_c),
         er_sst_c_sd=er_sst_c-mean(er_sst_c),
         had_sst_c_sd=had_sst_c-mean(had_sst_c))

# Confirm TB standardization
data.frame(tapply(final$tb_sd, final$stockid, max, na.rm=T)) # must be 1

# Confirm SST standardizations
data.frame(tapply(final$cobe_sst_c_sd, final$stockid, mean, na.rm=T)) # must be 0
data.frame(tapply(final$er_sst_c_sd, final$stockid, mean, na.rm=T)) # must be 0
data.frame(tapply(final$had_sst_c_sd, final$stockid, mean, na.rm=T)) # must be 0

# Convince yourself that SP ~ TB is really the same as SP(sd) ~ TB(sd)
stock <- unique(final$stockid)[4]
par(mfrow=c(1,2))
plot(sp ~ tb, final, subset=stockid==stock, bty="n", las=1)
plot(sp_sd ~ tb_sd, final, subset=stockid==stock, bty="n", las=1)

# Inspect completeness
complete(final)


# Export final dataset
################################################################################

# Export data
write.csv(final, paste(datadir, "ramldb_v3.8_production_data_clean.csv", sep="/"), row.names=F)
save(final, spfits, file=paste(datadir, "ramldb_v3.8_production_data_clean.Rdata", sep="/"))

