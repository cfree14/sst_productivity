
### Impact of historical warming on marine fisheries production

This GitHub repository contains the data, code, output, tables, and figures for:

Free CM, Thorson JT, Pinsky ML, Oken KL, Wiedenmann J, Jensen OP (2019) Impacts of historical warming on marine fisheries production. *Science* 363(6430): 979–983. [link](http://science.sciencemag.org/content/363/6430/979)

The majority of the files of interest are in the "revisions" folder in the base directory which is organized using the following architecture:

1. code......... folder containing code to fit models, perform hindcasts, build figures and tables
    + figures.......... folder containing R scripts to plot figures and build tables
    + tmb_code......... folder containing TMB models called by the R scripts in the "code" directory
    + 232_stocks....... IGNORE THIS FOLDER (GARBAGE RESULTS)
    + ... and R scripts for fitting the various models and performing the hindcasts
2. figures...... folder containing figures for the manuscript (and some garbage figures)
3. output....... folder containing output from model fits, hindcasts, and other results
4. tables....... folder containing tables for the manuscript (and some garbage tables)

The only files of interest that are not in the "revisions" folder are the following:

1. input....... folder containing analyzed data (RAM Legacy biomass + surplus production and observed and simulated SST time series)
2. code........ folder containing code to assemble data analyzed in manuscript
    + Step0a_build_stock_data.R
    + Step0b_build_life_history_data.R
    + Step1a_build_spmodel_data.R
    + Step1b_add_sst_null_models.R

Everything else in these folders is from the pre-review analysis. The analysis improved significantly through the revision process and the majority of the up-to-date code, data, output, tables, figures, etc are in the "revisions" folder.

The Rdata files containing the bootstrapped temperature-driven MSY trajectories are very large (>500MB) and are not included in the repository. Please email Chris Free (cfree14 AT gmail.com) to request access to these files.
