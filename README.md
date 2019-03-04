
### Impact of historical warming on marine fisheries production

This GitHub repository contains the data, code, output, tables, and figures for:

Free CM, Thorson JT, Pinsky ML, Oken KL, Wiedenmann J, Jensen OP (2019) Impacts of historical warming on marine fisheries production. *Science* 363(6430): 979–983. [[link]](http://science.sciencemag.org/content/363/6430/979)

I apologize for the slight disarray but the files of interest fall in two places in the repository. First, the analyzed dataset gets assembled by scripts in the "code" folder of the base directory and gets saved in the "input" folder of the base directory. The remainder of the analysis (i.e., model fitting, hindcasts, figure/table construction) is performed in the "revisions" folder of the base directory.

Here's a quick run down of the files that assemble the input data and where this data is saved:

1. code........ folder containing code to assemble data analyzed in manuscript
    + Step0a_build_stock_data.R
    + Step0b_build_life_history_data.R
    + Step1a_build_spmodel_data.R
    + Step1b_add_sst_null_models.R
2. input....... folder containing analyzed data (RAM Legacy biomass + surplus production and observed and simulated SST time series)

The remainder of the files are in the "revisions" folder which is organized according to the following architecture:

1. code......... folder containing code to fit models, perform hindcasts, build figures and tables
    + figures.......... folder containing R scripts to plot figures and build tables
    + tmb_code......... folder containing TMB models called by the R scripts in the "code" directory
    + ... and R scripts for fitting the various models and performing the hindcasts
2. figures...... folder containing figures for the manuscript (and some garbage figures)
3. output....... folder containing output from model fits, hindcasts, and other results
4. tables....... folder containing tables for the manuscript (and some garbage tables)

The Rdata files containing the bootstrapped temperature-driven MSY trajectories are very large (>500MB) and are not included in the repository. Please email Chris Free (cfree14 AT gmail.com) to request access to these files.
