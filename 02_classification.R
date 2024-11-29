###-###-###-###-###-###-###-###-###
#
# Classify time series
#
# 02_classification.R
#
###-###-###-###-###-###-###-###-###

# Load packages, functions, and data --------------------------------------
source("00_Rpackages.R")
source("R/fun01_classification.R")
source("R/fun02_analyses.R")

# Load RAMLDB v4.61:
load("data/RAMLDB v4.61/R Data/DBdata[asmt][v4.61].RData")

# Compute and add surplus production timeseries to the data:
add_surplus()

# Prepare surplus production time series --------------------------------

len <- 25 # minimal timeseries length
run_loo <- TRUE # run Leave-One-Out cross validation
start_year <- 1950 # start year

# Stocks to discard because of biomass likely estimated with deterministic model
determ_model_stk <-
  c("ALBASATL", "BIGEYEATL", "BLTILESATLC", "DVANGSASC", "DVANGSAWC",
    "KINGKLIPSASC", "KINGKLIPSAWC", "PATCODARGS", "REDFISHSPP3LN",
    "SPURDNEATL", "SWORDSATL", "YELL3LNO", "YFINATL", "YNOSESKASCH")

# List stocks with time series available for classification
stock_list <- list_available(c("SProd"), min_len=0, # with no missing year
                             min_year=start_year) %>% # from 1950 onwards
  dplyr::filter(phylum == "Chordata") # keep only fish stocks
nrow(stock_list) # 397

# Discard deterministic models
stock_list <- stock_list %>%
  dplyr::filter(!stockid %in% determ_model_stk)
nrow(stock_list) # 383

# Discard time series below 25 points
stock_list <- stock_list %>%
  dplyr::filter(length>=len)
nrow(stock_list) # 315
median(stock_list$length) # 41
max(stock_list$length) # 71

stock_names <- stock_list$stockid

# Get SProd time series
df_SProd <- timeseries_values_views %>%
  dplyr::filter(stockid %in% stock_names,
                year >= start_year) %>%
  dplyr::group_by(stockid) %>%
  tidyr::drop_na(SProd)

# Normalize SProd by average total biomass
df_SProd_norm <- df_SProd %>%
  dplyr::select(stockid, year, SProd) %>%
  dplyr::left_join(
    df_SProd %>%
      dplyr::summarise(mean_TB = mean(TBbest)),
    by="stockid") %>%
  dplyr::mutate(SProd = SProd/mean_TB)

# Split into list of timeseries
list_SProd_norm <- df_SProd_norm %>%
  dplyr::group_split() %>%
  stats::setNames(stock_names)

# Save SProd timeseries
dir.create("data/ts_data/", showWarnings=FALSE)
readr::write_rds(list_SProd_norm, paste0("data/ts_data/SProd_RAMLDBv4.61_minlen",
                                         len,"_start",start_year,".rds"))

# Classify surplus production time series --------------------------------
list_SProd_norm <-
  readr::read_rds(paste0("data/ts_data/SProd_RAMLDBv4.61_minlen",
                         len,"_start",start_year,".rds"))

# Create directory to store classification output:
dir.create("res/classif/", showWarnings=FALSE)
dir1 <- paste0("res/classif/SProd_RAMLDBv4.61_minlen",
               len,"_start",start_year,"/")
dir.create(dir1, showWarnings=FALSE)

# Run the classification:
classif_SProd <-
  run_classif_data(
    df_list=list_SProd_norm, min_len=len,
    group="stockid", time="year", variable=c("SProd"),
    str="aic_asd", run_loo=run_loo, two_bkps=TRUE,
    makeplots=TRUE, ind_plot=NULL,
    dirname=dir1, save_plot=TRUE)

path1 <- paste0("classif_v4.61_SProd_minlen", len ,
                "_normTBavg_loo", run_loo,"_aicasd_start", start_year)

saveRDS(classif_SProd,
        file=paste0("res/classif/", path1, ".rds"))


# Classify surplus production time series based on AICc only ------------------
list_SProd_norm <-
  readr::read_rds(paste0("data/ts_data/SProd_RAMLDBv4.61_minlen",
                         len,"_start",start_year,".rds"))

# Create directory to store classification output:
dir.create("res/classif/", showWarnings=FALSE)
dir2 <- paste0("res/classif/SProd_RAMLDBv4.61_minlen",
               len,"_start",start_year,"_AICconly/")
dir.create(dir2, showWarnings=FALSE)

# Run the classification:
classif_SProd <-
  run_classif_data(
    df_list=list_SProd_norm, min_len=len,
    group="stockid", time="year", variable=c("SProd"),
    str="aic", run_loo=run_loo, two_bkps=TRUE,
    makeplots=TRUE, ind_plot=NULL,
    dirname=dir2, save_plot=TRUE)

path2 <- paste0("classif_v4.61_SProd_minlen", len ,
                "_normTBavg_loo", run_loo,"_aicasd_start", start_year,"_AICconly")

saveRDS(classif_SProd,
        file=paste0("res/classif/", path2, ".rds"))

