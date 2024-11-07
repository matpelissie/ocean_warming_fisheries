###-###-###-###-###-###-###-###-###
#
# 04/09/2023 mathieu.pelissie@ens-lyon.fr
#
# Functions for joint pattern analyses
#
# functions_jointpattern.R
#
###-###-###-###-###-###-###-###-###


# I) Data preparation and selection -----------------------------------------

#' Compute surplus production in raw data frames and add it to the dataset
#'
#' @param
#'
#' @return No return value.
#'
#' @export

add_surplus <- function(){

  if(is.null(timeseries_values_views$SProd)){

    # Compute surplus production (SProd):
    timeseries_values_views <<- timeseries_values_views %>%
      dplyr::group_by(stockid) %>%
      dplyr::mutate(SProd = dplyr::lead(TBbest)-TBbest+TCbest) %>%
      dplyr::ungroup()

    # Compute exploitation rate (ER) because some are missing:
    timeseries_values_views <<- timeseries_values_views %>%
      dplyr::mutate(ERbest2 = TCbest/TBbest)

    # Add SProd to all time series data frame:
    timeseries <<- timeseries %>%
      dplyr::bind_rows(
        timeseries_values_views %>%
          dplyr::select(c(stockid, stocklong, year, SProd)) %>%
          dplyr::rename(tsyear = year) %>%
          tidyr::pivot_longer(cols=SProd,
                              names_to = "tsid",
                              values_to = "tsvalue",
                              values_drop_na = TRUE) %>%
          dplyr::mutate(tsid = "SProd-MT") %>%

          # Add best assessment for total biomass (TB):
          dplyr::left_join(timeseries_assessments_views %>%
                             dplyr::select(TBbest, stockid) %>%
                             dplyr::rename(assessid = TBbest),
                           by = "stockid")
      )

    # Add range of years available:
    SProd_years <<- timeseries_values_views %>%
      dplyr::select(stockid, stocklong, year, SProd) %>%
      dplyr::group_by(stockid) %>%
      tidyr::drop_na(SProd) %>%
      dplyr::summarise(SProd = paste0(min(year, na.rm=FALSE),"-",
                                      max(year, na.rm=FALSE)))

    timeseries_years_views <<- timeseries_years_views %>%
      dplyr::left_join(SProd_years, by="stockid")

  }

  return(invisible(NULL))

}



#' List stocks with time series available for classification
#'
#' @param ts_types List of types of empirical time series
#' (TCbest, TBbest, SProd...).
#' @param min_len Minimal length of time series to include (default is 25).
#' @param min_year Minimum year
#'
#' @return List of stocks with available timeseries as requested.
#'
#' @export

list_available <- function(ts_types, min_len=25, min_year=NULL) {

  if(!is.null(min_year) & is.numeric(min_year)){

    ts <- timeseries_values_views %>%
      dplyr::filter(year>=min_year)

  } else {ts <- timeseries_values_views}

  ts <- ts %>%
    # Add geographic and taxonomic info to allow filtering:
    dplyr::left_join(metadata, by="stockid") %>%
    dplyr::select(-c(FisheryType, taxGroup)) %>%
    dplyr::left_join(taxonomy %>%
                       dplyr::distinct(scientificname, .keep_all=TRUE),
                     by="scientificname") %>%

    dplyr::group_by(stockid) %>%
    dplyr::select(stockid, year, dplyr::all_of(ts_types)) %>%
    stats::na.omit(ts_types) %>%
    dplyr::summarise(first=min(year),
                     last=max(year)) %>%
    dplyr::mutate(length = last-first+1) %>%
    dplyr::mutate(scen = paste(paste(ts_types, collapse="."),
                               stockid, sep="_")) %>%
    dplyr::relocate(scen) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metadata, by="stockid") %>%
    dplyr::select(-c(FisheryType, taxGroup)) %>%
    dplyr::left_join(taxonomy %>%
                       dplyr::distinct(scientificname, .keep_all=TRUE),
                     by="scientificname") %>%

    # Discard stocks with timeseries with at least one gap:
    dplyr::left_join(gaps_info(ts_types) %>%
                       dplyr::select(stockid, gaps, missing),
                     by="stockid") %>%
    dplyr::mutate(gaps = ifelse(is.na(gaps), FALSE, TRUE)) %>%
    dplyr::filter(length>=min_len & !gaps)

  return(ts)
}



#' Specify stocks for which spatial boundaries are available
#' from https://chrismfree.com/ram-legacy-stock-boundary-database/
#'
#' @param
#'
#' @return No return value.
#'
#' @export

add_boundaries <- function(){

  if(is.null(timeseries_notes_views$boundaries_available)){

    bounds_meta <- readr::read_csv(
      paste0("data/_spatial/ramldb_boundaries/",
             "ramldb_v3.8_stock_boundary_table_v2_formatted.csv")) %>%
      dplyr::mutate(boundaries_available = TRUE)

    timeseries_notes_views <<- timeseries_notes_views %>%
      dplyr::left_join(bounds_meta %>%
                         dplyr::select(stockid, boundaries_available),
                       by="stockid") %>%
      dplyr::mutate(boundaries_available =
                      ifelse(is.na(boundaries_available), FALSE, TRUE))

    boundaries_added <<- TRUE
  }

}



#' Gives the number of time series with gaps
#'
#' @param ts_type Type(s) of empirical time series (TCbest, TBbest, SProd...).
#'
#' @return Data frame with timeseries with gaps and the number of missing years.
#' @export

gaps_info <- function(ts_type){

  df <- timeseries_values_views %>%
    dplyr::select(stockid, stocklong, year,
                  dplyr::all_of(ts_type)) %>%
    na.omit(ts_type) %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(first = min(year),
                  last = max(year),
                  length = last - first + 1) %>%
    dplyr::add_count(stockid) %>%
    dplyr::mutate(gaps = length > n)

  # is_gaps <- df %>%
  #   dplyr::summarise(is_gap = all(gaps)) %>%
  #   dplyr::summarise(with_gap = sum(is_gap))

  gap_details <- df %>%
    dplyr::filter(gaps) %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(missing = length - n) %>%
    dplyr::slice(1)

  return(gap_details)
}



#' Summarise the number of time series available by type
#'
#' @param min_len Minimal length of time series to include (default is 20).
#'
#' @return Data frame with the number of timeseries available.
#'
#' @export

ts_available_type <- function(min_len){

  ts_avail <- timeseries %>%
    tidyr::drop_na(tsvalue) %>%
    dplyr::group_by(stockid, tsid) %>%
    dplyr::mutate(n=n()) %>%
    dplyr::filter(n>=min_len) %>%
    dplyr::slice(1) %>%
    dplyr::group_by(tsid) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::arrange(desc(n))

    return(ts_avail)
}



#' Summarise the number of bioparameters available
#'
#' @return Data frame with the number of bioparameters available.
#'
#' @export

bioparam_available_type <- function(){

  bioparam_avail <- bioparams %>%
    dplyr::group_by(stockid, bioid) %>%
    dplyr::mutate(n=n()) %>%
    dplyr::slice_tail() %>%
    dplyr::group_by(bioid) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(bioshort = sub("-[^-]+$","",bioid)) %>%
    dplyr::left_join(biometrics %>%
                       dplyr::select(bioshort, biolong) %>%
                       dplyr::group_by(bioshort) %>%
                       dplyr::slice_tail(),
                     by="bioshort") %>%
    dplyr::select(-bioshort)

  return(bioparam_avail)
}


#' Extract all bioparameters of a given type
#'
#' @param param Name of the bioparameter to extract.
#'
#' @return Data frame with the bioparameter for all stocks available.
#'
#' @export

extract_bioparam <- function(param){

  bioparam_stocks <- bioparams %>%
    dplyr::filter(bioid == param) %>%
    dplyr::mutate(bioshort = sub("-[^-]+$","",bioid)) %>%
    dplyr::left_join(biometrics %>%
                       dplyr::select(bioshort, biolong) %>%
                       dplyr::group_by(bioshort) %>%
                       dplyr::slice_tail(),
                     by="bioshort")

  return(bioparam_stocks)
}



#' Define mature phase of fishery
#'
#' @param id Short stock ID character.
#'
#' @return Data frame with the first year of the mature phase, whether the
#' onset of the mature phase occurred before data are available (truncated),
#' or NA if information to define mature is missing.
#'
#' @export

mature_def <- function(id){

  # Criteria for mature fishery from Melnychuk et al. 2022
  mature_ts <- timeseries_values_views %>%
    dplyr::filter(stockid==id) %>%
    dplyr::mutate(mature_a = ifelse(BdivBmsypref<0.8 | BdivBmgtpref<0.8,
                                    TRUE, FALSE),
                  mature_b = ifelse(UdivUmsypref>1 | UdivUmgtpref>1,
                                    TRUE, FALSE),
                  mature_c = ifelse(CdivMSY>1, TRUE, FALSE),
                  mature_d = ifelse(CdivMEANC>1.25, TRUE, FALSE),
                  mature_f = ifelse(BdivBmsypref<1 | BdivBmgtpref<1,
                                    TRUE, FALSE)
    )

  # Preliminary computation for criterion (g):
  bef_assess <- mature_ts %>%
    dplyr::filter((dplyr::if_any(contains("BdivB"),
                                 function(x) !is.na(x)) &
                     dplyr::if_any(contains("UdivU"),
                                   function(x) !is.na(x)))) %>%
    dplyr::pull(year)

  if(length(bef_assess)>0){

    meanC_10y_df <- mature_ts %>%
      dplyr::filter(year < min(bef_assess)) %>%
      dplyr::slice_tail(n=10) %>%
      tidyr::drop_na(TCbest)

    if (nrow(meanC_10y_df)<10){
      meanC_10y <- NA
    } else {
      meanC_10y <- meanC_10y_df %>%
        dplyr::pull(TCbest) %>%
        mean()
    }

    mature_g <- mature_ts %>%
      dplyr::filter(year>=min(bef_assess)) %>%
      dplyr::slice_head(n=1) %>%
      dplyr::mutate(mature_g = ifelse(TCbest < 1.25*meanC_10y, TRUE, FALSE)) %>%
      dplyr::pull(mature_g)

  } else { mature_g <- NA}

  if(length(bef_assess)>0){

    mature_ts <- mature_ts %>%
      dplyr::left_join(
        mature_ts %>%
          dplyr::filter(year>=min(bef_assess)) %>%
          dplyr::select(year) %>%
          dplyr::mutate(mature_g = mature_g),
        by="year")

  } else {
    mature_ts$mature_g <- mature_g
  }

  mature_ts <- mature_ts %>%
    # truncated related to the onset of the biomass timeseries:
    dplyr::select(stockid, year, dplyr::contains("mature_"), TCbest) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(mature = any(mature_a, mature_b, mature_c, mature_d,
                               na.rm=TRUE) & (!is.na(mature_c) | !is.na(TCbest)),
                  # to exclude years without catch data
                  mature_trunc = any(mature_f, mature_g, na.rm=TRUE)) %>%
    dplyr::ungroup()

  if(length(bef_assess)>0){

  mature_trunc <- mature_ts %>%
    dplyr::filter(year>=min(bef_assess)) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::pull(mature_trunc)

  } else {
    mature_trunc <- NA
  }
  # or do not integrate truncation into the data frame

  # Get onset of mature phase if not from the beginning:
  if (is.na(mature_trunc) | mature_trunc==FALSE){

    start_mature_df <- mature_ts %>%
      dplyr::filter(mature == TRUE) %>%
      tidyr::drop_na(mature_a, mature_b) # To begin in years with B and U

    if (nrow(start_mature_df)==0){

      start_mature <- NA

    } else {

      start_mature <- start_mature_df %>%
        dplyr::pull(year) %>% min()
    }

  } else if (mature_trunc==TRUE){

    start_mature <- mature_ts %>%
      dplyr::filter(mature_trunc == TRUE & (!is.na(mature_c) | !is.na(TC))) %>%
      dplyr::slice_head(n=1) %>%
      dplyr::pull(year)
  }

  return(tibble::tibble("stockid"=id,
                        "start_mature"=start_mature,
                        "mature_trunc"=mature_trunc))
}

#' Define mature phase of fishery
#'
#' @param id Short stock ID character.
#'
#' @return Data frame with the first year of the mature phase, whether the
#' onset of the mature phase occurred before data are available (truncated),
#' or NA if information to define mature is missing.
#'
#' @export

mature_def_simpl <- function(id){

  # Criteria for mature fishery from Melnychuk et al. 2021
  mature_ts <- timeseries_values_views %>%
    dplyr::filter(stockid==id) %>%
    dplyr::mutate(mature_a = ifelse(BdivBmsypref<0.8 | BdivBmgtpref<0.8,
                                    TRUE, FALSE),
                  mature_b = ifelse(UdivUmsypref>1 | UdivUmgtpref>1,
                                    TRUE, FALSE),
                  mature_c = ifelse(CdivMSY>1, TRUE, FALSE),
                  mature_d = ifelse(CdivMEANC>1.25, TRUE, FALSE),
                  # criterion (e) on rebuilding plans (not used)
                  mature_f = ifelse(BdivBmsypref<1 | BdivBmgtpref<1,
                                    TRUE, FALSE)
    )

  # Preliminary computation for criterion (g):
  # Get first year with both B/Bref and U/Uref:
  bef_assess <- mature_ts %>%
    dplyr::filter((dplyr::if_any(contains("BdivB"),
                                 function(x) !is.na(x)) &
                     dplyr::if_any(contains("UdivU"),
                                   function(x) !is.na(x)))) %>%
    dplyr::pull(year)

  if(length(bef_assess)>0){
    # Catch from the 10y before 1st B/Bref & U/Uref year:
    meanC_10y_df <- mature_ts %>%
      dplyr::filter(year < min(bef_assess)) %>%
      dplyr::slice_tail(n=10) %>%
      tidyr::drop_na(TCbest)

    # Discard if not at least 10 years:
    if (nrow(meanC_10y_df)<10){
      meanC_10y <- NA

    # Compute mean if 10 years complete:
    } else {
      meanC_10y <- meanC_10y_df %>%
        dplyr::pull(TCbest) %>%
        mean()
    }

    # Compute criterion (g)
    mature_g <- mature_ts %>%
      dplyr::filter(year>=min(bef_assess)) %>%
      dplyr::slice_head(n=1) %>%
      dplyr::mutate(mature_g = ifelse(TCbest < 1.25*meanC_10y, TRUE, FALSE)) %>%
      dplyr::pull(mature_g)

  } else { mature_g <- NA}

  if(length(bef_assess)>0){

    mature_ts <- mature_ts %>%
      dplyr::left_join(
        mature_ts %>%
          dplyr::filter(year>=min(bef_assess)) %>%
          dplyr::select(year) %>%
          dplyr::mutate(mature_g = mature_g),
        by="year")

  } else {
    mature_ts$mature_g <- mature_g
  }

  mature_ts <- mature_ts %>%
    # truncated related to the onset of the biomass timeseries:
    dplyr::select(stockid, year, dplyr::contains("mature_"), TCbest) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(mature = any(mature_a, mature_b, mature_c, mature_d,
                               na.rm=TRUE) & (!is.na(mature_c) | TCbest>0),
                  # to exclude years without catch data
                  mature_trunc = any(mature_f, mature_g, na.rm=TRUE)) %>% ## HERE
    dplyr::ungroup()

  if(length(bef_assess)>0){

    mature_trunc <- mature_ts %>%
      dplyr::filter(year>=min(bef_assess)) %>%
      dplyr::slice_head(n=1) %>%
      dplyr::pull(mature_trunc)

  } else {
    mature_trunc <- NA
  }
  # or do not integrate truncation into the data frame

  # Get onset of mature phase if not from the beginning:
  if (is.na(mature_trunc) | mature_trunc==FALSE){

    start_mature_df <- mature_ts %>%
      dplyr::filter(mature == TRUE)
      # tidyr::drop_na(mature_a, mature_b) # To begin in years with B and U

    if (nrow(start_mature_df)==0){

      start_mature <- NA

    } else {

      start_mature <- start_mature_df %>%
        dplyr::pull(year) %>% min()
    }

  } else if (mature_trunc==TRUE){

    start_mature <- mature_ts %>%
      dplyr::filter(mature_trunc == TRUE & (!is.na(mature_c) | TCbest>0)) %>%
      dplyr::slice_head(n=1) %>%
      dplyr::pull(year)
  }

  return(tibble::tibble("stockid"=id,
                        "start_mature"=start_mature,
                        "mature_trunc"=mature_trunc))
}


#' Remove model output when catch not available
#'
#' @param
#'
#' @return Data frame with all timeseries (timeseries_values_views) truncated
#' when needed.
#'
#' @export

relevant_ts <- function(){


  timeseries_values_views_bis <- timeseries_values_views %>%
    dplyr::filter(!is.na(TC)|!is.na(CdivMSY)|!is.na(TCbest)|!is.na(TL))

  nrow(timeseries_values_views_bis)

  t <- timeseries_values_views_bis %>%
    dplyr::group_by(stockid) %>%
    dplyr::slice_tail(n=1) %>%
    pull(stockid)

  length(t)
  nrow(timeseries_values_views_bis)
    # pull(year) %>%
    # mean()

  if(is.null(timeseries_values_views$SProd)){

    timeseries_values_views <<- timeseries_values_views %>%
      dplyr::group_by(stockid) %>%
      dplyr::mutate(SProd = dplyr::lead(TBbest)-TBbest+TCbest) %>%
      dplyr::ungroup()


    timeseries <<- timeseries %>%
      dplyr::bind_rows(
        timeseries_values_views %>%
          dplyr::select(c(stockid, stocklong, year, SProd)) %>%
          dplyr::rename(tsyear = year) %>%
          tidyr::pivot_longer(cols=SProd,
                              names_to = "tsid",
                              values_to = "tsvalue",
                              values_drop_na = TRUE) %>%
          dplyr::mutate(tsid = "SProd-MT") %>%
          # Add best assessment for TB:
          dplyr::left_join(timeseries_assessments_views %>%
                             dplyr::select(TBbest, stockid) %>%
                             dplyr::rename(assessid = TBbest),
                           by = "stockid")
      )

    SProd_years <<- timeseries_values_views %>%
      dplyr::select(stockid, stocklong, year, SProd) %>%
      dplyr::group_by(stockid) %>%
      tidyr::drop_na(SProd) %>%
      dplyr::summarise(SProd = paste0(min(year, na.rm=FALSE),"-",
                                      max(year, na.rm=FALSE)))

    timeseries_years_views <<- timeseries_years_views %>%
      dplyr::left_join(SProd_years, by="stockid")

  }

}





# II) Timeseries exploration ----------------------------------------------

#' List the types of time series available for a given stock
#'
#' @param id Short stock ID character.
#'
#' @return List of types of time series available for this stock.
#' @export

ts_available <- function(id){

  tsmetr_light <- tsmetrics %>%
    dplyr::arrange(tsshort) %>%
    tidyr::separate(tsshort, c("tsshort.1","tsshort.2"),"-") %>%
    dplyr::select(tsshort.1, tslong) %>%
    dplyr::rename(ts_type = tsshort.1) %>%
    tidyr::separate(tslong, c("tslong.1","tslong.2"),"[(]") %>%
    dplyr::select(ts_type, tslong.1) %>%
    dplyr::rename(tslong = tslong.1) %>%
    dplyr::group_by(ts_type) %>%
    dplyr::slice(1) %>%
    dplyr::group_by(tslong = stringi::stri_trans_totitle(tslong)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    suppressWarnings() %>%
    dplyr::bind_rows(c(ts_type="SProd", tslong="Surplus Production"))

  types <- timeseries_years_views %>%
    dplyr::filter(stockid == id)

  ts_type <- names(types[,!is.na(types)])[-c(1,2)]
  ts_type_df <- data.frame(ts_type) %>%
    dplyr::left_join(tsmetr_light, by="ts_type") %>%
    dplyr::rename(meaning = tslong)

  return(ts_type_df)
}



#' Perform breakpoint analysis
#'
#' @param id Short stock ID character.
#' @param ts_type Type of time series (TBbest, ERbest...) character.
#'
#' @return Output from breakpoint analysis.
#' @export

shifts_RAM <- function(id, ts_type){

  ts <- extract_RAM(id, ts_type)

  incomplete <- FALSE
  if(nrow(ts) != tail(ts$year,1) - head(ts$year,1) + 1){
    warning("Time series might be incomplete")
    incomplete <- TRUE
  }

  brk_out <- shifts(ts, abr_mtd=c("asd","chg"), asd_thr=0.15, asd_chk=TRUE,
                    lowwl=5, highwl="default", mad_thr=3, mad_cst=1.4826)

  return(brk_out)
}




#' Extract a stock RAMLBD time series and plot it
#'
#' @param id Short stock ID character.
#' @param ts_type Type of time series (TBbest, ERbest...) character.
#'
#' @return List of three objects:
#' - Timeseries plot
#' - Timeseries itself
#' @export

plot_RAM_bp <- function(id, ts_type){

  ts <- extract_RAM(id, ts_type)

  stock_ts <- ts %>%
    pull(ts_type) %>%
    ts(start=ts$year[1], end=tail(ts$year,1))

  metadt <- metadata %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(c(2,3,5,6,10))

  taxo <- taxonomy %>%
    dplyr::filter(scientificname == metadt[,3]) %>%
    dplyr::select(c(2,4,6,7,15))

  ts_info <- metadt %>%
    dplyr::left_join(taxo, by = "scientificname")

  stocklong <- metadt$stocklong
  sciname <- ts_info$scientificname
  famname <- ts_info$family
  reg <- ts_info$region
  ts <- ts %>% dplyr::rename(iter=scen)

  p <- plot_simu_simple(lib_ts=ts, var=ts_type)

  title <- paste0(stocklong, " [", id, "] – (",
                  ts$year[1], "–", tail(ts$year, 1), ")\n",
                  sciname, " – ", famname, " – ", reg, "\n",
                  ts_type)

  pop <- p + ggtitle(title)

  return(list("pop"=pop, "ts"=ts))
}






#' Plot a given stock time series searched first by region
#'
#' Prints a time series plot and the data frame used to plot it, and
#' upon request leading indicators plots and the data frame used to plot them.
#'
#' @param
#'
#' @return No return value.
#' @export

search_region <- function(){

  reg_list <- stock %>% as.data.frame() %>%
    dplyr::pull(region) %>% unique() %>% sort()
  print(reg_list)
  reg_nb <- readline(prompt="Choose region line: ") %>% as.integer()
  reg <- reg_list[reg_nb]

  taxa_list <- stock %>% as.data.frame() %>% dplyr::filter(region == reg)
  print(taxa_list[,c("commonname","scientificname","stocklong","stockid")])
  print(paste("Region: ", reg))
  id_nb <- readline(prompt="Choose taxa line: ") %>% as.integer()
  id <- taxa_list[id_nb,"stockid"]

  ts_type_list <- ts_available(id)
  rep <- TRUE
  while(rep){
    print(ts_type_list)
    print(paste("Region:", reg))
    print(paste("Stock:", taxa_list[id_nb,"stocklong"]))
    ts_type_nb <- readline(prompt="Choose time series type: ") %>% as.integer()
    ts_type <- ts_type_list[ts_type_nb, "ts_type"]
    dev.new()
    print(plot_RAM_bp(id,ts_type)$pop)
    print(paste("stockid:",id,"–",ts_type))
    rep <- ifelse(readline(prompt="Want to see another time series?
                           \n 1: Yes\n 0: No ")=="1", TRUE, FALSE)
  }

  return(invisible(NULL))
}


#' Plot a given stock time series searched first by taxa
#'
#' Prints a time series plot and the data frame used to plot it, and
#' upon request leading indicators plots and the data frame used to plot them.
#'
#' @param
#'
#' @return No return value.
#' @export

search_taxa <- function(){

  taxa_list <- taxonomy %>%
    dplyr::select(classname, family, scientificname, commonname1) %>%
    dplyr::arrange(classname, family, scientificname)
  print(unique(taxa_list$classname))
  class_nb <- readline(prompt="Choose class line: ") %>% as.integer()
  class <- unique(taxa_list$classname)[class_nb]

  fams <- taxa_list %>%
    dplyr::filter(classname == class) %>%
    dplyr::select(-classname)
  print(fams)
  print(paste("Class: ", class))
  sp_nb <- readline(prompt="Choose taxa line: ") %>% as.integer()
  sp <- fams[sp_nb,"scientificname"]

  reg_list <- stock %>%
    as.data.frame() %>%
    dplyr::filter(scientificname==sp) %>%
    dplyr::select(region, stocklong, stockid) %>%
    dplyr::arrange(region)
  print(reg_list)
  print(paste("Species: ", sp))
  id_nb <- readline(prompt="Choose stock line: ") %>% as.integer()
  id <- reg_list[id_nb,"stockid"]

  ts_type_list <- ts_available(id)
  rep <- TRUE
  while(rep){
    print(ts_type_list)
    print(paste("Species: ", sp))
    print(paste("Region:", reg_list[id_nb,"stocklong"]))
    ts_type_nb <- readline(prompt="Choose time series type: ") %>% as.integer()
    ts_type <- ts_type_list[ts_type_nb,"ts_type"]
    # dev.new()
    print(plot_RAM_bp(id,ts_type)$pop)
    print(paste("stockid:",id,"–",ts_type))

    rep <- ifelse(readline(prompt="Want to see another time series?
                           \n 1: Yes\n 0: No ")=="1", TRUE, FALSE)
  }

  return(invisible(NULL))
}



# III) Classification ------------------------------------------------------


#' Run time series classification (for empirical data)
#'
#' @param df_list List of data frames with timeseries to classify.
#' @param min_len Minimal timeseries length to consider.
#' @param group Column name containing the timeseries identification name (e.g.
#' species or population name).
#' @param time Column name containing timeseries time units.
#' @param variable Column name containing timeseries values.
#' @param str Character specifying how the best trajectory is chosen, based on
#' lowest AICc only ("aic") or also the confirmation with asdetect ("aic_asd").
#' @param run_loo Logical, whether to perform leave-one-out process.
#' @param two_bkps Logical, if true, looks for more than one breakpoints.
#' @param makeplots Logical that indicate whether to generate plots
#' (slows down computation if true).
#' @param ind_plot Character to display the plot for a given trajectory
#' (either "nch", "lin", "pol", "abr", or "best") or all four if NULL.
#' @param dirname Directory name where to save plots.
#' @param save_plot Logical, if plots are made, whether to save plots
#' (either the four trajectory fits or one of a given fit).
#' @param ... Additional arguments to be passed.
#'
#' @return List of three objects:
#' - Data frame with classification output with one row for each timeseries
#' - List of detailed classification output for each timeseries
#' - List summarizing the parameters used for the classification
#'
#' @export

run_classif_data <- function(df_list, min_len=20, group, time, variable,
                             str, run_loo, two_bkps,
                             makeplots=FALSE, ind_plot=NULL,
                             dirname=NULL, save_plot=TRUE, ...){
  # Load arguments:
  l <- list(...)

  # List of timeseries meeting the length criterion:
  df_list <- df_list[df_list %>% lapply(nrow)>=min_len]

  # Classify for all timeseries of this type:
  traj_ts <- data.frame()
  outlist <- list()

  for (i in 1:length(df_list)){ # by group or population

    set <- df_list[[i]] %>%
      dplyr::rename(scen = dplyr::all_of(group)) %>%
      dplyr::select(scen, dplyr::all_of(time), all_of(variable)) %>%
      tidyr::drop_na() %>%
      prep_data_simpl()

    set_length <- nrow(set$ts[[1]])

    if(str == "aic") abr_mtd <- c("chg")
    if(str == "aic_asd") abr_mtd <- c("chg", "asd")

    trajs <- traj_class(sets=set, str=str, abr_mtd=abr_mtd,
                        run_loo=run_loo, two_bkps=two_bkps,
                        makeplots=makeplots, ind_plot=ind_plot,
                        dirname=dirname, save_plot=save_plot, ...)

    traj_ts <- traj_ts %>% dplyr::bind_rows(trajs$best_traj)
    outlist[[names(df_list)[i]]] <- trajs

    if (i%%10 == 0) print(paste0(i,"/",length(df_list)))
  }

  traj_ts_full <- traj_ts %>%
    dplyr::mutate(species = simu_id %>%
                    sub("_iter01","", .))

  return(list("traj_ts_full"=traj_ts_full, "outlist"=outlist,
              "param_list"=trajs$param_list))

}


#' Run time series classification on multiple time series types (for empirical data)
#'
#' @param df_list List of data frames with timeseries to classify.
#' @param min_len Minimal timeseries length to consider.
#' @param group Column name containing the timeseries identification name (e.g.
#' species or population name).
#' @param time Column name containing timeseries time units.
#' @param variable Column name containing timeseries values.
#' @param str Character specifying how the best trajectory is chosen, based on
#' lowest AICc only ("aic") or also the confirmation with asdetect ("aic_asd").
#' @param run_loo Logical, whether to perform leave-one-out process.
#' @param two_bkps Logical, if true, looks for more than one breakpoints.
#' @param makeplots Logical that indicate whether to generate plots
#' (slows down computation if true).
#' @param ind_plot Character to display the plot for a given trajectory
#' (either "nch", "lin", "pol", "abr", or "best") or all four if NULL.
#' @param detection_plot Logical to plot detection plot when abrupt.
#' @param dirname Directory name where to save plots.
#' @param save_plot Logical, if plots are made, whether to save plots
#' (either the four trajectory fits or one of a given fit).
#' @param save_detplot Logical, if plots are made, whether to save plot_bis
#' (best fit with asdetect detection curve).
#' @param showlastplot Logical, last plot in return value.
#' @param ... Additional arguments to be passed.
#'
#' @return List of two objects:
#' - Data frame with classification output with one row for each timeseries
#' - List of detailed classification output for each timeseries
#' - List summarizing the parameters used for the classification
#'
#' @export

run_classif_data_multi <- function(df_list, min_len=20, group, time, variable,
                                   str, run_loo, two_bkps,
                                   makeplots=FALSE, ind_plot=NULL,
                                   dirname=NULL, save_plot=TRUE,
                                   showlastplot=FALSE, ...){

  # Load arguments:
  l <- list(...)

  # List of timeseries meeting the length criterion:
  df_list <- df_list[df_list %>% lapply(nrow)>=min_len]

  # Classify for all timeseries of this type:
  traj_ts <- data.frame()
  outlist <- list()

  for (i in 1:length(df_list)){ # by group or population

    fig_list <- list() # List of figures with different time series types

    for (j in 1:length(variable)){ # by time series type

      set <- df_list[[i]] %>%
        dplyr::rename(scen = dplyr::all_of(group),
                      year = dplyr::all_of(time)) %>%
        dplyr::select(scen, year, all_of(variable[j])) %>%
        tidyr::drop_na() %>%
        prep_data(thr=NULL, type="data", apriori=FALSE)

      set_length <- nrow(set$ts[[1]])

      if(str == "aic") abr_mtd <- c("chg")
      if(str == "aic_asd") abr_mtd <- c("chg", "asd")

      trajs <- traj_class(sets=set, str=str, abr_mtd=abr_mtd,
                          run_loo=run_loo, two_bkps=two_bkps,
                          makeplots=makeplots, ind_plot=ind_plot,
                          dirname=dirname, save_plot=save_plot, ...)

      traj_ts <- traj_ts %>%
        dplyr::bind_rows(trajs$best_traj %>%
                           dplyr::mutate(ts_types = variable[j]))
      outlist[[names(df_list)[i]]] <- trajs

      fig_list[[paste0(names(df_list)[i],"_",variable[j])]] <- trajs$class_plot$plots[[1]]

    }

    if (showlastplot==TRUE){

      min_x <- lapply(fig_list, function(x) x$data$X %>% min()) %>%
        unlist() %>% min()
      max_x <- lapply(fig_list, function(x) x$data$X %>% max()) %>%
        unlist() %>% max()

      fig_list <- lapply(fig_list, function(x) x+
                           scale_x_continuous(limits = c(min_x, max_x)))

      all_plots <- cowplot::plot_grid(plotlist = fig_list, ncol = 1)

    }

    # png(filename = paste0("analyses/class_plot/", names(df_list)[i], "_",
    #                       paste(variable, collapse="."),"_",
    #                       str,"_coeff",smooth_signif,"_loo",run_loo,".png"),
    #     width = 12, height=25, units = "cm", res=300)
    # print(all_plots)
    # dev.off()

    if (i%%10 == 0) print(paste0(i,"/",length(df_list)))
  }

  # Add start, end and length:
  traj_ts <- traj_ts %>%
    dplyr::rename(!!group := simu_id) %>%
    dplyr::left_join(lapply(1:length(df_list), function(i){
      data.frame(group = names(df_list)[i],
                 first = df_list[[i]] %>% pull(dplyr::all_of(time)) %>% min(),
                 last = df_list[[i]] %>% pull(dplyr::all_of(time)) %>% max(),
                 length = df_list[[i]] %>% nrow())
    }) %>% do.call(what=dplyr::bind_rows) %>%
      dplyr::rename(!!group := group),
    by=group)

  return(list("traj_ts_full"=traj_ts, "outlist"=outlist,
              "param_list"=trajs$param_list))

}


#' Add collapse to classification summary according to a given definition
#'
#' @param traj_SProd A classification summary.
#' @param coll_def A vector defining the threshold of collapsed state.
#' The 1st value indicates the parameter considered ("bavg", "bmax", or "bmsy",
#' for average, maximum biomass, and biomass at MSY respectively). The second
#' value indicates the fraction of biomass parameter considered as threshold.
#' @param csec_coll An integer specifying the minimal number of consecutive
#' years of under collapse threshold to be considered as collapsed.
#'
#' @return The classification summary with information about collapse.
#' @export
#'
#' @examples

add_collapsed <- function(traj_SProd, coll_def=c("bavg", 0.25), csec_coll=5){

  # Filter biomass timeseries:
  tb <- timeseries_values_views %>%
    dplyr::filter(stockid %in% (traj_SProd %>% dplyr::pull(stockid))) %>%
    dplyr::select(stockid, year, TBbest) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(traj_SProd %>%
                       dplyr::mutate(length = last-first+1) %>%
                       dplyr::select(stockid, first, last, length),
                     by="stockid") %>%
    dplyr::group_by(stockid) %>%
    # Average and maximum biomass estimated over all data available
    dplyr::mutate(
      tb_tbmean = TBbest/mean(TBbest),
      tb_tbmax = TBbest/max(TBbest),
      first_max = year[which(TBbest == max(TBbest))][1]
      # assuming that collapse should not precede biomass maximum
    ) %>%
    dplyr::filter(year>=first & year<=last)

  # Specify when stocks collapsed:
  if (coll_def[1]=="bavg"){
    tb <- tb %>%
      dplyr::mutate(collapsed = ifelse(tb_tbmean<coll_def[2] &
                                         year>first_max,
                                       TRUE, FALSE)) %>%
      dplyr::group_by(stockid) %>%
      # Add length of collapsed state:
      dplyr::mutate(run_len = rep(rle(collapsed)[1]$length,
                                  rle(collapsed)[1]$length),
                    run_len = ifelse(collapsed==TRUE, run_len, NA))

  } else if (coll_def[1]=="bmax"){
    tb <- tb %>%
      dplyr::mutate(collapsed = ifelse(tb_tbmax<coll_def[2] &
                                         year>first_max,
                                       TRUE, FALSE)) %>%
      dplyr::group_by(stockid) %>%
      # Add length of collapsed state:
      dplyr::mutate(run_len = rep(rle(collapsed)[1]$length,
                                  rle(collapsed)[1]$length),
                    run_len = ifelse(collapsed==TRUE, run_len, NA))

  } else if (coll_def[1]=="bmsy"){
    tb <- tb %>%
      dplyr::left_join(
        timeseries_values_views %>%
          dplyr::filter(stockid %in% (traj_SProd %>% dplyr::pull(stockid))) %>%
          dplyr::select(stockid, year, BdivBmsypref) %>%
          tidyr::drop_na() %>%
          dplyr::group_by() %>%
          dplyr::mutate(collapsed = ifelse(BdivBmsypref<coll_def[2],
                                           TRUE, FALSE)),
        by=c("stockid", "year")
      ) %>%
      dplyr::group_by(stockid) %>%
      # Add length of collapsed state:
      dplyr::mutate(run_len = rep(rle(collapsed)[1]$length,
                                  rle(collapsed)[1]$length),
                    run_len = ifelse(collapsed==TRUE, run_len, NA))
  }

  # Get collapse year:
  coll_year <- tb %>%
    dplyr::filter(collapsed==TRUE & run_len>=csec_coll) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::select(stockid, year) %>%
    dplyr::rename(first_coll_y = year)

  tb <- tb %>%
    dplyr::left_join(coll_year, by="stockid") %>%
    dplyr::select(stockid, year, tb_tbmax, collapsed, first_coll_y)

  traj_SProd_coll <- traj_SProd %>%
    dplyr::left_join(coll_year, by="stockid") %>%
    dplyr::mutate(did_collapse = ifelse(!is.na(first_coll_y), "Yes", "No"),
                  coll_def = paste(coll_def, collapse = ""),
                  traj = dplyr::case_when(
                    traj == "1_breakpoint" &
                      trend == "decrease" ~ "decrease_abrupt",
                    traj == "1_breakpoint" &
                      trend == "increase" ~ "increase_abrupt",
                    traj != "1_breakpoint" ~ traj),
                  traj = plyr::revalue(
                    traj, c("increase_constant"="increase_linear",
                            "decrease_constant"="decrease_linear")),
                  traj = sub("_", " ", traj),
                  traj = factor(
                    traj,
                    levels = c("stable constant",
                               "stable concave","stable convex",
                               "increase linear","increase decelerated",
                               "increase accelerated","increase abrupt",
                               "decrease linear","decrease decelerated",
                               "decrease accelerated","decrease abrupt")))

  return(traj_SProd_coll)
}



# IIIb) Trait retrieval ---------------------------------------------------

traits_fun <- function(traj_SProd){

  # Add names to match databases:
  name_table <- traj_SProd %>%
    dplyr::select(stockid, scientificname) %>%
    dplyr::mutate(

      # Names for FishBase population growth traits (tm, Lm):
      fishbase_name = scientificname,
      fishbase_name =
        plyr::revalue(fishbase_name,
                      c("Theragra chalcogramma" = "Gadus chalcogrammus",
                        "Zearaja chilensis" = "Dipturus chilensis",
                        "Tetrapturus albidus" = "Kajikia albida",
                        "Sardinops melanostictus" = "Sardinops sagax",
                        "Raja rhina" = "Beringraja rhina",
                        "Neoplatycephalus richardsoni" = "Platycephalus richardsoni",
                        "Limanda ferruginea" = "Myzopsetta ferruginea",
                        "Etrumeus teres" = "Etrumeus sadina",
                        "Epinephelus niveatus" = "Hyporthodus niveatus",
                        "Epinephelus flavolimbatus" = "Hyporthodus flavolimbatus",
                        "Dentex tumifrons" = "Evynnis tumifrons",
                        "Chrysophrys auratus" = "Pagrus auratus",
                        "Bathyraja parmifera" = "Arctoraja parmifera"
                      )),

      # Names for FishLife population growth traits (tm, Lm):
      fishlife_name = scientificname,
      fishlife_name =
        plyr::revalue(fishlife_name,
                      c(
                        "Ammodytes spp"= "Ammodytes",
                        "Merluccius spp" = "Merluccius",
                        "Sebastes spp" = "Sebastes",
                        "Tetrapturus albidus" = "Kajikia albida",
                        "Sardinops melanostictus" = "Sardinops sagax",
                        "Neoplatycephalus richardsoni" = "Platycephalus richardsoni",
                        "Etrumeus teres" = "Etrumeus sadina",
                        "Epinephelus niveatus" = "Hyporthodus niveatus",
                        "Epinephelus flavolimbatus" = "Hyporthodus flavolimbatus",
                        "Chrysophrys auratus" = "Pagrus auratus",
                        "Bathyraja parmifera" = "Bathyraja",
                        "Centroberyx gerrardi" = "Centroberyx",
                        "Clupea pallasii" = "Clupea pallasii pallasii",
                        "Hexagrammos decagrammus" = "Hexagrammos",
                        "Merluccius gayi" = "Merluccius gayi gayi",
                        "Platycephalus conatus" = "Platycephalus",
                        "Scomber colias" = "Scomber",
                        "Sebastes auriculatus" = "Sebastes",
                        "Sebastes aurora" = "Sebastes",
                        "Sebastes carnatus" = "Sebastes",
                        "Sebastes chlorostictus" = "Sebastes",
                        "Sebastes diploproa" = "Sebastes",
                        "Sebastes entomelas" = "Sebastes",
                        "Sebastes nebulosus" = "Sebastes",
                        "Sebastolobus altivelis" = "Sebastolobus",
                        "Strangomera bentincki" = "Clupea bentincki")))

  traj_SProd <- traj_SProd %>%
    dplyr::left_join(name_table, by=c("stockid", "scientificname"))


  ## FishBase:
  # Maximum weight, growth rate, lifespan...
  LH_grw <- rfishbase::popgrowth(
    traj_SProd %>% dplyr::pull(fishbase_name)) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(mean_Loo = mean(Loo, na.rm = TRUE),
                     mean_K = mean(K, na.rm = TRUE),
                     mean_Winfinity = mean(Winfinity, na.rm = TRUE),
                     mean_tmax = mean(tmax, na.rm = TRUE),
                     mean_tm = mean(tm, na.rm = TRUE),
                     mean_M = mean(M, na.rm = TRUE),
                     mean_Lm = mean(Lm, na.rm = TRUE))

  # Trophic level
  LH_troph <- rfishbase::ecology(traj_SProd %>%
                                   dplyr::pull(fishbase_name)) %>%
    dplyr::select(Species, FoodTroph) %>%
    # https://github.com/ropensci/rfishbase/issues/199
    # "FoodTroph gives a MonteCarlo estimate of trophic level based on known food items.
    # DietTroph uses the mean or median of trophic levels derived from actual diet composition studies.
    # While in theory troph from diet should be more reliable,
    # many diet studies are incomplete or biased and I often find FoodTroph more reasonable." Rainer Froese
    dplyr::group_by(Species) %>%
    dplyr::slice(1) # make sure to keep only one row by species

  # Type of habitat (demersal, pelagic, reef)
  LH_hab <- rfishbase::fb_tbl("species") %>%
    dplyr::mutate(Species = paste(Genus, Species)) %>%
    dplyr::filter(Species %in% (traj_SProd %>% dplyr::pull(fishbase_name))) %>%
    dplyr::select(Species, DemersPelag) %>%
    dplyr::bind_rows(
      data.frame(Species=c("Ammodytes spp", "Merluccius spp", "Sebastes spp"),
                 DemersPelag = c("demersal", "bathydemersal", "bathypelagic"))) %>%
    dplyr::mutate(DemersPelag = dplyr::case_when(
      DemersPelag %in% c("bathydemersal", "benthopelagic", "demersal", "reef-associated") ~ "demersal",
      DemersPelag %in% c("pelagic-neritic", "pelagic-oceanic", "bathypelagic") ~ "pelagic")
      # DemersPelag %in% c("bathydemersal", "bathypelagic", "demersal", "reef-associated") ~ "demersal",
      # DemersPelag %in% c("pelagic-neritic", "pelagic-oceanic", "benthopelagic") ~ "pelagic")
      # https://fishbase.se/manual/English/FishBaseThe_Species_Table.htm
    )

  # Habitat RAMLDB:
  hab_ram <- traj_SProd %>%
    dplyr::left_join(bioparams %>%
                       dplyr::filter(bioid=="Habitat-Habitat") %>%
                       dplyr::rename(habitat_raw=biovalue) %>%
                       dplyr::select(stockid, habitat_raw),
                     by="stockid") %>%
    dplyr::select(stockid, scientificname, habitat_raw) %>%
    dplyr::mutate(habitat_raw = stringr::str_to_lower(habitat_raw),
                  habitat = dplyr::case_when(
      grepl("pelagic|pekagic|pleagic|diadromous", habitat_raw) ~ "pelagic",
      grepl("demersal", habitat_raw) ~ "demersal",
      is.na(habitat_raw) ~ NA
    )) %>%
    dplyr::mutate(
      # Names for FishBase:
      scientificname =
        plyr::revalue(scientificname,
                      c(
                        "Theragra chalcogramma" = "Gadus chalcogrammus",
                        "Zearaja chilensis" = "Dipturus chilensis",
                        "Tetrapturus albidus" = "Kajikia albida",
                        "Sardinops melanostictus" = "Sardinops sagax",
                        "Raja rhina" = "Beringraja rhina",
                        "Neoplatycephalus richardsoni" = "Platycephalus richardsoni",
                        "Limanda ferruginea" = "Myzopsetta ferruginea",
                        "Etrumeus teres" = "Etrumeus sadina",
                        "Epinephelus niveatus" = "Hyporthodus niveatus",
                        "Epinephelus flavolimbatus" = "Hyporthodus flavolimbatus",
                        "Dentex tumifrons" = "Evynnis tumifrons",
                        "Chrysophrys auratus" = "Pagrus auratus",
                        "Bathyraja parmifera" = "Arctoraja parmifera"
                      ), warn_missing=FALSE)) %>%
    dplyr::group_by(scientificname, habitat) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::group_by(scientificname) %>%
    dplyr::mutate(dup = n()) %>%
    dplyr::left_join(LH_hab, by=c("scientificname"="Species")) %>%
    tidyr::drop_na(habitat) %>%
    dplyr::mutate(no_match = ifelse(habitat != DemersPelag, "mismatch", NA))


  LH <-
    traj_SProd %>%
    dplyr::select(stockid, fishbase_name, classname, ordername, family) %>%
    dplyr::group_by(fishbase_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-stockid) %>%
    dplyr::rename(Species=fishbase_name) %>%
    dplyr::left_join(LH_grw, by="Species") %>%
    dplyr::left_join(LH_troph, by="Species") %>%
    dplyr::left_join(LH_hab, by="Species") %>%
    dplyr::select(Species, classname, ordername, family, FoodTroph, DemersPelag)

  ## FishLife:
  traits_fishlife <- FishLife::FishBase$beta_gv %>%
    tibble::as_tibble(rownames="taxa") %>%
    dplyr::mutate(taxa_type = dplyr::case_when(
      stringr::str_count(taxa, "_")==4 ~"species",
      stringr::str_count(taxa, "_")==3 ~"genus",
      stringr::str_count(taxa, "_")==2 ~"family",
      stringr::str_count(taxa, "_")==1 ~"order",
      stringr::str_count(taxa, "_")==0 ~"class")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Species = dplyr::case_when(
      taxa_type=="species" ~paste(strsplit(taxa, "_")[[1]][4:5], collapse=" "),
      taxa_type=="genus" ~strsplit(taxa, "_")[[1]][4],
      taxa_type=="family" ~strsplit(taxa, "_")[[1]][3],
      taxa_type=="order" ~strsplit(taxa, "_")[[1]][2],
      taxa_type=="class" ~strsplit(taxa, "_")[[1]][1]))

  fishlife_traits <- traj_SProd %>%
    dplyr::select(stockid, fishlife_name, phylum, classname, ordername, family) %>%
    dplyr::group_by(fishlife_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-stockid) %>%
    dplyr::rename(Species=fishlife_name) %>%
    dplyr::mutate(Species = sub(" spp", "", Species)) %>%
    dplyr::left_join(traits_fishlife %>%
                       dplyr::select(-taxa), by="Species")

  LH_fishlife_complete <-
    traj_SProd %>%
    dplyr::select(stockid, scientificname, fishbase_name, fishlife_name) %>%
    dplyr::group_by(scientificname) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-stockid) %>%
    dplyr::rename(Species=scientificname) %>%
    dplyr::left_join(fishlife_traits %>%
                       dplyr::select(-c(phylum, classname, ordername, family)),
                     by=c("fishlife_name"="Species")) %>%
    dplyr::left_join(LH, by=c("fishbase_name"="Species"))


  ## Fishing-related traits:
  fish_Umsy <-
    extract_RAM(traj_SProd %>% dplyr::pull(stockid),
                ts_type="UdivUmsypref", drop_na=TRUE) %>%
    dplyr::mutate(stockid = sub("UdivUmsypref_","",scen)) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, loc_brk_chg) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, loc_brk_chg),
      by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(
      max_U = max(UdivUmsypref),
      mean_U = mean(UdivUmsypref),
      U_change =
        summary(lm(UdivUmsypref~year))$coefficients["year","Estimate"],
      sd_U = sd(UdivUmsypref),
      detrsd_U = sd(pracma::detrend(UdivUmsypref)),
      len_U = length(UdivUmsypref),
      D_U = (1/(len_U-1))*sum(abs(log((dplyr::lead(UdivUmsypref)+0.01*mean_U)/
                                        (UdivUmsypref+0.01*mean_U))), na.rm=TRUE)
    ) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>%
    dplyr::mutate(
      max_U_prshf = max(UdivUmsypref),
      mean_U_prshf = mean(UdivUmsypref),
      U_change_prshf =
        summary(lm(UdivUmsypref~year))$coefficients["year","Estimate"],
      sd_U_prshf = sd(UdivUmsypref),
      detrsd_U_prshf = sd(pracma::detrend(UdivUmsypref)),
      len_U_prshf = length(UdivUmsypref),
      D_U_prshf = (1/(len_U-1))*sum(abs(log((dplyr::lead(UdivUmsypref)+0.01*mean_U)/
                                              (UdivUmsypref+0.01*mean_U))), na.rm=TRUE)
    ) %>%
    dplyr::slice(1) %>%
    tidyr::drop_na(mean_U) %>%
    dplyr::select(stockid,
                  max_U, mean_U, U_change,
                  sd_U, detrsd_U, D_U, len_U,
                  max_U_prshf, mean_U_prshf, U_change_prshf,
                  sd_U_prshf, detrsd_U_prshf, D_U_prshf, len_U_prshf
    )
  # There are two stocks ("COD2J3KL" and "RKCRABBB") for which the shift occurred
  # before values of UdivUmsypref were available.

  fish_ER <- traj_SProd %>%
    dplyr::pull(stockid) %>%
    extract_RAM(ts_type="ERbest2", drop_na=TRUE) %>%
    dplyr::mutate(stockid = sub("ERbest2_","",scen)) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, loc_brk_chg, first, last) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, loc_brk_chg, first, last),
      by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(mean_ER = mean(ERbest2),
                  len_ER = length(ERbest2)) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>%
    dplyr::filter(year>=first | year<=last) %>%
    dplyr::mutate(mean_ER_prshf = mean(ERbest2),
                  ER_change_prshf =
                    summary(lm(ERbest2~year))$coefficients["year","Estimate"],
                  len_ER_prshf = length(ERbest2)) %>%
    dplyr::slice(1) %>%
    tidyr::drop_na(mean_ER_prshf) %>%
    dplyr::select(stockid,
                  mean_ER, mean_ER_prshf, ER_change_prshf, len_ER, len_ER_prshf)

  # Minimum & Average B/Bmsy:
  fish_Bmsy <-
    extract_RAM(traj_SProd %>% dplyr::pull(stockid),
                ts_type="BdivBmsypref", drop_na=TRUE) %>%
    dplyr::mutate(stockid = sub("BdivBmsypref_","",scen)) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, loc_brk_chg, first, last) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, loc_brk_chg, first, last),
      by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(
      min_B = min(BdivBmsypref),
      mean_B = mean(BdivBmsypref),
      sd_B = sd(BdivBmsypref),
      detrsd_B = sd(pracma::detrend(BdivBmsypref)),
      len_B = length(BdivBmsypref),
      D_B = (1/(len_B-1))*sum(abs(log((dplyr::lead(BdivBmsypref)+0.01*mean_B)/
                                        (BdivBmsypref+0.01*mean_B))), na.rm=TRUE)
    ) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>%
    dplyr::mutate(
      min_B_prshf = min(BdivBmsypref),
      mean_B_prshf = mean(BdivBmsypref),
      sd_B_prshf = sd(BdivBmsypref),
      detrsd_B_prshf = sd(pracma::detrend(BdivBmsypref)),
      len_B_prshf = length(BdivBmsypref),
      D_B_prshf = (1/(len_B-1))*sum(log((dplyr::lead(BdivBmsypref)+0.01*mean_B)/
                                          (BdivBmsypref+0.01*mean_B)), na.rm=TRUE)
    ) %>%
    dplyr::slice(1) %>%
    tidyr::drop_na(mean_B) %>%
    dplyr::select(stockid,
                  min_B, mean_B,
                  sd_B, detrsd_B, D_B, len_B,
                  min_B_prshf, mean_B_prshf,
                  sd_B_prshf, detrsd_B_prshf, D_B_prshf, len_B_prshf
    )

  # Total exploitation time:
  fish_expl_time <- timeseries %>%
    dplyr::group_by(stockid) %>%
    tidyr::drop_na(tsvalue) %>%
    dplyr::mutate(first_year_avail = min(tsyear)) %>%
    dplyr::filter(tsyear==first_year_avail) %>%
    dplyr::slice(1) %>%
    dplyr::select(stockid, first_year_avail) %>%
    dplyr::filter(stockid %in% (traj_SProd %>%
                                  dplyr::pull(stockid))) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, first, last, loc_brk_chg) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, first, last, loc_brk_chg),
      by="stockid") %>%
    dplyr::mutate(total_time = last-first_year_avail+1,
                  time_prshf = ifelse(!is.na(loc_brk_chg),
                                      loc_brk_chg-first_year_avail+1,
                                      total_time))

  # Absolute catch (mean or cumulative):
  fish_catch <- extract_RAM(traj_SProd %>% dplyr::pull(stockid),
                            ts_type="TCbest", drop_na=TRUE) %>%
    dplyr::mutate(stockid = sub("TCbest_","",scen)) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, loc_brk_chg) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, loc_brk_chg),
      by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(
      sum_C = sum(TCbest, na.rm=TRUE),
      mean_C = mean(TCbest)
    ) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>%
    dplyr::mutate(
      sum_C_prshf = sum(TCbest, na.rm=TRUE),
      mean_C_prshf = mean(TCbest)
    ) %>%
    dplyr::slice(1) %>%
    tidyr::drop_na(mean_C) %>%
    dplyr::select(stockid,
                  sum_C, sum_C_prshf,
                  mean_C, mean_C_prshf)

  # Management approach
  fish_param <- traj_SProd %>%
    dplyr::select(stockid, scientificname) %>%
    dplyr::left_join(fish_expl_time, by="stockid") %>%
    dplyr::left_join(fish_Bmsy, by="stockid") %>%
    dplyr::left_join(fish_Umsy, by="stockid") %>%
    dplyr::left_join(fish_ER, by="stockid") %>%
    dplyr::left_join(fish_catch, by="stockid")


  # Environmental data:
  SST_had <-
    readr::read_csv("data/spatial_data/yearly_had_sst.csv") %>%
    dplyr::filter(stockid %in% (traj_SProd %>% dplyr::pull(stockid)))

  # Trim SST timeseries to match period available for productivity
  SST_had_SProd <- SST_had %>%
    dplyr::left_join(traj_SProd %>%
                       dplyr::select(stockid, first, last, loc_brk_chg, class), by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::filter(year>=first & year<=last) %>%  # to match period available
    tidyr::drop_na(sst_c) %>%
    dplyr::mutate(sst_avg = mean(sst_c),
                  sst_change =
                    summary(lm(sst_c~year))$coefficients["year","Estimate"]) %>%
    dplyr::mutate(loc_brk_chg = ifelse(class=="abrupt", loc_brk_chg, NA)) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) # to match preshift period

  sst_stock_list <- SST_had_SProd %>%
    tidyr::drop_na(sst_c) %>%
    dplyr::distinct(stockid) %>%
    dplyr::pull(stockid) %>%
    sort()

  SST_had_SProd_df <- SST_had_SProd %>%
    tidyr::drop_na(sst_c) %>%
    dplyr::mutate(sst_avg_prshf = mean(sst_c),
                  sst_cv = sd(sst_c)/mean(sst_c),
                  sst_sd = sd(sst_c),
                  sst_detrsd = sd(pracma::detrend(sst_c)),
                  sst_len = length(sst_c),
                  sst_change_prshf =
                    summary(lm(sst_c~year))$coefficients["year","Estimate"]) %>%
    dplyr::slice(1) %>%
    dplyr::select(stockid, sst_avg, sst_avg_prshf, sst_change, sst_change_prshf)

  ## Polygons:
  polygons <- traj_SProd %>%
    dplyr::left_join(
      readxl::read_xlsx(
        paste0("data/spatial_data/stock_boundaries/",
             "ramldb_v3.8_stock_boundary_table_v2_formatted.xlsx")) %>%
        dplyr::mutate(batch="assess1") %>%
        dplyr::bind_rows(
          readr::read_csv(
            paste0("data/spatial_data/stock_boundaries/",
                   "additional_ramldb_stock_boundary_formatted.csv")) %>%
            dplyr::mutate(batch="assess2")) %>%
        dplyr::rename(assessid_polygon_avail = assessid) %>%
        dplyr::group_by(stockid) %>%
        dplyr::mutate(n=n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n==1 | batch=="assess2") %>%
        dplyr::select(stockid, assessid_polygon_avail, batch),
      by="stockid")

  stock_pol <- polygons %>%
    tidyr::drop_na(assessid_polygon_avail) %>%
    dplyr::select(stockid, assessid_polygon_avail)

  centroids_available <- mapply(function(x, y){

    batch <- polygons %>%
      dplyr::filter(assessid_polygon_avail==x) %>%
      dplyr::pull(batch)

    if(batch=="assess1"){
      path <- "ramldb_boundaries"
    } else {
      path <- "additional_boundaries"
    }

    sf::read_sf(
      paste0("data/spatial_data/stock_boundaries/",path,"/",x,".shp")) %>%
      # transform to a projection with projected coordinates (World Mercator):
      sf::st_transform(3395) %>%
      sf::st_centroid() %>%
      # transform back to WGS 84 (with geographic coordinates):
      sf::st_transform(4326) %>%
      sf::st_coordinates() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(stockid=y) %>%
      suppressWarnings()
    },
    stock_pol$assessid_polygon_avail, stock_pol$stockid, SIMPLIFY = FALSE) %>%
    do.call(what=dplyr::bind_rows)

  polygons <- polygons %>%
    dplyr::left_join(centroids_available, by="stockid")

  # Merge all traits
  traits <- LH_fishlife_complete %>%
    dplyr::select(-c(taxa_type, fishbase_name, fishlife_name)) %>%
    dplyr::left_join(polygons %>%
                       dplyr::select(scientificname, stockid, X, Y) %>%
                       dplyr::rename(Species=scientificname), by="Species") %>%
    dplyr::left_join(fish_param %>%
                       dplyr::select(stockid, mean_B_prshf, mean_U_prshf, U_change_prshf,
                                     mean_ER_prshf, ER_change_prshf, sum_C_prshf),
                     by="stockid") %>%
    dplyr::left_join(SST_had_SProd_df %>%
                       dplyr::select(stockid, sst_avg, sst_avg_prshf, sst_change, sst_change_prshf),
                     by="stockid") %>%
    dplyr::relocate(stockid)

  return(traits)
}


show_traits <- function(traj_SProd){

  # Add names to match databases:
  name_table <- traj_SProd %>%
    dplyr::select(stockid, scientificname) %>%
    dplyr::mutate(

      # Names for FishBase population growth traits (tm, Lm):
      fishbase_name = scientificname,
      fishbase_name =
        plyr::revalue(fishbase_name,
                      c("Theragra chalcogramma" = "Gadus chalcogrammus",
                        "Zearaja chilensis" = "Dipturus chilensis",
                        "Tetrapturus albidus" = "Kajikia albida",
                        "Sardinops melanostictus" = "Sardinops sagax",
                        "Raja rhina" = "Beringraja rhina",
                        "Neoplatycephalus richardsoni" = "Platycephalus richardsoni",
                        "Limanda ferruginea" = "Myzopsetta ferruginea",
                        "Etrumeus teres" = "Etrumeus sadina",
                        "Epinephelus niveatus" = "Hyporthodus niveatus",
                        "Epinephelus flavolimbatus" = "Hyporthodus flavolimbatus",
                        "Dentex tumifrons" = "Evynnis tumifrons",
                        "Chrysophrys auratus" = "Pagrus auratus",
                        "Bathyraja parmifera" = "Arctoraja parmifera"
                      )),

      # Names for FishBase population growth traits (tm, Lm):
      fishlife_name = scientificname,
      fishlife_name =
        plyr::revalue(fishlife_name,
                      c(
                        "Ammodytes spp"= "Ammodytes",
                        "Merluccius spp" = "Merluccius",
                        "Sebastes spp" = "Sebastes",
                        "Tetrapturus albidus" = "Kajikia albida",
                        "Sardinops melanostictus" = "Sardinops sagax",
                        "Neoplatycephalus richardsoni" = "Platycephalus richardsoni",
                        "Etrumeus teres" = "Etrumeus sadina",
                        "Epinephelus niveatus" = "Hyporthodus niveatus",
                        "Epinephelus flavolimbatus" = "Hyporthodus flavolimbatus",
                        "Chrysophrys auratus" = "Pagrus auratus",
                        "Bathyraja parmifera" = "Bathyraja",
                        "Centroberyx gerrardi" = "Centroberyx",
                        "Clupea pallasii" = "Clupea pallasii pallasii",
                        "Hexagrammos decagrammus" = "Hexagrammos",
                        "Merluccius gayi" = "Merluccius gayi gayi",
                        "Platycephalus conatus" = "Platycephalus",
                        "Scomber colias" = "Scomber",
                        "Sebastes auriculatus" = "Sebastes",
                        "Sebastes aurora" = "Sebastes",
                        "Sebastes carnatus" = "Sebastes",
                        "Sebastes chlorostictus" = "Sebastes",
                        "Sebastes diploproa" = "Sebastes",
                        "Sebastes entomelas" = "Sebastes",
                        "Sebastes nebulosus" = "Sebastes",
                        "Sebastolobus altivelis" = "Sebastolobus",
                        "Strangomera bentincki" = "Clupea bentincki")))

  traj_SProd <- traj_SProd %>%
    dplyr::left_join(name_table, by=c("stockid", "scientificname"))


  ## FishBase:
  # Maximum weight, growth rate, lifespan...
  LH_grw <- rfishbase::popgrowth(
    traj_SProd %>% dplyr::pull(fishbase_name)) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(mean_Loo = mean(Loo, na.rm = TRUE),
                     mean_K = mean(K, na.rm = TRUE),
                     mean_Winfinity = mean(Winfinity, na.rm = TRUE),
                     mean_tmax = mean(tmax, na.rm = TRUE),
                     mean_tm = mean(tm, na.rm = TRUE),
                     mean_M = mean(M, na.rm = TRUE),
                     mean_Lm = mean(Lm, na.rm = TRUE))

  # Trophic level
  LH_troph <- rfishbase::ecology(traj_SProd %>%
                                   dplyr::pull(fishbase_name)) %>%
    dplyr::select(Species, FoodTroph) %>%
    # https://github.com/ropensci/rfishbase/issues/199
    # "FoodTroph gives a MonteCarlo estimate of trophic level based on known food items.
    # DietTroph uses the mean or median of trophic levels derived from actual diet composition studies.
    # While in theory troph from diet should be more reliable,
    # many diet studies are incomplete or biased and I often find FoodTroph more reasonable." Rainer Froese
    dplyr::group_by(Species) %>%
    dplyr::slice(1) # make sure to keep only one row by species

  # Type of habitat (demersal, pelagic, reef)
  LH_hab <- rfishbase::fb_tbl("species") %>%
    dplyr::mutate(Species = paste(Genus, Species)) %>%
    dplyr::filter(Species %in% (traj_SProd %>% dplyr::pull(fishbase_name))) %>%
    dplyr::select(Species, DemersPelag) %>%
    dplyr::bind_rows(
      data.frame(Species=c("Ammodytes spp", "Merluccius spp", "Sebastes spp"),
                 DemersPelag = c("demersal", "bathydemersal", "bathypelagic"))) %>%
    dplyr::mutate(DemersPelag = dplyr::case_when(
      DemersPelag %in% c("bathydemersal", "benthopelagic", "demersal", "reef-associated") ~ "demersal",
      DemersPelag %in% c("pelagic-neritic", "pelagic-oceanic", "bathypelagic") ~ "pelagic")
      # DemersPelag %in% c("bathydemersal", "bathypelagic", "demersal", "reef-associated") ~ "demersal",
      # DemersPelag %in% c("pelagic-neritic", "pelagic-oceanic", "benthopelagic") ~ "pelagic")
      # https://fishbase.se/manual/English/FishBaseThe_Species_Table.htm
    )

  LH <-
    traj_SProd %>%
    dplyr::select(stockid, fishbase_name, classname, ordername, family) %>%
    dplyr::group_by(fishbase_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-stockid) %>%
    dplyr::rename(Species=fishbase_name) %>%
    dplyr::left_join(LH_grw, by="Species") %>%
    dplyr::left_join(LH_troph, by="Species") %>%
    dplyr::left_join(LH_hab, by="Species") %>%
    dplyr::select(Species, classname, ordername, family, FoodTroph, DemersPelag)

  ## FishLife:
  traits_fishlife <- FishLife::FishBase$beta_gv %>%
    tibble::as_tibble(rownames="taxa") %>%
    dplyr::mutate(taxa_type = dplyr::case_when(
      stringr::str_count(taxa, "_")==4 ~"species",
      stringr::str_count(taxa, "_")==3 ~"genus",
      stringr::str_count(taxa, "_")==2 ~"family",
      stringr::str_count(taxa, "_")==1 ~"order",
      stringr::str_count(taxa, "_")==0 ~"class")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Species = dplyr::case_when(
      taxa_type=="species" ~paste(strsplit(taxa, "_")[[1]][4:5], collapse=" "),
      taxa_type=="genus" ~strsplit(taxa, "_")[[1]][4],
      taxa_type=="family" ~strsplit(taxa, "_")[[1]][3],
      taxa_type=="order" ~strsplit(taxa, "_")[[1]][2],
      taxa_type=="class" ~strsplit(taxa, "_")[[1]][1]))

  fishlife_traits <- traj_SProd %>%
    dplyr::select(stockid, fishlife_name, phylum, classname, ordername, family) %>%
    dplyr::group_by(fishlife_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-stockid) %>%
    dplyr::rename(Species=fishlife_name) %>%
    dplyr::mutate(Species = sub(" spp", "", Species)) %>%
    dplyr::left_join(traits_fishlife %>%
                       dplyr::select(-taxa), by="Species")

  LH_fishlife_complete <-
    traj_SProd %>%
    dplyr::select(stockid, scientificname, fishbase_name, fishlife_name) %>%
    dplyr::group_by(scientificname) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-stockid) %>%
    dplyr::rename(Species=scientificname) %>%
    dplyr::left_join(fishlife_traits %>%
                       dplyr::select(-c(phylum, classname, ordername, family)),
                     by=c("fishlife_name"="Species")) %>%
    dplyr::left_join(LH, by=c("fishbase_name"="Species"))


  ## Fishing-related traits:
  fish_Umsy <-
    extract_RAM(traj_SProd %>% dplyr::pull(stockid),
                ts_type="UdivUmsypref", drop_na=TRUE) %>%
    dplyr::mutate(stockid = sub("UdivUmsypref_","",scen)) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, loc_brk_chg) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, loc_brk_chg),
      by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(
      max_U = max(UdivUmsypref),
      mean_U = mean(UdivUmsypref),
      sd_U = sd(UdivUmsypref),
      detrsd_U = sd(pracma::detrend(UdivUmsypref)),
      len_U = length(UdivUmsypref),
      D_U = (1/(len_U-1))*sum(abs(log((dplyr::lead(UdivUmsypref)+0.01*mean_U)/
                                        (UdivUmsypref+0.01*mean_U))), na.rm=TRUE)
    ) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>%
    dplyr::mutate(
      max_U_prshf = max(UdivUmsypref),
      mean_U_prshf = mean(UdivUmsypref),
      sd_U_prshf = sd(UdivUmsypref),
      detrsd_U_prshf = sd(pracma::detrend(UdivUmsypref)),
      len_U_prshf = length(UdivUmsypref),
      D_U_prshf = (1/(len_U-1))*sum(abs(log((dplyr::lead(UdivUmsypref)+0.01*mean_U)/
                                              (UdivUmsypref+0.01*mean_U))), na.rm=TRUE)
    ) %>%
    dplyr::slice(1) %>%
    tidyr::drop_na(mean_U) %>%
    dplyr::select(stockid,
                  max_U, mean_U,
                  sd_U, detrsd_U, D_U, len_U,
                  max_U_prshf, mean_U_prshf,
                  sd_U_prshf, detrsd_U_prshf, D_U_prshf, len_U_prshf
    )
  # There are two stocks ("COD2J3KL" and "RKCRABBB") for which the shift occurred
  # before values of UdivUmsypref were available.

  fish_ER <- traj_SProd %>%
    dplyr::pull(stockid) %>%
    extract_RAM(ts_type="ERbest2", drop_na=TRUE) %>%
    dplyr::mutate(stockid = sub("ERbest2_","",scen)) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, loc_brk_chg) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, loc_brk_chg),
      by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(mean_ER = mean(ERbest2),
                  len_ER = length(ERbest2)) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>%
    dplyr::mutate(mean_ER_prshf = mean(ERbest2),
                  ER_change_prshf =
                    summary(lm(ERbest2~year))$coefficients["year","Estimate"],
                  len_ER_prshf = length(ERbest2)) %>%
    dplyr::slice(1) %>%
    tidyr::drop_na(mean_ER_prshf) %>%
    dplyr::select(stockid,
                  mean_ER, mean_ER_prshf, ER_change_prshf, len_ER, len_ER_prshf)

  # Minimum & Average B/Bmsy:
  fish_Bmsy <-
    extract_RAM(traj_SProd %>% dplyr::pull(stockid),
                ts_type="BdivBmsypref", drop_na=TRUE) %>%
    dplyr::mutate(stockid = sub("BdivBmsypref_","",scen)) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, loc_brk_chg) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, loc_brk_chg),
      by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(
      min_B = min(BdivBmsypref),
      mean_B = mean(BdivBmsypref),
      sd_B = sd(BdivBmsypref),
      detrsd_B = sd(pracma::detrend(BdivBmsypref)),
      len_B = length(BdivBmsypref),
      D_B = (1/(len_B-1))*sum(abs(log((dplyr::lead(BdivBmsypref)+0.01*mean_B)/
                                        (BdivBmsypref+0.01*mean_B))), na.rm=TRUE)
    ) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>%
    dplyr::mutate(
      min_B_prshf = min(BdivBmsypref),
      mean_B_prshf = mean(BdivBmsypref),
      sd_B_prshf = sd(BdivBmsypref),
      detrsd_B_prshf = sd(pracma::detrend(BdivBmsypref)),
      len_B_prshf = length(BdivBmsypref),
      D_B_prshf = (1/(len_B-1))*sum(log((dplyr::lead(BdivBmsypref)+0.01*mean_B)/
                                          (BdivBmsypref+0.01*mean_B)), na.rm=TRUE)
    ) %>%
    dplyr::slice(1) %>%
    tidyr::drop_na(mean_B) %>%
    dplyr::select(stockid,
                  min_B, mean_B,
                  sd_B, detrsd_B, D_B, len_B,
                  min_B_prshf, mean_B_prshf,
                  sd_B_prshf, detrsd_B_prshf, D_B_prshf, len_B_prshf
    )

  # Total exploitation time:
  fish_expl_time <- timeseries %>%
    dplyr::group_by(stockid) %>%
    tidyr::drop_na(tsvalue) %>%
    dplyr::mutate(first_year_avail = min(tsyear)) %>%
    dplyr::filter(tsyear==first_year_avail) %>%
    dplyr::slice(1) %>%
    dplyr::select(stockid, first_year_avail) %>%
    dplyr::filter(stockid %in% (traj_SProd %>%
                                  dplyr::pull(stockid))) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, first, last, loc_brk_chg) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, first, last, loc_brk_chg),
      by="stockid") %>%
    dplyr::mutate(total_time = last-first_year_avail+1,
                  time_prshf = ifelse(!is.na(loc_brk_chg),
                                      loc_brk_chg-first_year_avail+1,
                                      total_time))

  # Absolute catch (mean or cumulative):
  fish_catch <- extract_RAM(traj_SProd %>% dplyr::pull(stockid),
                            ts_type="TCbest", drop_na=TRUE) %>%
    dplyr::mutate(stockid = sub("TCbest_","",scen)) %>%
    dplyr::left_join(
      traj_SProd %>%
        dplyr::select(stockid, scientificname, class, trend, loc_brk_chg) %>%
        dplyr::mutate(loc_brk_chg = ifelse(class == "abrupt", loc_brk_chg, NA)) %>%
        dplyr::select(stockid, loc_brk_chg),
      by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(
      sum_C = sum(TCbest, na.rm=TRUE),
      mean_C = mean(TCbest)
    ) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>%
    dplyr::mutate(
      sum_C_prshf = sum(TCbest, na.rm=TRUE),
      mean_C_prshf = mean(TCbest)
    ) %>%
    dplyr::slice(1) %>%
    tidyr::drop_na(mean_C) %>%
    dplyr::select(stockid,
                  sum_C, sum_C_prshf,
                  mean_C, mean_C_prshf)

  # Management approach
  fish_param <- traj_SProd %>%
    dplyr::select(stockid, scientificname) %>%
    dplyr::left_join(fish_expl_time, by="stockid") %>%
    dplyr::left_join(fish_Bmsy, by="stockid") %>%
    dplyr::left_join(fish_Umsy, by="stockid") %>%
    dplyr::left_join(fish_ER, by="stockid") %>%
    dplyr::left_join(fish_catch, by="stockid")


  # Environmental data:
  SST_had <-
    readr::read_csv("data/spatial_data/yearly_had_sst.csv") %>%
    dplyr::left_join(assessment %>%
                       dplyr::select(assessid, stockid), by="assessid") %>%
    dplyr::filter(stockid %in% (traj_SProd %>% dplyr::pull(stockid)))


  # Trim SST timeseries to match period available for productivity
  SST_had_SProd <- SST_had %>%
    dplyr::left_join(traj_SProd %>%
                       dplyr::select(stockid, first, last, loc_brk_chg, class), by="stockid") %>%
    dplyr::group_by(stockid) %>%
    dplyr::mutate(loc_brk_chg = ifelse(class=="abrupt", loc_brk_chg, NA)) %>%
    dplyr::filter(is.na(loc_brk_chg) | year<=loc_brk_chg) %>% # to match preshift period
    dplyr::filter(year>=first & year<=last) # to match period available

  sst_stock_list <- SST_had_SProd %>%
    tidyr::drop_na(sst_c) %>%
    dplyr::distinct(stockid) %>%
    dplyr::pull(stockid) %>%
    sort()

  SST_had_SProd_df <- SST_had_SProd %>%
    tidyr::drop_na(sst_c) %>%
    dplyr::mutate(sst_avg = mean(sst_c),
                  sst_cv = sd(sst_c)/mean(sst_c),
                  sst_sd = sd(sst_c),
                  sst_detrsd = sd(pracma::detrend(sst_c)),
                  sst_len = length(sst_c),
                  # sst_D = (1/(length(sst_c)-1))*sum(abs(log((dplyr::lead(sst_c)+0.01*sst_avg)/
                  #                                             (sst_c+0.01*sst_avg))), na.rm=TRUE),
                  sst_change =
                    summary(lm(sst_c~year))$coefficients["year","Estimate"]) %>%
    dplyr::slice(1) %>%
    dplyr::select(stockid, sst_avg, sst_change
                  # , sst_cv, sst_sd, sst_detrsd, sst_D
    )

  ## Polygons:
  polygons <- traj_SProd %>%
    dplyr::left_join(
      readr::read_csv(
        paste0("data/_spatial/ramldb_boundaries/",
               "ramldb_v3.8_stock_boundary_table_v2_formatted.csv")) %>%
        dplyr::mutate(batch="assess1") %>%
        dplyr::bind_rows(
          readr::read_csv(
            paste0("data/_spatial/ramldb_boundaries/",
                   "ramldb_stock_boundary_formatted_bis.csv")) %>%
            dplyr::mutate(batch="assess2")) %>%
        dplyr::rename(assessid_polygon_avail = assessid) %>%
        dplyr::group_by(stockid) %>%
        dplyr::mutate(n=n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n==1 | batch=="assess2") %>%
        dplyr::select(stockid, assessid_polygon_avail, batch),
      by="stockid")

  stocks_available_pol <- polygons %>%
    tidyr::drop_na(assessid_polygon_avail) %>%
    dplyr::pull(stockid)


  centroids_available <- lapply(stocks_available_pol, function(x){

    batch <- polygons %>%
      dplyr::filter(stockid==x) %>%
      dplyr::pull(batch)

    if(batch=="assess1"){
      path <- "shapes"
    } else {
      path <- "add_shapes"
    }

    sf::read_sf(paste0("data/_spatial/ramldb_boundaries/",path,"/",
                       polygons$assessid_polygon_avail[polygons$stockid==x],
                       ".shp")) %>%
      # transform to a projection with projected coordinates (World Mercator):
      sf::st_transform(3395) %>%
      sf::st_centroid() %>%
      # transform back to WGS 84 (with geographic coordinates):
      sf::st_transform(4326) %>%
      sf::st_coordinates() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(stockid=x) %>%
      suppressWarnings()

  }) %>%
    do.call(what=dplyr::bind_rows)

  polygons <- polygons %>%
    dplyr::left_join(centroids_available, by="stockid")

  # Merge all traits
  traits <- LH_fishlife_complete %>%
    dplyr::select(-c(taxa_type, fishbase_name, fishlife_name)) %>%
    dplyr::left_join(polygons %>%
                       dplyr::select(scientificname, stockid, X, Y) %>%
                       dplyr::rename(Species=scientificname), by="Species") %>%
    dplyr::left_join(fish_param %>%
                       dplyr::select(stockid, mean_B_prshf, mean_U_prshf,
                                     mean_ER_prshf, ER_change_prshf, sum_C_prshf),
                     by="stockid") %>%
    dplyr::left_join(SST_had_SProd_df %>%
                       dplyr::select(stockid, sst_avg, sst_change
                                     # , sst_sd, sst_D
                       ), by="stockid") %>%
    dplyr::relocate(stockid)

  return(traits)
}


# IV) Kobe plots ----------------------------------------------------------

#' Make the Kobe plot for a given stock
#'
#' @param id Short stock ID character.
#' @param start Start year to plot (default first year available).
#' @param end End year to plot (default last year available).
#' @param save Logical, whether to save Kobe plot or not.
#'
#' @return List of two objects:
#' - Kobe plot
#' - Data frame used to plot it
#' @export

kobe_plot <- function(id, start=NULL, end=NULL, save=FALSE){

  ts <- timeseries_values_views %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(stockid, year, BdivBmsypref, UdivUmsypref) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(kobe_zone = dplyr::case_when(
      BdivBmsypref >= 1 & UdivUmsypref < 1 ~ 1,
      BdivBmsypref >= 1 & UdivUmsypref >= 1 ~ 2,
      BdivBmsypref < 1 & UdivUmsypref < 1 ~ 3,
      BdivBmsypref < 1 & UdivUmsypref >= 1 ~ 4))

  if(!is.null(start)) ts <- ts %>% dplyr::filter(year>=start)
  if(!is.null(end)) ts <- ts %>% dplyr::filter(year<=end)

  metadt <- metadata %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(stockid, stocklong, scientificname, commonname, region)

  taxo <- taxonomy %>%
    dplyr::filter(scientificname == metadt$scientificname) %>%
    dplyr::select(scientificname, phylum, ordername, family, FisheryType)

  ts_info <- left_join(metadt, taxo, by = "scientificname")

  stocklong <- metadt$stocklong
  sciname <- ts_info$scientificname
  famname <- ts_info$family
  reg <- ts_info$region
  max_x <- ts %>% dplyr::select(BdivBmsypref) %>% max() %>% ceiling()
  max_y <- ts %>% dplyr::select(UdivUmsypref) %>% max() %>% ceiling()

  kobe <- ggplot(ts, aes(x=BdivBmsypref, y=UdivUmsypref
                         , color=year
                         ))+
    coord_cartesian(xlim=c(0, max_x), ylim=c(0, max_y))+
    annotate("rect", xmin=c(1,1,-Inf,-Inf), xmax=c(Inf,Inf,1,1),
             ymin=c(-Inf,1,-Inf,1) , ymax=c(1,Inf,1,Inf),
             fill=c("green","orange","yellow","red"), alpha=0.1)+
    geom_hline(yintercept=1) + geom_vline(xintercept=1)+
    scale_x_continuous(breaks=seq(0,max_x,1))+
    scale_y_continuous(breaks=seq(0,max_y,1))+
    geom_path(linewidth=1)+
    # geom_point(data=ts %>% dplyr::slice_tail(),
    #            aes(x=BdivBmsypref, y=UdivUmsypref),
    #            size=4, shape=23, fill="orange")+
    # geom_point(size=1, alpha=0.6)+
    theme_light()+
    scale_colour_gradient2(low = "blue", mid = "yellow" , high = "red",
                           midpoint=median(ts$year))+
    xlab(bquote("B/B"[MSY])) + ylab(bquote("U/U"[MSY]))+
    ggtitle(paste0(stocklong, " [",id,"] – (",ts$year[1],"–",
                   tail(ts$year, 1),")\n", sciname," – ",famname," – ",reg))

  if(save) ggsave(filename=paste0("analyses/kobe_plot/",id,"_kobe.png"),
                  plot = kobe, width=20, height = 15, units="cm")


  return(list("kobe" = kobe, "ts" = ts))
}



#' Display several stocks on the same Kobe plot
#'
#' @param ids List of short stock ID characters.
#' @param save Logical, whether to save Kobe plot or not.
#'
#' @return List of two objects:
#' - Kobe plot
#' - Data frame used to plot it
#' @export

kobe_plot_batch <- function(ids, save=FALSE){

  ts <- timeseries_values_views %>%
    dplyr::filter(stockid %in% ids) %>%
    dplyr::select(stockid, year, BdivBmsypref, UdivUmsypref) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(kobe_zone = dplyr::case_when(
      BdivBmsypref >= 1 & UdivUmsypref < 1 ~ 1,
      BdivBmsypref >= 1 & UdivUmsypref >= 1 ~ 2,
      BdivBmsypref < 1 & UdivUmsypref < 1 ~ 3,
      BdivBmsypref < 1 & UdivUmsypref >= 1 ~ 4))

  metadt <- metadata %>%
    dplyr::filter(stockid %in% ids) %>%
    dplyr::select(stockid, stocklong, scientificname, commonname, region)

  taxo <- taxonomy %>%
    dplyr::filter(scientificname %in% metadt$scientificname) %>%
    dplyr::select(scientificname, phylum, ordername, family, FisheryType)

  ts_info <- left_join(metadt, taxo, by = "scientificname")

  # stocklong <- metadt$stocklong
  # sciname <- ts_info$scientificname
  # famname <- ts_info$family
  # reg <- ts_info$region
  max_x <- ts %>% dplyr::select(BdivBmsypref) %>% max() %>% ceiling()
  max_y <- ts %>% dplyr::select(UdivUmsypref) %>% max() %>% ceiling()

  kobe <- ggplot(ts, aes(x=BdivBmsypref, y=UdivUmsypref, color=year, group=stockid))+
    coord_cartesian(xlim=c(0, 5), ylim=c(0, 10))+
    # coord_cartesian(xlim=c(0, max_x), ylim=c(0, max_y))+
    annotate("rect", xmin=c(1,1,-Inf,-Inf), xmax=c(Inf,Inf,1,1),
             ymin=c(-Inf,1,-Inf,1) , ymax=c(1,Inf,1,Inf),
             fill=c("green","orange","yellow","red"), alpha=0.1)+
    geom_hline(yintercept=1) + geom_vline(xintercept=1)+
    # scale_x_continuous(breaks=seq(0,max_x,1))+
    # scale_y_continuous(breaks=seq(0,max_y,1))+
    geom_path(alpha=0.5)+
    # geom_point(size=1, alpha=0.6)+
    theme_light()+
    scale_colour_gradient2(low = "blue", mid = "yellow" , high = "red",
                           midpoint=median(ts$year))+
    xlab(bquote("B/B"[MSY])) + ylab(bquote("U/U"[MSY]))

  if(save) ggsave(filename=paste0("analyses/kobe_plot/00_batch_kobe.png"),
                  plot = kobe, width=20, height = 15, units="cm")


  return(list("kobe" = kobe, "ts" = ts))
}



#' Make the Kobe plot for a given stock and pinpoint the location of a shift
#'
#' @param id Short stock ID characters.
#' @param traj_type Type of timeseries for which the shifts are searched for.
#' @param run_loo Logical, whether to perform leave-one-out process.
#' @param save Logical, whether to save Kobe plot or not.
#'
#' @return List of two objects:
#' - Kobe plot
#' - Data frame used to plot it
#' @export

kobe_plot_shift <- function(id, traj_type, run_loo=FALSE, save=FALSE){

  ts <- timeseries_values_views %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(stockid, year, BdivBmsypref, UdivUmsypref) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(kobe_zone = dplyr::case_when(
      BdivBmsypref >= 1 & UdivUmsypref < 1 ~ 1,
      BdivBmsypref >= 1 & UdivUmsypref >= 1 ~ 2,
      BdivBmsypref < 1 & UdivUmsypref < 1 ~ 3,
      BdivBmsypref < 1 & UdivUmsypref >= 1 ~ 4))

  metadt <- metadata %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(stockid, stocklong, scientificname, commonname, region)

  taxo <- taxonomy %>%
    dplyr::filter(scientificname == metadt$scientificname) %>%
    dplyr::select(scientificname, phylum, ordername, family, FisheryType)

  ts_info <- left_join(metadt, taxo, by = "scientificname")

  stocklong <- metadt$stocklong
  sciname <- ts_info$scientificname
  famname <- ts_info$family
  reg <- ts_info$region
  max_x <- ts %>% dplyr::select(BdivBmsypref) %>% max() %>% ceiling()
  max_y <- ts %>% dplyr::select(UdivUmsypref) %>% max() %>% ceiling()

  # Traj classif:
  sets <- extract_RAM(id, traj_type) %>%
    prep_data(thr=NULL, type="RAM", apriori=FALSE)

  traj <- traj_class(sets, str="aic_asd", abr_mtd=c("chg", "asd"),
                      asd_thr=0.15, type="RAM", showplots=TRUE,
                      apriori=FALSE, run_loo=run_loo, two_bkps=TRUE,
                      smooth_signif=TRUE, lowwl=5, highwl="default",
                      save_plot=FALSE, save_detplot=FALSE,
                      detection_plot=FALSE, edge_lim=0, congr_brk=5,
                      dirname="analyses/class_plot/",
                      outplot = TRUE, ind_plot="best")

  kobe <- ggplot(ts, aes(x=BdivBmsypref, y=UdivUmsypref, color=year))+
    coord_cartesian(xlim=c(0, max_x), ylim=c(0, max_y))+
    annotate("rect", xmin=c(1,1,-Inf,-Inf), xmax=c(Inf,Inf,1,1),
             ymin=c(-Inf,1,-Inf,1) , ymax=c(1,Inf,1,Inf),
             fill=c("green","orange","yellow","red"), alpha=0.1)+
    geom_hline(yintercept=1) + geom_vline(xintercept=1)+
    scale_x_continuous(breaks=seq(0,max_x,1))+
    scale_y_continuous(breaks=seq(0,max_y,1))+
    geom_path()+
    # geom_point(size=1, alpha=0.6)+
    theme_light()+
    scale_colour_gradient2(low = "blue", mid = "yellow" , high = "red",
                           midpoint=median(ts$year))+
    xlab(bquote("B/B"[MSY])) + ylab(bquote("U/U"[MSY]))+
    ggtitle(paste0(stocklong, " [",id,"] – (",ts$year[1],"–",
                   tail(ts$year, 1),")\n", sciname," – ",famname," – ",reg))



  if(traj$best_traj$class == "abrupt"){

    pred_chg <- traj$res_detail$res_abr$shifts_res[[1]]$chg_outlist$pred_chg

    sft_asd <- data.frame(year = as.numeric(traj$res_detail$res_abr$abr_res$chg$loc_brk),
                          mag = pred_chg$bp[length(pred_chg$bp)] - pred_chg$bp[1]) %>%
      dplyr::mutate(shift_sign = sign(mag))
    # %>% dplyr::slice(which.min(mag))
    kp_pan_asd <- ts %>%
      dplyr::filter(year %in% sft_asd$year) %>%
      dplyr::left_join(sft_asd, by = "year")

    kobe <- kobe +
      ggnewscale::new_scale_color()+
      geom_point(data=kp_pan_asd,
                 aes(x=BdivBmsypref, y=UdivUmsypref,
                     size=abs(mag), colour=as.factor(shift_sign)),
                 alpha=0.6)+
      scale_colour_manual(breaks=c("1", "-1"), values=c("blue","red"))+
      guides(colour="none", size="none")

  }

  kobe_traj <- cowplot::plot_grid(kobe, traj$class_plot$plots[[1]],
                                  ncol=1, rel_heights = c(2, 1))

  # } else {
  #
  #   kobe_traj <- cowplot::plot_grid(kobe, traj$class_plot[[1]],
  #                                   ncol=1, rel_heights = c(2, 1))
  # }

  if(save) ggsave(filename=paste0("analyses/kobe_plot_shift/",id,"_",traj_type,"_kobe_loo",run_loo,".png"),
                  plot = kobe_traj, width=15, height = 15, units="cm")

  return(list("kobe" = kobe, "traj" = traj, "kobe_traj" = kobe_traj, "ts" = ts))
}




# /!\ Compare with BdivBmsypref directly
kobe_plot_depr <- function(id){

  biop <- bioparams_values_views %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(stockid, TBmsybest, ERmsybest) %>%
    na.omit()

  if(nrow(biop)==0) stop("No parameters available for this stock")

  ts <- timeseries_values_views %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(stockid, year, TBbest, ERbest, BdivBmsypref, UdivUmsypref) %>%
    na.omit() %>%
    dplyr::mutate(TBbest_TBmsybest = TBbest/biop$TBmsybest,
                  ERbest_ERmsybest = ERbest/biop$ERmsybest,
                  kobe_zone = dplyr::case_when(TBbest_TBmsybest >= 1 & ERbest_ERmsybest < 1 ~ 1,
                                               TBbest_TBmsybest >= 1 & ERbest_ERmsybest >= 1 ~ 2,
                                               TBbest_TBmsybest < 1 & ERbest_ERmsybest < 1 ~ 3,
                                               TBbest_TBmsybest < 1 & ERbest_ERmsybest >= 1 ~ 4))

  metadt <- metadata %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(c(2,3,5,6,10))

  taxo <- taxonomy %>%
    dplyr::filter(scientificname == metadt[,3]) %>%
    dplyr::select(c(2,4,6,7,15))

  ts_info <- left_join(metadt, taxo, by = "scientificname")

  stocklong <- metadt$stocklong
  sciname <- ts_info$scientificname
  famname <- ts_info$family
  reg <- ts_info$region
  max_x <- ts %>% dplyr::select(TBbest_TBmsybest) %>% max() %>% ceiling()
  max_y <- ts %>% dplyr::select(ERbest_ERmsybest) %>% max() %>% ceiling()

  kobe <- ggplot(ts, aes(x=TBbest_TBmsybest, y=ERbest_ERmsybest, color=year))+
    coord_cartesian(xlim=c(0, max_x), ylim=c(0, max_y))+
    annotate("rect", xmin=c(1,1,-Inf,-Inf), xmax=c(Inf,Inf,1,1),
             ymin=c(-Inf,1,-Inf,1) , ymax=c(1,Inf,1,Inf),
             fill=c("green","orange","yellow","red"), alpha=0.1)+
    geom_hline(yintercept=1) + geom_vline(xintercept=1)+
    scale_x_continuous(breaks=seq(0,max_x,1))+
    scale_y_continuous(breaks=seq(0,max_y,1))+
    geom_path()+
    # geom_point(size=1, alpha=0.6)+
    theme_light()+
    scale_colour_gradient2(low = "blue", mid = "yellow" , high = "red", midpoint=median(ts$year))+
    xlab(bquote("TB / TB"[MSY])) + ylab(bquote("ER / ER"[MSY]))+
    ggtitle(paste0(stocklong, " [",id,"] – (",ts$year[1],"–",tail(ts$year, 1),")\n",
                   sciname," – ",famname," – ",reg))


  return(list("kobe" = kobe, "ts" = ts))
}



# V) Data display ---------------------------------------------------------

#' Plot sunburst pie chart
#'
#' @param df_sum a data frame with sums for trend and class
#' @param title
#' @param size
#' @param no_text
#' @param no_caption
#'
#' @return the data frame used to build the Kobe plot
#' @export


sunburst_plot <- function(df_sum, title=NULL, size=7,
                          no_text=FALSE, no_caption=FALSE){

  df_sum <- df_sum %>%
    dplyr::mutate(p = n/sum(n)*100)

  innerCircleData <- df_sum %>%
    dplyr::group_by(trend) %>%
    dplyr::summarise(tot_rev = sum(p)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot_rev),
      ymin = lag(ymax, n = 1, default = 0),
      csum = rev(cumsum(rev(tot_rev))),
      pos = tot_rev/2 + lead(csum, 1),
      pos = if_else(is.na(pos), tot_rev/2, pos)
    )

  innerCircle <- ggplot(innerCircleData)+
    geom_rect(aes(xmin = 0, xmax = 3, ymin = ymin, ymax = ymax, fill = trend),
              color = "white")+
    scale_fill_manual(values =
                        c("stable" = "#C7C860", "decrease" = "#F94042", "increase" = "#81B0FF",
                          "sta_no_change"="#E8DECB", "sta_quadratic"="#E8DE9C",
                          "dec_linear"="#FFBCBC", "dec_quadratic"="#FF6666", "dec_abrupt"="#660000",
                          "inc_linear"="#99CCFF", "inc_quadratic"="#6161F3", "inc_abrupt"="#000066"))+
    xlim(0, 4)+
    theme_bw()+
    coord_polar(theta = "y")

  if(!no_text){
    innerCircle <- innerCircle+
      geom_text(aes(x=1.7, y=tot_rev,
                    label = paste0(stringr::str_to_title(trend),"\n",
                                   round(tot_rev,digits=1),"%")),
                position = position_stack(vjust = 0.5),
                size=size)
  }

  outerCircleData <- df_sum %>%
    dplyr::group_by(trend, spe_class) %>%
    dplyr::summarise(tot = sum(p)) %>%
    dplyr::left_join(innerCircleData %>% select(trend, tot_rev), by="trend") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot),
      ymin = lag(ymax, n = 1, default = 0),
      pos = (ymax+ymin)/2
    ) %>%
    dplyr::filter(!is.na(spe_class))

  outerCircle <- innerCircle +
    geom_rect(
      data = outerCircleData,
      aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = spe_class),
      color = "white")

  if(!no_text){
    outerCircle <- outerCircle +
      ggrepel::geom_label_repel(data = outerCircleData,
                       aes(x=4, y = pos, label = paste0(stringr::str_to_title(
                         stringi::stri_sub(spe_class,5)),"\n",
                         round(tot, digits=1),"%")),
                       size = size,
                       show.legend = FALSE,
                       label.size=NA, fill=NA)
  }



  sunburst <- outerCircle+
    theme_void()+
    theme(legend.position="none",
          # plot.background = element_rect(colour=NA, fill="white"),
          # plot.margin = margin(1,4,1,4, "cm")
          plot.caption = element_text(vjust = 1, size=15),
          plot.title = element_text(vjust = 1, size=15))+
    ggtitle(title)

  if(!no_caption){

    sunburst <- sunburst+
      labs(caption = paste("N =", sum(df_sum$n)))
  }

  return(sunburst)

}


#' Plot sunburst pie chart by trajectory
#'
#' @param df_sum Data frame with sums for trend and class
#' @param title Character to display as title
#' @param size Numeric, text size
#' @param size_caption Numeric, caption text size
#' @param no_text Logical, to hide trajectory names and percentages
#' @param no_caption Logical, to hide caption with number of time series
#'
#' @return a sunburst plot
#' @export


sunburst_traj_plot <- function(df_sum, title=NULL, size=7, size_caption=13,
                          no_text=FALSE, no_caption=FALSE){

  df_sum <- df_sum %>%
    dplyr::mutate(p = n/sum(n)*100)

  innerCircleData <- df_sum %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(tot_rev = sum(p)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot_rev),
      ymin = lag(ymax, n = 1, default = 0),
      csum = rev(cumsum(rev(tot_rev))),
      pos = tot_rev/2 + lead(csum, 1),
      pos = if_else(is.na(pos), tot_rev/2, pos)
    )

  innerCircle <- ggplot(innerCircleData)+
    geom_rect(aes(xmin = 0, xmax = 3, ymin = ymin, ymax = ymax, fill = class),
              color = "white")+
    scale_fill_manual(
      values =
        c("no_change" = "grey80", "linear" = "grey60",
          "quadratic" = "grey40", "abrupt" = "grey30",
          "stable_constant"="#E8DECB", "stable_quadratic"="#E8DE9C",
          "decrease_linear"="#FFBCBC", "decrease_quadratic"="#FF6666", "decrease_abrupt"="#660000",
          "increase_linear"="#99CCFF", "increase_quadratic"="#6161F3", "increase_abrupt"="#000066"))+
    xlim(0, 4)+
    theme_bw()+
    coord_polar(theta = "y")

  if(!no_text){
    innerCircle <- innerCircle+
      geom_text(aes(x=1.7, y=tot_rev,
                    label = paste0(stringr::str_to_title(class),"\n",
                                   round(tot_rev,digits=1),"%")),
                position = position_stack(vjust = 0.5),
                size=size)
  }

  outerCircleData <- df_sum %>%
    dplyr::group_by(class, spe_class) %>%
    dplyr::summarise(tot = sum(p)) %>%
    dplyr::left_join(innerCircleData %>% select(class, tot_rev), by="class") %>%
    dplyr::left_join(df_sum %>% dplyr::select(trend, spe_class),  by="spe_class") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot),
      ymin = lag(ymax, n = 1, default = 0),
      pos = (ymax+ymin)/2
    ) %>%
    dplyr::filter(!is.na(spe_class))

  outerCircle <- innerCircle +
    geom_rect(
      data = outerCircleData,
      aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = spe_class),
      color = "white")

  if(!no_text){
    outerCircle <- outerCircle +
      ggrepel::geom_label_repel(data = outerCircleData,
                       aes(x=4, y = pos, label = paste0(stringr::str_to_title(
                         trend),"\n",
                         round(tot, digits=1),"%")),
                       size = size,
                       show.legend = FALSE,
                       label.size=NA, fill=NA)
  }



  sunburst <- outerCircle+
    theme_void()+
    theme(legend.position="none",
          # plot.background = element_rect(colour=NA, fill="white"),
          # plot.margin = margin(1,4,1,4, "cm")
          plot.caption = element_text(hjust = 0.8, vjust = 1, size=size_caption),
          plot.title = element_text(vjust = 1, size=15))+
    ggtitle(title)

  if(!no_caption){

    sunburst <- sunburst+
      labs(caption = paste("N =", sum(df_sum$n)))
  }

  return(sunburst)

}


#' Plot sunburst pie chart by trajectory (quadratic detailed)
#'
#' @param df_sum Data frame with sums for trend and class
#' @param title Character to display as title
#' @param size Numeric, text size
#' @param size_caption Numeric, caption text size
#' @param no_text Logical, to hide trajectory names and percentages
#' @param no_caption Logical, to hide caption with number of time series
#'
#' @return a sunburst plot
#' @export


sunburst_traj_plot_det <- function(df_sum, title=NULL, size=7, size_caption=15,
                               no_text=FALSE, no_caption=FALSE){

  df_sum <- df_sum %>%
    dplyr::mutate(p = n/sum(n)*100)

  innerCircleData <- df_sum %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(tot_rev = sum(p)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot_rev),
      ymin = lag(ymax, n = 1, default = 0),
      csum = rev(cumsum(rev(tot_rev))),
      pos = tot_rev/2 + lead(csum, 1),
      pos = if_else(is.na(pos), tot_rev/2, pos)
    )

  innerCircle <- ggplot(innerCircleData)+
    geom_rect(aes(xmin = 0, xmax = 3, ymin = ymin, ymax = ymax, fill = class),
              color = "white")+
    scale_fill_manual(values =
                        c("no_change" = "grey80", "linear" = "grey60",
                          "quadratic" = "grey40", "abrupt" = "grey30",
                          "stable_constant"="#E8DECB",
                          "stable_concave"="#BDAC28","stable_convex"="#F0E67D",
                          "decrease_constant"="#FFBCBC", "decrease_decelerated"="#FF6666",
                          "decrease_accelerated"="#FF6606", "decrease_abrupt"="#660000",
                          "increase_constant"="#99CCFF", "increase_decelerated"="#6161F3",
                          "increase_accelerated"="#3254CF", "increase_abrupt"="#000066"))+
    xlim(0, 4)+
    theme_bw()+
    coord_polar(theta = "y")

  if(!no_text){
    innerCircle <- innerCircle+
      geom_text(aes(x=1.7, y=tot_rev,
                    label = paste0(stringr::str_to_title(class),"\n",
                                   round(tot_rev,digits=1),"%")),
                position = position_stack(vjust = 0.5),
                size=size)
  }

  outerCircleData <- df_sum %>%
    dplyr::group_by(class, spe_class) %>%
    dplyr::summarise(tot = sum(p)) %>%
    dplyr::left_join(innerCircleData %>% select(class, tot_rev), by="class") %>%
    dplyr::left_join(df_sum %>% dplyr::select(trend, spe_class),  by="spe_class") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot),
      ymin = lag(ymax, n = 1, default = 0),
      pos = (ymax+ymin)/2
    ) %>%
    dplyr::filter(!is.na(spe_class)) %>%
    dplyr::mutate(traj_name = dplyr::case_when(
      class!="quadratic" ~ stringr::str_remove(spe_class, pattern = "\\_.*"),
      class=="quadratic" ~ stringr::str_replace(spe_class, "_", " "))
    )

  outerCircle <- innerCircle +
    geom_rect(
      data = outerCircleData,
      aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = spe_class),
      color = "white")

  if(!no_text){
    outerCircle <- outerCircle +
      ggrepel::geom_label_repel(data = outerCircleData,
                       aes(x=4, y = pos, label = paste0(stringr::str_to_title(
                         traj_name),"\n",
                         round(tot, digits=1),"%")),
                       size = size,
                       show.legend = FALSE,
                       label.size=NA, fill=NA)
  }



  sunburst <- outerCircle+
    theme_void()+
    theme(legend.position="none",
          # plot.background = element_rect(colour=NA, fill="white"),
          # plot.margin = margin(1,4,1,4, "cm")
          plot.caption = element_text(hjust = 0.8, vjust = 1, size=size_caption),
          plot.title = element_text(vjust = 1, size=15))+
    ggtitle(title)

  if(!no_caption){

    sunburst <- sunburst+
      labs(caption = paste("N =", sum(df_sum$n)))
  }

  return(sunburst)

}


#' Plot sunburst pie chart by trajectory (abrupt)
#'
#' @param df_sum Data frame with sums for trend and class
#' @param title Character to display as title
#' @param size Numeric, text size
#' @param size_caption Numeric, caption text size
#' @param no_text Logical, to hide trajectory names and percentages
#' @param no_caption Logical, to hide caption with number of time series
#'
#' @return a sunburst plot
#' @export


sunburst_traj_abr_plot <- function(df_sum, title=NULL, size=7, size_caption=15,
                               no_text=FALSE, no_caption=FALSE){

  df_sum <- df_sum %>%
    dplyr::mutate(p = n/sum(n)*100)

  innerCircleData <- df_sum %>%
    dplyr::group_by(traj_abr) %>%
    dplyr::summarise(tot_rev = sum(p)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot_rev),
      ymin = lag(ymax, n = 1, default = 0),
      csum = rev(cumsum(rev(tot_rev))),
      pos = tot_rev/2 + lead(csum, 1),
      pos = if_else(is.na(pos), tot_rev/2, pos)
    )

  innerCircle <- ggplot(innerCircleData)+
    geom_rect(aes(xmin = 0, xmax = 3, ymin = ymin, ymax = ymax, fill = traj_abr),
              color = "white")+
    scale_fill_manual(values =
                        c("abrupt" = "grey30", "not abrupt" = "grey70",
                          "sta_not abrupt"="#E3E1CC",
                          "dec_not abrupt"="#FABAB6", "dec_abrupt"="#660000",
                          "inc_not abrupt"="#ADD3FF", "inc_abrupt"="#000066"))+
    xlim(0, 4)+
    theme_bw()+
    coord_polar(theta = "y")

  if(!no_text){
    innerCircle <- innerCircle+
      geom_text(aes(x=1.7, y=tot_rev,
                    label = paste0(stringr::str_to_title(traj_abr),"\n",
                                   round(tot_rev,digits=1),"%")),
                position = position_stack(vjust = 0.5),
                size=size)
  }

  outerCircleData <- df_sum %>%
    dplyr::group_by(traj_abr, spe_class) %>%
    dplyr::summarise(tot = sum(p)) %>%
    dplyr::left_join(innerCircleData %>% select(traj_abr, tot_rev), by="traj_abr") %>%
    dplyr::left_join(df_sum %>% dplyr::select(trend, spe_class),  by="spe_class") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot),
      ymin = lag(ymax, n = 1, default = 0),
      pos = (ymax+ymin)/2
    ) %>%
    dplyr::filter(!is.na(spe_class))

  outerCircle <- innerCircle +
    geom_rect(
      data = outerCircleData,
      aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = spe_class),
      color = "white")

  if(!no_text){
    outerCircle <- outerCircle +
      ggrepel::geom_label_repel(data = outerCircleData,
                       aes(x=4, y = pos, label = paste0(stringr::str_to_title(
                         trend),"\n",
                         round(tot, digits=1),"%")),
                       size = size,
                       show.legend = FALSE,
                       label.size=NA, fill=NA)
  }



  sunburst <- outerCircle+
    theme_void()+
    theme(legend.position="none",
          # plot.background = element_rect(colour=NA, fill="white"),
          # plot.margin = margin(1,4,1,4, "cm")
          plot.caption = element_text(hjust = 0.8, vjust = 1, size=size_caption),
          plot.title = element_text(vjust = 1, size=15))+
    ggtitle(title)

  if(!no_caption){

    sunburst <- sunburst+
      labs(caption = paste("N =", sum(df_sum$n)))
  }

  return(sunburst)

}



#' Plot sunburst pie chart
#'
#' @param df_sum Data frame with sums for trend and class
#' @param title Character to display as title
#' @param size Numeric, text size
#' @param no_text Logical, to hide trajectory names and percentages
#' @param no_caption Logical, to hide caption with number of time series
#'
#' @return a sunburst plot
#' @export


pie_plot <- function(df_sum, title=NULL, size=7,
                          no_text=FALSE, no_caption=FALSE){

  library(ggrepel)

  df_sum <- df_sum %>%
    dplyr::mutate(p = n/sum(n)*100)

  innerCircleData <- df_sum %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(tot_rev = sum(p)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot_rev),
      ymin = lag(ymax, n = 1, default = 0),
      csum = rev(cumsum(rev(tot_rev))),
      pos = tot_rev/2 + lead(csum, 1),
      pos = if_else(is.na(pos), tot_rev/2, pos)
    )

  innerCircle <- ggplot(innerCircleData)+
    geom_rect(aes(xmin = 0, xmax = 3, ymin = ymin, ymax = ymax, fill = class),
              color = "white")+
    scale_fill_manual(values =
                        c("abrupt" = "yellow", "nonabrupt" = "grey",
                          "abrupt_dec"="indianred", "abrupt_inc"="lightblue",
                          " "="grey"))+
    xlim(0, 4)+
    theme_bw()+
    coord_polar(theta = "y")

  if(!no_text){
    innerCircle <- innerCircle+
      geom_text(aes(x=1.7, y=tot_rev,
                    label = paste0(stringr::str_to_title(class),"\n",
                                   round(tot_rev,digits=1),"%")),
                position = position_stack(vjust = 0.5),
                size=size)
  }

  outerCircleData <- df_sum %>%
    dplyr::group_by(class, spe_class) %>%
    dplyr::summarise(tot = sum(p)) %>%
    dplyr::left_join(innerCircleData %>% select(class, tot_rev), by="class") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ymax = cumsum(tot),
      ymin = lag(ymax, n = 1, default = 0),
      pos = (ymax+ymin)/2
    ) %>%
    dplyr::filter(!is.na(spe_class))

  outerCircle <- innerCircle +
    geom_rect(
      data = outerCircleData,
      aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = spe_class),
      color = "white")

  if(!no_text){
    outerCircle <- outerCircle +
      geom_label_repel(data = outerCircleData,
                       aes(x=4, y = pos, label = paste0(stringr::str_to_title(
                         stringi::stri_sub(spe_class,1)),"\n",
                         round(tot, digits=1),"%")),
                       size = size,
                       show.legend = FALSE,
                       label.size=NA, fill=NA)
  }



  sunburst <- outerCircle+
    theme_void()+
    theme(legend.position="none",
          # plot.background = element_rect(colour=NA, fill="white"),
          # plot.margin = margin(1,4,1,4, "cm")
          plot.caption = element_text(vjust = 1, size=15),
          plot.title = element_text(vjust = 1, size=15))+
    ggtitle(title)

  if(!no_caption){

    sunburst <- sunburst+
      labs(caption = paste("N =", sum(df_sum$n)))
  }

  return(sunburst)

}



#' Plot all summary plots
#'
#' @param traj a data frame output from traj_classification function.
#' @param df_list the list of timeseries classified.
#' @param var_name a character giving the name of the variable of interest.
#' @param loo a logical to specify whether LOO used.
#' @param filename a character to name the output plot.
#' @param title a string to display at the top of the figure.
#' @param caption a string specifying parameters used.
#'
#' @return No return value.
#' @export

# traj <-
#   readRDS("analyses/classif/classif_v4.61_SProd_minlen20_normFALSE_looFALSE.rds")
# var_name <- "Surplus production"

summary_plot <- function(traj, df_list, var_name, loo, filename, title, caption){

  # Class & trend proportions:
  traj_sum_det <- traj %>%
    dplyr::mutate(traj = dplyr::case_when(
      class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
      class=="abrupt" & trend=="increase" ~"increase_abrupt",
      class!="abrupt" ~traj
    )) %>%
    # dplyr::filter(stockid %in% kept_stocks) %>%
    dplyr::group_by(traj, trend, class) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::mutate(spe_class =
                    factor(traj,
                           levels = c("stable_constant",
                                      "stable_concave","stable_convex",
                                      "decrease_constant","decrease_decelerated",
                                      "decrease_accelerated","decrease_abrupt",
                                      "increase_constant","increase_decelerated",
                                      "increase_accelerated","increase_abrupt"))) %>%
    dplyr::ungroup()

  sunburst_traj <- sunburst_traj_plot_det(traj_sum_det, var_name, size=5)


  # Quality scores by class:
  quality_scores_df <- traj %>%
    dplyr::select(simu_id, dplyr::contains("weight_aic")) %>%
    tidyr::pivot_longer(cols=-simu_id, names_to = "score_class",
                        values_to = "wAICc") %>%
    dplyr::mutate(score_class = sub("weight_aic_","", score_class)) %>%
    dplyr::left_join(
      traj %>%
        dplyr::select(simu_id, dplyr::contains("nrmse")) %>%
        tidyr::pivot_longer(cols=-simu_id, names_to = "score_class",
                            values_to = "nrmse") %>%
        dplyr::mutate(score_class = sub("nrmse_","", score_class)),
      by=c("simu_id", "score_class")
    )
  score_types <- c("wAICc", "nrmse")

  if (loo == TRUE){
    quality_scores_df <- quality_scores_df %>%
      dplyr::left_join(
        traj %>%
          dplyr::select(simu_id, dplyr::contains("loo")) %>%
          tidyr::pivot_longer(cols=-simu_id, names_to = "score_class",
                              values_to = "loo") %>%
          dplyr::mutate(score_class = sub("loo_","", score_class)),
        by=c("simu_id", "score_class")
      )
    score_types <- c("wAICc", "nrmse", "loo")
  }

  quality_scores_df <- quality_scores_df %>%
    dplyr::left_join(traj %>%
                       dplyr::select(simu_id, class),
                     by="simu_id") %>%
    dplyr::filter(class==score_class) %>%
    dplyr::select(-class) %>%
    tidyr::pivot_longer(dplyr::all_of(score_types),
                        names_to = "score_type", values_to="score") %>%
    dplyr::mutate(score_class = factor(score_class,
                                       levels=c("no_change", "linear",
                                                "quadratic","abrupt")))

  quality_scores <- quality_scores_df %>%
    ggplot()+
    geom_boxplot(aes(x=score_class, y=score, fill=score_class),
                 outlier.shape = NA)+
    scale_fill_manual(values =
                        c("no_change" = "grey80", "linear" = "grey60",
                          "quadratic" = "grey40", "abrupt" = "grey30"))+
    geom_text(data = quality_scores_df %>%
                dplyr::filter(score_type=="nrmse") %>%
                dplyr::group_by(score_class) %>%
                dplyr::summarise(n=n()),
              aes(score_class, 0.05, label = paste0("n = ",n)),
              vjust = 1, size=3)+
    expand_limits(y=0)+
    labs(x="Class")+
    facet_wrap(vars(score_type), nrow=1)+
    theme_light()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
    ggtitle("Quality scores")


  # Quality scores by trajectory:
  spe_col <-
    c("sta_no_change"="#E8DECB", "sta_quadratic"="#E8DE9C",
      "dec_linear"="#FFBCBC", "dec_quadratic"="#FF6666", "dec_abrupt"="#660000",
      "inc_linear"="#99CCFF", "inc_quadratic"="#6161F3", "inc_abrupt"="#000066")

  quality_scores_traj_df <- traj %>%
    dplyr::select(simu_id, trend, traj, dplyr::contains("weight_aic")) %>%
    tidyr::pivot_longer(cols=-c(simu_id, trend, traj),
                        names_to = "score_class",
                        values_to = "wAICc") %>%
    dplyr::mutate(score_class = sub("weight_aic_","", score_class)) %>%
    dplyr::left_join(
      traj %>%
        dplyr::select(simu_id, dplyr::contains("nrmse")) %>%
        tidyr::pivot_longer(cols=-simu_id,
                            names_to = "score_class",
                            values_to = "nrmse") %>%
        dplyr::mutate(score_class = sub("nrmse_","", score_class)),
      by=c("simu_id", "score_class")
    )
  score_types <- c("wAICc", "nrmse")

  if (loo == TRUE){
    quality_scores_traj_df <- quality_scores_traj_df %>%
      dplyr::left_join(
        traj %>%
          dplyr::select(simu_id, dplyr::contains("loo")) %>%
          tidyr::pivot_longer(cols=-simu_id, names_to = "score_class",
                              values_to = "loo") %>%
          dplyr::mutate(score_class = sub("loo_","", score_class)),
        by=c("simu_id", "score_class")
      )
    score_types <- c("wAICc", "nrmse", "loo")
  }

  quality_scores_traj_df <- quality_scores_traj_df %>%
    dplyr::left_join(traj %>%
                       dplyr::select(simu_id, class),
                     by="simu_id") %>%
    dplyr::filter(class==score_class) %>%
    dplyr::select(-class) %>%
    tidyr::pivot_longer(dplyr::all_of(score_types),
                        names_to = "score_type", values_to="score") %>%
    dplyr::mutate(spe_class = factor(paste0(
      stringr::str_extract(trend, "^.{3}"), "_", score_class),
      levels = c("sta_no_change","sta_quadratic","dec_linear","dec_quadratic",
                 "dec_abrupt","inc_linear","inc_quadratic","inc_abrupt")))

  quality_scores_traj <- quality_scores_traj_df %>%
    ggplot()+
    geom_boxplot(aes(x=spe_class, y=score, fill=spe_class),
                 outlier.shape = NA)+
    scale_fill_manual(values = spe_col)+
    geom_text(data = quality_scores_traj_df %>%
                dplyr::filter(score_type=="nrmse") %>%
                dplyr::group_by(spe_class) %>%
                dplyr::summarise(n=n()),
              aes(spe_class, 0.05, label = paste0("n = ",n)),
              vjust = 1, size=3)+
    labs(x="Trajectory")+
    facet_wrap(vars(score_type), nrow=1)+
    theme_light()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1))

  # Temporal coverage of timeseries:
  ts_cover_df <- df_list %>%
    lapply(function(x)
      data.frame(first = unlist(x[1,1], use.names=FALSE),
                 last = unlist(x[nrow(x),1], use.names=FALSE),
                 gaps = ifelse(
                   unlist(x[nrow(x),1], use.names=FALSE) -
                     unlist(x[1,1], use.names=FALSE) + 1 == nrow(x),
                   FALSE, TRUE)
    )) %>% do.call(what=dplyr::bind_rows) %>%
    dplyr::arrange(first, -last) %>%
    dplyr::mutate(id = row_number())

  if(sapply(df_list[[1]], lubridate::is.Date)[[1]]==TRUE){
    ts_cover_df <- ts_cover_df %>%
      dplyr::mutate(first = lubridate::as_date(first),
                    last = lubridate::as_date(last))}

  ts_cover <- ts_cover_df %>%
    ggplot() +
    geom_segment(aes(x=id, xend=id, y=first, yend=last, color = gaps)) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "orangered")) +
    theme_light() +
    coord_flip() +
    xlab("Timeseries") +
    ylab(names(df_list[[1]])[1]) +
    theme(
      # panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      # axis.ticks.y = element_blank(),
      # axis.text.y = element_blank()
      ) +
    guides(color = "none")

  # Temporal distribution of abrupt shifts:
  if ("abrupt" %in% traj$class){

    temp_abr_sft_df <- traj %>%
      dplyr::filter(class=="abrupt") %>%
      dplyr::group_by(loc_brk_chg, trend) %>%
      dplyr::summarise(n_brk=n())

    if(sapply(df_list[[1]], lubridate::is.Date)[[1]]==TRUE){
      temp_abr_sft_df <- temp_abr_sft_df %>%
        dplyr::mutate(loc_brk_chg = lubridate::as_date(loc_brk_chg))
      }

    temp_abr_sft <- temp_abr_sft_df %>%
      ggplot()+
      geom_col(aes(x=loc_brk_chg, y=n_brk, fill=trend),
               inherit.aes = FALSE)+
      facet_wrap(vars(trend), ncol=1)+
      scale_fill_manual(values = c("decrease"="#660000", "increase"="#000066"))+
      labs(x="Breakdate", y="count")+
      theme_light()+
      theme(legend.position = "none")

    if(sapply(df_list[[1]], lubridate::is.Date)[[1]]==FALSE){
      temp_abr_sft <- temp_abr_sft +
        scale_x_continuous(breaks = scales::breaks_width(
          max(floor((max(traj$loc_brk_chg)-min(traj$loc_brk_chg))/8),1)
        ))
    }

  } else { temp_abr_sft <- ggplot() }



  # Distribution of abruptness of abrupt shifts:
  if ("abrupt" %in% traj$class){

    abrptnss_distr <- traj %>%
      dplyr::filter(class=="abrupt") %>%
      dplyr::select(abruptness) %>%
      dplyr::arrange(desc(abruptness)) %>%
      dplyr::mutate(sign = ifelse(abruptness>0, "increase", "decrease"),
                    x = row_number()/nrow(.)) %>%
      ggplot()+
      geom_segment( aes(x=x, xend=x, y=0, yend=abruptness, color=sign),
                    linewidth=1.3, alpha=0.9) +
      theme_light() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
      ) +
      scale_x_continuous(breaks = scales::breaks_width(0.1))+
      scale_colour_manual(values = c("decrease"="#660000", "increase"="#000066"))+
      expand_limits(y=-1)+
      labs(x="Proportion of abrupt timeseries", y="Abruptness")+
      coord_flip()

  } else { abrptnss_distr <- ggplot() }


  # Combine all plots:
  p <- (sunburst_traj + (quality_scores / quality_scores_traj)) /
    (ts_cover + abrptnss_distr + temp_abr_sft) +
    patchwork::plot_layout(nrow = 2, height = c(2, 1)) +
    plot_annotation(title = title,
                    caption = caption,
                    theme = theme(plot.title = element_text(size = 20)))


  ggsave(filename=paste0("analyses/datasets/",filename,"_summary_fig.png"),
         plot=p,
         width=22, height=14)

  return(invisible(NULL))
}


#' Plot all summary plots
#'
#' @param traj a data frame output from traj_classification function.
#' @param df_list the list of timeseries classified.
#' @param var_name a character giving the name of the variable of interest.
#' @param loo a logical to specify whether LOO used.
#' @param filename a character to name the output plot.
#' @param title a string to display at the top of the figure.
#' @param caption a string specifying parameters used.
#'
#' @return No return value.
#' @export

# traj <-
#   readRDS("analyses/classif/classif_v4.61_SProd_minlen20_normFALSE_looFALSE.rds")
# var_name <- "Surplus production"

summary_plot_light <- function(traj, df_list, var_name, loo, filename, title, caption){

  # Class & trend proportions:
  traj_sum_det <- traj %>%
    dplyr::mutate(traj = dplyr::case_when(
      class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
      class=="abrupt" & trend=="increase" ~"increase_abrupt",
      class!="abrupt" ~traj
    )) %>%
    # dplyr::filter(stockid %in% kept_stocks) %>%
    dplyr::group_by(traj, trend, class) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::mutate(spe_class =
                    factor(traj,
                           levels = c("stable_constant",
                                      "stable_concave","stable_convex",
                                      "decrease_constant","decrease_decelerated",
                                      "decrease_accelerated","decrease_abrupt",
                                      "increase_constant","increase_decelerated",
                                      "increase_accelerated","increase_abrupt"))) %>%
    dplyr::ungroup()

  sunburst_traj <- sunburst_traj_plot_det(traj_sum_det, var_name, size=5)


  # Temporal coverage of timeseries:
  ts_cover_df <- df_list %>%
    lapply(function(x)
      data.frame(first = unlist(x[1,1], use.names=FALSE),
                 last = unlist(x[nrow(x),1], use.names=FALSE),
                 gaps = ifelse(
                   unlist(x[nrow(x),1], use.names=FALSE) -
                     unlist(x[1,1], use.names=FALSE) + 1 == nrow(x),
                   FALSE, TRUE)
      )) %>% do.call(what=dplyr::bind_rows) %>%
    dplyr::arrange(first, -last) %>%
    dplyr::mutate(id = row_number())

  if(sapply(df_list[[1]], lubridate::is.Date)[[1]]==TRUE){
    ts_cover_df <- ts_cover_df %>%
      dplyr::mutate(first = lubridate::as_date(first),
                    last = lubridate::as_date(last))}

  ts_cover <- ts_cover_df %>%
    ggplot() +
    geom_segment(aes(x=id, xend=id, y=first, yend=last, color = gaps)) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "orangered")) +
    theme_light() +
    coord_flip() +
    xlab("Timeseries") +
    ylab(names(df_list[[1]])[1]) +
    theme(
      # panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      # axis.ticks.y = element_blank(),
      # axis.text.y = element_blank()
    ) +
    guides(color = "none")

  # Temporal distribution of abrupt shifts:
  if ("abrupt" %in% traj$class){

    temp_abr_sft_df <- traj %>%
      dplyr::filter(class=="abrupt") %>%
      dplyr::group_by(loc_brk_chg, trend) %>%
      dplyr::summarise(n_brk=n())

    if(sapply(df_list[[1]], lubridate::is.Date)[[1]]==TRUE){
      temp_abr_sft_df <- temp_abr_sft_df %>%
        dplyr::mutate(loc_brk_chg = lubridate::as_date(loc_brk_chg))
    }

    temp_abr_sft <- temp_abr_sft_df %>%
      ggplot()+
      geom_col(aes(x=loc_brk_chg, y=n_brk, fill=trend),
               inherit.aes = FALSE)+
      facet_wrap(vars(trend), ncol=1)+
      scale_fill_manual(values = c("decrease"="#660000", "increase"="#000066"))+
      labs(x="Breakdate", y="count")+
      expand_limits(x=c(min(ts_cover_df$first), max(ts_cover_df$last)))+
      theme_light()+
      theme(legend.position = "none")

    if(sapply(df_list[[1]], lubridate::is.Date)[[1]]==FALSE){
      temp_abr_sft <- temp_abr_sft +
        scale_x_continuous(breaks = scales::breaks_width(
          max(floor((max(traj$loc_brk_chg)-min(traj$loc_brk_chg))/8),1)
        ))
    }

  } else { temp_abr_sft <- ggplot() }


  # Combine all plots:
  p <- sunburst_traj / (temp_abr_sft / ts_cover)+
    plot_annotation(title = title,
                    caption = caption,
                    theme = theme(plot.title = element_text(size = 20)))

  return(p)
}



#' Summary plots for a given stock
#'
#' @param stk Stock name.
#' @param drop_SProd Display years without productivity?
#'
#' @return A patchwork of plots.
#' @export

stock_summary <- function(stk, drop_SProd=FALSE){

  stock_ts <-

    # Start from SProd and biomass timeseries:
    timeseries_values_views %>%
    dplyr::filter(stockid==stk) %>%
    dplyr::select(stockid, year, SProd, TBbest) %>%

    # Add exploitation rate timeseries:
    dplyr::left_join(
      timeseries_values_views %>%
        dplyr::filter(stockid==stk) %>%
        dplyr::select(stockid, year, ERbest2, UdivUmsypref, BdivBmsypref) %>%
        dplyr::rename(ERbest = ERbest2),
      by=c("stockid", "year")
    ) %>%

    # Add SST timeseries:
    dplyr::right_join(
      readr::read_csv(paste0("data/SST_share_Mathieu/Annual_SST_code/Output/",
                             "1e_yearly_had_sst_polMP.csv")) %>%
        dplyr::left_join(assessment %>%
                           dplyr::select(assessid, stockid), by="assessid") %>%
        dplyr::filter(stockid==stk) %>%
        dplyr::left_join(traj_SProd %>%
                           dplyr::select(stockid, first, last, length),
                         by="stockid") %>%
        dplyr::filter(year>=first & year<=last) %>%
        dplyr::select(stockid, year, sst_c),
      by=c("stockid", "year")
    ) %>%

    # Add trajectory:
    dplyr::left_join(
      traj_SProd %>%
        dplyr::mutate(traj = dplyr::case_when(
          class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
          class=="abrupt" & trend=="increase" ~"increase_abrupt",
          class!="abrupt" ~traj),
          traj = plyr::revalue(traj,
                               c("increase_constant"="increase_linear",
                                 "decrease_constant"="decrease_linear"))) %>%
        dplyr::mutate(traj = sub("_", " ", traj),
                      traj =
                        factor(traj,
                               levels =
                                 c("stable constant",
                                   "stable concave","stable convex",
                                   "increase linear","increase decelerated",
                                   "increase accelerated","increase abrupt",
                                   "decrease linear","decrease decelerated",
                                   "decrease accelerated","decrease abrupt"
                                 ))) %>%
        dplyr::mutate(loc_brk_chg_nospe = loc_brk_chg,
                      loc_brk_chg=ifelse(class=="abrupt", loc_brk_chg, NA),
                      loc_brk_asd_1=ifelse(class=="abrupt", loc_brk_asd_1, NA),
                      loc_brk_asd_2=ifelse(class=="abrupt", loc_brk_asd_2, NA)) %>%
        dplyr::select(stockid, traj, trend, loc_brk_chg, loc_brk_chg_nospe,
                      loc_brk_asd_1, loc_brk_asd_2) %>%
        dplyr::filter(stockid==stk),
      by="stockid"
    ) %>%

    # Add year of first collapse:
    # dplyr::left_join(coll %>%
    #                    dplyr::select(stockid, first_coll_y), by="stockid") %>%

    # Transform in years to collapse or to productivity abrupt shift (PAS):
    # dplyr::mutate(year_to_coll = year-first_coll_y,
    #               year_to_PAS = year-loc_brk_chg_nospe) %>%

    # SST compared to SST in the first year:
    dplyr::mutate(sst_norm = sst_c-dplyr::first(sst_c)) %>%

    # Pinpoint shift year:
    dplyr::mutate(shift_y = ifelse(is.na(loc_brk_chg) | loc_brk_chg!=year, FALSE, TRUE))

  if(drop_SProd==TRUE){
    stock_ts <- stock_ts %>%
      tidyr::drop_na(SProd)
    }

  mean_TB <- mean(stock_ts %>%
                    tidyr::drop_na(SProd) %>%
                    pull(TBbest))

  traits <- read.csv("analyses/classif/mod_df01_default.csv") %>%
    dplyr::mutate(that_stock = ifelse(stockid==stk, TRUE, FALSE))


  # Productivity vs. time:
  if(!all(is.na(stock_ts$SProd))){

    SProd_traj <- stock_ts %>%
      dplyr::select(stockid, year, SProd) %>%
      dplyr::mutate(SProd=SProd/mean_TB*100) %>%
      prep_data(thr=NULL, type="data", apriori=FALSE) %>%
      traj_class(sets=., str="aic_asd", abr_mtd=c("chg", "asd"), asd_thr=0.15,
                 type="data", showplots=FALSE,
                 apriori=FALSE, run_loo=FALSE, two_bkps=FALSE,
                 smooth_signif=TRUE, outplot=FALSE,
                 ind_plot=NULL, lowwl=5, highwl="default",
                 congr_brk=5, dirname=NULL,
                 save_plot=FALSE)

    p_sprod <- stock_ts %>%
      ggplot(aes(x=year, y=SProd/mean_TB*100))+
      geom_line()+
      geom_point()+
      geom_line(aes(x=X, y=fit), col="blue", data=SProd_traj$best_fit, linewidth=1)+
      scale_y_continuous(
        name = "Normalized productivity\n(% mean total biomass)",
        sec.axis = sec_axis(transform=~.*mean_TB/100,
                            name="Raw productivity"))+
      geom_hline(yintercept = 0, lty="dashed")+
      geom_vline(xintercept =
                   if(SProd_traj$best_traj$class=="abrupt")
                     SProd_traj$res$loc_brk_asd %>%
                   strsplit(";") %>% `[[`(1) %>%
                   as.numeric() %>%
                   suppressWarnings()
                 else as.numeric(NA),
                 col="red", lty="dotted", linewidth=1)+
      geom_vline(xintercept =
                   if(SProd_traj$best_traj$class=="abrupt")
                     as.numeric(SProd_traj$res$loc_brk_chg)
                 else as.numeric(NA),
                 col="blue", lty="dashed", linewidth=1)+
      theme_light()+
      labs(x="")+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_sprod <- stock_ts %>%
      ggplot(aes(x=year, y=SProd/mean_TB*100))+
      scale_y_continuous(
        name = "Normalized productivity\n(% mean total biomass)",
        sec.axis = sec_axis(transform=~.*mean_TB/100,
                            name="Raw productivity"))+
      labs(x="")+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  }

  # Productivity vs. biomass
  p_sprod_TB <- stock_ts %>%
    ggplot(aes(x=TBbest, y=SProd, color=year))+
    geom_vline(xintercept = 0, lty="dashed")+
    geom_hline(yintercept = 0, lty="dashed")+
    geom_path(linewidth=2)+
    geom_point(aes(shape=shift_y, size=shift_y), fill="red", stroke=1)+
    scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 16)) +
    scale_size_manual(values = c("TRUE" = 8, "FALSE" = 3)) +
    theme_light()+
    labs(x="Total biomass", y="Raw productivity")+
    scale_color_continuous(type = "viridis")+
    # scale_color_distiller(palette = "Spectral")+
    guides(shape="none", size="none")+
    theme(panel.grid=element_blank(),
          panel.background = element_rect(fill="grey90"),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          legend.position="bottom")

  # Biomass relative to MSY vs. time:
  if(!all(is.na(stock_ts$BdivBmsypref))){

    BdivBmsypref_traj <- stock_ts %>%
      dplyr::select(stockid, year, BdivBmsypref) %>%
      tidyr::drop_na(BdivBmsypref) %>%
      prep_data(thr=NULL, type="data", apriori=FALSE) %>%
      traj_class(sets=., str="aic_asd", abr_mtd=c("chg", "asd"), asd_thr=0.15,
                 type="data", showplots=FALSE,
                 apriori=FALSE, run_loo=FALSE, two_bkps=FALSE,
                 smooth_signif=TRUE, outplot=FALSE,
                 ind_plot=NULL, lowwl=5, highwl="default",
                 congr_brk=5, dirname=NULL,
                 save_plot=FALSE)

    p_bbmsy <- stock_ts %>%
      ggplot(aes(x=year, y=BdivBmsypref))+
      geom_line()+
      geom_point()+
      labs(y="B/Bmsy", x="")+
      geom_line(aes(x=X, y=fit), col="blue", data=BdivBmsypref_traj$best_fit, linewidth=1)+
      geom_hline(yintercept = 1, lty="dashed")+
      geom_vline(xintercept =
                   if(BdivBmsypref_traj$best_traj$class=="abrupt")
                     BdivBmsypref_traj$res$loc_brk_asd %>%
                   strsplit(";") %>% `[[`(1) %>%
                   as.numeric() %>%
                   suppressWarnings()
                 else as.numeric(NA),
                 col="red", lty="dotted", linewidth=1)+
      geom_vline(xintercept =
                   if(BdivBmsypref_traj$best_traj$class=="abrupt")
                     as.numeric(BdivBmsypref_traj$res$loc_brk_chg)
                 else as.numeric(NA),
                 col="blue", lty="dashed", linewidth=1)+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_bbmsy <- stock_ts %>%
      ggplot(aes(x=year, y=BdivBmsypref))+
      labs(y="B/Bmsy", x="")+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  }

  # Fishing intensity relative to MSY vs. time:
  if(!all(is.na(stock_ts$UdivUmsypref))){

    UdivUmsypref_traj <- stock_ts %>%
      dplyr::select(stockid, year, UdivUmsypref) %>%
      tidyr::drop_na(UdivUmsypref) %>%
      prep_data(thr=NULL, type="data", apriori=FALSE) %>%
      traj_class(sets=., str="aic_asd", abr_mtd=c("chg", "asd"), asd_thr=0.15,
                 type="data", showplots=FALSE,
                 apriori=FALSE, run_loo=FALSE, two_bkps=FALSE,
                 smooth_signif=TRUE, outplot=FALSE,
                 ind_plot=NULL, lowwl=5, highwl="default",
                 congr_brk=5, dirname=NULL,
                 save_plot=FALSE)

    p_uumsy <- stock_ts %>%
      ggplot(aes(x=year, y=UdivUmsypref))+
      geom_line()+
      geom_point()+
      labs(y="U/Umsy", x="")+
      geom_line(aes(x=X, y=fit), col="blue", data=UdivUmsypref_traj$best_fit, linewidth=1)+
      geom_hline(yintercept = 1, lty="dashed")+
      geom_vline(xintercept =
                   if(UdivUmsypref_traj$best_traj$class=="abrupt")
                     UdivUmsypref_traj$res$loc_brk_asd %>%
                   strsplit(";") %>% `[[`(1) %>%
                   as.numeric() %>%
                   suppressWarnings()
                 else as.numeric(NA),
                 col="red", lty="dotted", linewidth=1)+
      geom_vline(xintercept =
                   if(UdivUmsypref_traj$best_traj$class=="abrupt")
                     as.numeric(UdivUmsypref_traj$res$loc_brk_chg)
                 else as.numeric(NA),
                 col="blue", lty="dashed", linewidth=1)+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_uumsy <- stock_ts %>%
      ggplot(aes(x=year, y=UdivUmsypref))+
      labs(y="U/Umsy", x="")+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))
  }

  # Exploitation rate vs. time:
  if(!all(is.na(stock_ts$ERbest))){

    ERbest_traj <- stock_ts %>%
      dplyr::select(stockid, year, ERbest) %>%
      tidyr::drop_na(ERbest) %>%
      prep_data(thr=NULL, type="data", apriori=FALSE) %>%
      traj_class(sets=., str="aic_asd", abr_mtd=c("chg", "asd"), asd_thr=0.15,
                 type="data", showplots=FALSE,
                 apriori=FALSE, run_loo=FALSE, two_bkps=FALSE,
                 smooth_signif=TRUE, outplot=FALSE,
                 ind_plot=NULL, lowwl=5, highwl="default",
                 congr_brk=5, dirname=NULL,
                 save_plot=FALSE)

    p_er <- stock_ts %>%
      ggplot(aes(x=year, y=ERbest))+
      geom_line()+
      geom_point()+
      labs(y="Exploitation rate", x="")+
      geom_line(aes(x=X, y=fit), col="blue", data=ERbest_traj$best_fit, linewidth=1)+
      geom_vline(xintercept =
                   if(ERbest_traj$best_traj$class=="abrupt")
                     ERbest_traj$res$loc_brk_asd %>%
                   strsplit(";") %>% `[[`(1) %>%
                   as.numeric() %>%
                   suppressWarnings()
                 else as.numeric(NA),
                 col="red", lty="dotted", linewidth=1)+
      geom_vline(xintercept =
                   if(ERbest_traj$best_traj$class=="abrupt")
                     as.numeric(ERbest_traj$res$loc_brk_chg)
                 else as.numeric(NA),
                 col="blue", lty="dashed", linewidth=1)+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_er <- stock_ts %>%
      ggplot(aes(x=year, y=ERbest))+
      labs(y="Exploitation rate", x="")+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))
  }

  # Sea surface temperature vs. time:
  if(!all(is.na(stock_ts$sst_c))){

    SST_traj <- stock_ts %>%
      dplyr::select(stockid, year, sst_c) %>%
      tidyr::drop_na(sst_c) %>%
      prep_data(thr=NULL, type="data", apriori=FALSE) %>%
      traj_class(sets=., str="aic_asd", abr_mtd=c("chg", "asd"), asd_thr=0.15,
                 type="data", showplots=FALSE,
                 apriori=FALSE, run_loo=FALSE, two_bkps=FALSE,
                 smooth_signif=TRUE, outplot=FALSE,
                 ind_plot=NULL, lowwl=5, highwl="default",
                 congr_brk=5, dirname=NULL,
                 save_plot=FALSE)

    p_sst <- stock_ts %>%
      ggplot(aes(x=year, y=sst_c))+
      geom_line()+
      geom_point()+
      labs(y="Sea surface temperature", x="")+
      geom_line(aes(x=X, y=fit), col="blue", data=SST_traj$best_fit, linewidth=1)+
      geom_vline(xintercept =
                   if(SST_traj$best_traj$class=="abrupt")
                     SST_traj$res$loc_brk_asd %>%
                   strsplit(";") %>% `[[`(1) %>%
                   as.numeric() %>%
                   suppressWarnings()
                 else as.numeric(NA),
                 col="red", lty="dotted", linewidth=1)+
      geom_vline(xintercept =
                   if(SST_traj$best_traj$class=="abrupt")
                     as.numeric(SST_traj$res$loc_brk_chg)
                 else as.numeric(NA),
                 col="blue", lty="dashed", linewidth=1)+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_sst <- stock_ts %>%
      ggplot(aes(x=year, y=sst_c))+
      labs(y="Sea surface temperature", x="")+
      theme_light()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))
  }

  # Lm vs. tm:
  p_lh <- traits %>%
    ggplot(aes(x=exp(Lm), y=exp(tm), colour=that_stock, size=that_stock, shape=that_stock))+
    geom_point(alpha=0.7)+
    labs(x="Length at maturity (cm)", y="Age at maturity (years)")+
    scale_colour_manual(values = c("TRUE" = "blue", "FALSE" = "grey20"))+
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
    scale_size_manual(values = c("TRUE" = 10, "FALSE" = 2)) +
    theme_light()+
    guides(colour="none", shape="none", size="none")+
    theme(panel.grid=element_blank(),
          panel.background = element_rect(fill="grey90"),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          legend.position="bottom")

  # Polygon:
  shpdir1 <- "data/SST_share_Mathieu/Stock_boundary/ramldb_boundaries"
  shpdir2 <- "data/SST_share_Mathieu/Stock_boundary/add_shapes"

  pol_info <- assessment %>%
    dplyr::filter(
      assessid %in%
        (list.files(shpdir1, pattern=".shp") %>%
        sub(".shp", "", .))) %>%
    dplyr::filter(stockid==stk) %>%
    dplyr::mutate(path = shpdir1) %>%
    dplyr::bind_rows(
      assessment %>%
        dplyr::filter(
          assessid %in%
            (list.files(shpdir2, pattern=".shp") %>%
               sub(".shp", "", .))) %>%
        dplyr::filter(stockid==stk) %>%
        dplyr::mutate(path = shpdir2)
    )

  polyg <- sf::read_sf(dsn=pol_info$path,
              layer=pol_info$assessid)

  world <- map_data('world') %>% fortify()
  p_pol <- ggplot() +
    geom_map(data = world,
             map = world,
             aes(group = group, map_id=region),
             fill = "#7f7f7f", colour = "#7f7f7f", size=0.5)+
    geom_sf(data = polyg, fill="blue")+
    coord_sf(xlim = c(sf::st_bbox(polyg)[1]-30, sf::st_bbox(polyg)[3]+30),
             ylim = c(sf::st_bbox(polyg)[2]-10, sf::st_bbox(polyg)[4]+10),
             expand = FALSE)+
    theme_light()+
    theme(panel.grid=element_blank(),
          axis.text = element_text(size=15))


  # Fish silhouette:
  sp_name <- as.data.frame(stock) %>%
    dplyr::filter(stockid==stk) %>%
    dplyr::pull(scientificname)

  sp_silh_file <- grep(sp_name,
                       list.files(paste0("../../../res/silhouettes/")),
                       value=TRUE) %>%
    grep( ".png", ., value=TRUE)

  if(length(sp_silh_file)==0){
    sp_silh_file <- "fish_silh.png"
  }

  path_to_shape <- paste0("../../../res/silhouettes/", sp_silh_file)

  img <- png::readPNG(path_to_shape)
  g <- grid::rasterGrob(img, interpolate=TRUE)

  spname <- stock %>%
    as.data.frame() %>%
    dplyr::filter(stockid==stk) %>%
    dplyr::slice(1) %>%
    dplyr::left_join(taxonomy, by="scientificname")

  silh <- ggplot()+
    ggplot2::annotation_custom(g)+
    theme_void()

  text <- ggplot()+
    theme_void()+
    annotate("text", x = 0, y = 0, size=8,
             label=paste0(spname$scientificname, "\n", spname$family))

  p_silh <- silh / text + patchwork::plot_layout(height = c(3, 1))

  # Combine everything:

  design <- "
  12
  34
  54
  67
  89
"
  fig_stock_summary <-
    p_sprod +
    p_silh +
    p_bbmsy +
    p_pol+
    p_uumsy +
    p_er +
    p_lh +
    p_sst +
    p_sprod_TB +
    patchwork::plot_layout(design=design)+
    patchwork::plot_annotation(
      title=paste0(assessment %>%
                     dplyr::filter(stockid==stk) %>%
                     dplyr::slice(1) %>%
                     dplyr::pull(stocklong),
                   " (", stk, ")"),
      theme=theme(plot.title = element_text(hjust = 0.5, size=35)))

  return(fig_stock_summary)
}





stock_sprob_biom <- function(stk, drop_SProd=FALSE){

  stock_ts <-

    # Start from SProd and biomass timeseries:
    timeseries_values_views %>%
    dplyr::filter(stockid==stk) %>%
    dplyr::select(stockid, year, SProd, TBbest) %>%

    # Add exploitation rate timeseries:
    dplyr::left_join(
      timeseries_values_views %>%
        dplyr::filter(stockid==stk) %>%
        dplyr::select(stockid, year, ERbest2, UdivUmsypref, BdivBmsypref) %>%
        dplyr::rename(ERbest = ERbest2),
      by=c("stockid", "year")
    ) %>%

    # Add SST timeseries:
    dplyr::right_join(
      readr::read_csv(paste0("data/SST_share_Mathieu/Annual_SST_code/Output/",
                             "1e_yearly_had_sst_polMP.csv")) %>%
        dplyr::left_join(assessment %>%
                           dplyr::select(assessid, stockid), by="assessid") %>%
        dplyr::filter(stockid==stk) %>%
        dplyr::left_join(traj_SProd %>%
                           dplyr::select(stockid, first, last, length),
                         by="stockid") %>%
        dplyr::filter(year>=first & year<=last) %>%
        dplyr::select(stockid, year, sst_c),
      by=c("stockid", "year")
    ) %>%

    # Add trajectory:
    dplyr::left_join(
      traj_SProd %>%
        dplyr::mutate(traj = dplyr::case_when(
          class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
          class=="abrupt" & trend=="increase" ~"increase_abrupt",
          class!="abrupt" ~traj),
          traj = plyr::revalue(traj,
                               c("increase_constant"="increase_linear",
                                 "decrease_constant"="decrease_linear"))) %>%
        dplyr::mutate(traj = sub("_", " ", traj),
                      traj =
                        factor(traj,
                               levels =
                                 c("stable constant",
                                   "stable concave","stable convex",
                                   "increase linear","increase decelerated",
                                   "increase accelerated","increase abrupt",
                                   "decrease linear","decrease decelerated",
                                   "decrease accelerated","decrease abrupt"
                                 ))) %>%
        dplyr::mutate(loc_brk_chg_nospe = loc_brk_chg,
                      loc_brk_chg=ifelse(class=="abrupt", loc_brk_chg, NA),
                      loc_brk_asd_1=ifelse(class=="abrupt", loc_brk_asd_1, NA),
                      loc_brk_asd_2=ifelse(class=="abrupt", loc_brk_asd_2, NA)) %>%
        dplyr::select(stockid, traj, trend, loc_brk_chg, loc_brk_chg_nospe,
                      loc_brk_asd_1, loc_brk_asd_2) %>%
        dplyr::filter(stockid==stk),
      by="stockid"
    ) %>%

    # Add year of first collapse:
    # dplyr::left_join(coll %>%
    #                    dplyr::select(stockid, first_coll_y), by="stockid") %>%

    # Transform in years to collapse or to productivity abrupt shift (PAS):
    # dplyr::mutate(year_to_coll = year-first_coll_y,
    #               year_to_PAS = year-loc_brk_chg_nospe) %>%

    # SST compared to SST in the first year:
    dplyr::mutate(sst_norm = sst_c-dplyr::first(sst_c)) %>%

    # Pinpoint shift year:
    dplyr::mutate(shift_y = ifelse(is.na(loc_brk_chg) | loc_brk_chg!=year, FALSE, TRUE))

  if(drop_SProd==TRUE){
    stock_ts <- stock_ts %>%
      tidyr::drop_na(SProd)
  }

  mean_TB <- mean(stock_ts %>%
                    tidyr::drop_na(SProd) %>%
                    pull(TBbest))

  traits <- read.csv("analyses/classif/mod_df01_default.csv") %>%
    dplyr::mutate(that_stock = ifelse(stockid==stk, TRUE, FALSE))


  # Productivity vs. time:
  if(!all(is.na(stock_ts$SProd))){

    SProd_traj <- stock_ts %>%
      dplyr::select(stockid, year, SProd) %>%
      dplyr::mutate(SProd=SProd/mean_TB*100) %>%
      prep_data(thr=NULL, type="data", apriori=FALSE) %>%
      traj_class(sets=., str="aic_asd", abr_mtd=c("chg", "asd"), asd_thr=0.15,
                 type="data", showplots=FALSE,
                 apriori=FALSE, run_loo=FALSE, two_bkps=FALSE,
                 smooth_signif=TRUE, outplot=TRUE,
                 ind_plot="best", lowwl=5, highwl="default",
                 congr_brk=5, dirname=NULL, showlastplot=TRUE,
                 save_plot=FALSE)

    p_sprod <- stock_ts %>%
      ggplot(aes(x=year, y=SProd/mean_TB*100))+
      geom_line()+
      geom_point()+
      geom_line(aes(x=X, y=fit), col="blue", data=SProd_traj$best_fit, linewidth=1)+
      scale_y_continuous(
        name = "Normalized productivity\n(% mean total biomass)",
        sec.axis = sec_axis(transform=~.*mean_TB/100,
                            name="Raw productivity"))+
      geom_hline(yintercept = 0, lty="dashed")+
      geom_vline(xintercept =
                   if(SProd_traj$best_traj$class=="abrupt")
                     SProd_traj$res$loc_brk_asd %>%
                   strsplit(";") %>% `[[`(1) %>%
                   as.numeric() %>%
                   suppressWarnings()
                 else as.numeric(NA),
                 col="red", lty="dotted", linewidth=1)+
      geom_vline(xintercept =
                   if(SProd_traj$best_traj$class=="abrupt")
                     as.numeric(SProd_traj$res$loc_brk_chg)
                 else as.numeric(NA),
                 col="blue", lty="dashed", linewidth=1)+
      theme_classic()+
      labs(x="")+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_sprod <- stock_ts %>%
      ggplot(aes(x=year, y=SProd/mean_TB*100))+
      scale_y_continuous(
        name = "Normalized productivity\n(% mean total biomass)",
        sec.axis = sec_axis(transform=~.*mean_TB/100,
                            name="Raw productivity"))+
      labs(x="")+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  }


  # Biomass vs. time:
  if(!all(is.na(stock_ts$TBbest))){

    TBbest_traj <- stock_ts %>%
      dplyr::select(stockid, year, TBbest) %>%
      tidyr::drop_na(TBbest) %>%
      prep_data(thr=NULL, type="data", apriori=FALSE) %>%
      traj_class(sets=., str="aic_asd", abr_mtd=c("chg", "asd"), asd_thr=0.15,
                 type="data", showplots=FALSE,
                 apriori=FALSE, run_loo=FALSE, two_bkps=FALSE,
                 smooth_signif=TRUE, outplot=FALSE,
                 ind_plot=NULL, lowwl=5, highwl="default",
                 congr_brk=5, dirname=NULL,
                 save_plot=FALSE)

    p_b <- stock_ts %>%
      dplyr::mutate(coll_thr = mean(TBbest)*0.25) %>%
      ggplot()+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=coll_thr), fill="pink")+
      scale_y_continuous(expand = c(0, 0))+
      geom_line(aes(x=year, y=TBbest))+
      geom_point(aes(x=year, y=TBbest))+
      labs(y="Biomass", x="")+
      # geom_line(aes(x=X, y=fit), col="blue", data=TBbest_traj$best_fit, linewidth=1)+
      # geom_hline(yintercept = 1, lty="dashed")+
      # geom_vline(xintercept =
      #              if(TBbest_traj$best_traj$class=="abrupt")
      #                TBbest_traj$res$loc_brk_asd %>%
      #              strsplit(";") %>% `[[`(1) %>%
      #              as.numeric() %>%
      #              suppressWarnings()
      #            else as.numeric(NA),
      #            col="red", lty="dotted", linewidth=1)+
      # geom_vline(xintercept =
      #              if(TBbest_traj$best_traj$class=="abrupt")
      #                as.numeric(TBbest_traj$res$loc_brk_chg)
      #            else as.numeric(NA),
      #            col="blue", lty="dashed", linewidth=1)+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_b <- stock_ts %>%
      ggplot(aes(x=year, y=TBbest))+
      labs(y="Biomass", x="")+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  }

  # Combine everything:

  design <- "
  1
  2
"
  fig_stock_summary <-
    p_sprod +
    p_b +
    patchwork::plot_layout(design=design)+
    patchwork::plot_annotation(
      title=paste0(assessment %>%
                     dplyr::filter(stockid==stk) %>%
                     dplyr::slice(1) %>%
                     dplyr::pull(stocklong),
                   " (", stk, ")"),
      theme=theme(plot.title = element_text(hjust = 0.5, size=15)))

  return(fig_stock_summary)
}


stock_mature <- function(stk, drop_SProd=FALSE){

  stock_ts <-

    # Start from SProd and biomass timeseries:
    timeseries_values_views %>%
    dplyr::filter(stockid==stk) %>%
    dplyr::select(stockid, year, SProd, TBbest, TCbest) %>%

    # Add exploitation rate timeseries:
    dplyr::left_join(mature_def_simpl(stk), by="stockid")

  if(drop_SProd==TRUE){
    stock_ts <- stock_ts %>%
      tidyr::drop_na(SProd)
  }

  mean_TB <- mean(stock_ts %>%
                    tidyr::drop_na(SProd) %>%
                    pull(TBbest))

  # Productivity vs. time:
  if(!all(is.na(stock_ts$SProd))){

    p_sprod <- stock_ts %>%
      ggplot(aes(x=year, y=SProd/mean_TB*100))+
      geom_rect(aes(xmin=-Inf, xmax=start_mature, ymin=-Inf, ymax=Inf), fill="yellow", alpha=0.1)+
      geom_rect(aes(xmin=start_mature, xmax=Inf, ymin=-Inf, ymax=Inf), fill="lightblue", alpha=0.1)+
      geom_line()+
      geom_point()+
      scale_y_continuous(
        name = "Normalized productivity\n(% mean total biomass)",
        sec.axis = sec_axis(transform=~.*mean_TB/100,
                            name="Raw productivity"))+
      geom_hline(yintercept = 0, lty="dashed")+
      geom_vline(xintercept =
                   if(SProd_traj$best_traj$class=="abrupt")
                     SProd_traj$res$loc_brk_asd %>%
                   strsplit(";") %>% `[[`(1) %>%
                   as.numeric() %>%
                   suppressWarnings()
                 else as.numeric(NA),
                 col="red", lty="dotted", linewidth=1)+
      geom_vline(xintercept =
                   if(SProd_traj$best_traj$class=="abrupt")
                     as.numeric(SProd_traj$res$loc_brk_chg)
                 else as.numeric(NA),
                 col="blue", lty="dashed", linewidth=1)+
      theme_classic()+
      labs(x="")+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_sprod <- stock_ts %>%
      ggplot(aes(x=year, y=SProd/mean_TB*100))+
      scale_y_continuous(
        name = "Normalized productivity\n(% mean total biomass)",
        sec.axis = sec_axis(transform=~.*mean_TB/100,
                            name="Raw productivity"))+
      labs(x="")+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))
  }

  # Biomass vs. time:
  if(!all(is.na(stock_ts$TBbest))){

    p_b <- stock_ts %>%
      dplyr::mutate(coll_thr = mean(TBbest)*0.25) %>%
      ggplot()+
      geom_rect(aes(xmin=-Inf, xmax=start_mature, ymin=-Inf, ymax=Inf), fill="yellow", alpha=0.1)+
      geom_rect(aes(xmin=start_mature, xmax=Inf, ymin=-Inf, ymax=Inf), fill="lightblue", alpha=0.1)+
      scale_y_continuous(expand = c(0, 0))+
      geom_line(aes(x=year, y=TBbest))+
      geom_point(aes(x=year, y=TBbest))+
      labs(y="Biomass", x="")+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_b <- stock_ts %>%
      ggplot(aes(x=year, y=TBbest))+
      labs(y="Biomass", x="")+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  }


  # Catch vs. time:
  if(!all(is.na(stock_ts$TCbest))){

    p_c <- stock_ts %>%
      ggplot()+
      geom_rect(aes(xmin=-Inf, xmax=start_mature, ymin=-Inf, ymax=Inf), fill="yellow", alpha=0.1)+
      geom_rect(aes(xmin=start_mature, xmax=Inf, ymin=-Inf, ymax=Inf), fill="lightblue", alpha=0.1)+
      scale_y_continuous(expand = c(0, 0))+
      geom_line(aes(x=year, y=TCbest))+
      geom_point(aes(x=year, y=TCbest))+
      labs(y="Catch", x="")+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  } else {

    p_c <- stock_ts %>%
      ggplot(aes(x=year, y=TCbest))+
      labs(y="Biomass", x="")+
      theme_classic()+
      theme(panel.grid=element_blank(),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20))

  }

  # Combine everything:

  design <- "
  1
  2
  3
"
  mature_stock_summary <-
    p_sprod +
    p_b +
    p_c +
    patchwork::plot_layout(design=design)+
    patchwork::plot_annotation(
      title=paste0(assessment %>%
                     dplyr::filter(stockid==stk) %>%
                     dplyr::slice(1) %>%
                     dplyr::pull(stocklong),
                   " (", stk, ")"),
      theme=theme(plot.title = element_text(hjust = 0.5, size=15)))

  return(mature_stock_summary)
}
