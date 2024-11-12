###-###-###-###-###-###-###-###-###
#
# Functions for joint pattern analyses
#
# functions_jointpattern.R
#
###-###-###-###-###-###-###-###-###


# Data preparation and selection -----------------------------------------

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



#' Extract stocks RAMLBD timeseries
#' (the RAMLBD should be loaded in the environment)
#'
#' @param id Short stock ID character.
#' @param ts_type Type of timeseries (TBbest, ERbest...) character.
#' @param drop_na Keep only years with all types of timeseries available.
#'
#' @return Data frame of the selected timeseries.
#'
#' @export

extract_RAM <- function(id, ts_type, drop_na=TRUE){

  ts <- timeseries_values_views %>%
    dplyr::filter(stockid %in% id) %>%
    dplyr::mutate(scen = paste(ts_type, stockid, sep="_")) %>%
    dplyr::select(scen, year, dplyr::all_of(ts_type)) %>%
    dplyr::relocate(scen)

  if (drop_na){
    ts <- ts %>% na.omit()
  }

  return(ts)

}



# Traits -----------------------------------------------

#' Extract all stock traits
#'
#' @param traj_SProd Data frame with sums for trend and class
#'
#' @return a data frame with traits as columns and as many rows as stocks
#' @export

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



# Plots ---------------------------------------------------------


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



#' Plot relationship between warming rate and a given proportion at LME level
#'
#' @param df
#' @param main_lme
#' @param asstchange_LME
#' @param filter
#' @param weight
#' @param min_tot
#'
#' @return
#' @export
#'
#' @examples

prop_warm_plot <- function(df, main_lme, asstchange_LME, filter,
                           weight=NULL, min_tot=0){

  prop_sstchange_by_lme <- dplyr::left_join(df, main_lme, by="stockid") %>%
    dplyr::group_by(LME_NAME) %>%
    dplyr::summarise(n_tot = n()) %>%
    dplyr::left_join(
      dplyr::left_join(df, main_lme, by="stockid") %>%
        dplyr::filter(eval(filter)) %>%
        dplyr::group_by(LME_NAME) %>%
        dplyr::summarise(count = n()),
      by="LME_NAME") %>%
    dplyr::mutate(prop = count/n_tot) %>%
    left_join(asstchange_LME, by=c("LME_NAME"="lme"))

  if (weight=="n_tot"){
    mod <- lm(prop~sst_change,
              weights=n_tot,
              prop_sstchange_by_lme %>%
                dplyr::filter(n_tot>min_tot))
  } else {
    mod <- lm(prop~sst_change,
              prop_sstchange_by_lme %>%
                dplyr::filter(n_tot>min_tot))
  }

  summ <- summary(mod)
  coeff <- paste0("p = ", format(round(summ$coefficients[2,4],3), nsmall=3),
                  ", R2 = ", round(summ$adj.r.squared, 2))

  if (summ$coefficients[2,4] < 0.05){
    line <- "solid"
  } else {
    line <- "dashed"
  }

  plot <- prop_sstchange_by_lme %>%
    dplyr::filter(n_tot>min_tot) %>%
    tidyr::drop_na(prop) %>%
    ggplot(aes(x=sst_change, y=prop))+
    geom_point(aes(size=n_tot), alpha=0.6)+
    geom_abline(intercept=summ$coefficients[1,1],
                slope=summ$coefficients[2,1],
                col="grey30", lty=line, linewidth=1.5)+
    geom_text(x = 0.007, y = 0.9, label = coeff, size=5, check_overlap = TRUE)+
    coord_cartesian(ylim=c(0,1))+
    scale_y_continuous(expand=c(0,0))+
    theme_classic()+
    guides(size="none")+
    labs(x="SST rate of change (1950-2020) [°C/y]")

  return(plot)
}



#' Render plots combining abrupt shifts and collapse
#'
#' @param traj_SProd A classification summary.
#' @param coll_def A vector defining the threshold of collapsed state.
#' The 1st value indicates the parameter considered ("bavg", "bmax", or "bmsy",
#' for average, maximum biomass, and biomass at MSY respectively). The second
#' value indicates the fraction of biomass parameter considered as threshold.
#' @param csec_coll An integer specifying the minimal number of consecutive
#' years of under collapse threshold to be considered as collapsed.
#'
#' @return
#' @export
#'
#' @examples

# coll_def=c("bavg", 0.25) # Essington et al. 2015
# coll_def=c("bmsy", 0.2) # Costello et al. 2012

shift_coll_plot <- function(traj_SProd, coll_def=c("bavg", 0.25), csec_coll=1){

  # Add info about collapse:
  coll <- add_collapsed(traj_SProd, coll_def=coll_def, csec_coll=csec_coll)

  # Prepare data frame for plot:
  prop_coll_shift_df <- coll %>%
    dplyr::mutate(traj_simpl = dplyr::case_when(
      traj %in% c("stable constant", "stable convex", "stable concave") ~ "stable",
      traj %in% c("increase linear", "increase decelerated", "increase accelerated") ~ "increase gradual",
      traj %in% c("increase abrupt") ~ "positive PAS",
      traj %in% c("decrease linear", "decrease decelerated") ~ "decrease gradual",
      traj %in% c("decrease abrupt") ~ "negative PAS"),
      traj_simpl = factor(traj_simpl,
                          levels=c("stable", "increase gradual", "positive PAS",
                                   "decrease gradual", "negative PAS"))) %>%
    dplyr::group_by(traj_simpl, did_collapse, coll_def) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::group_by(did_collapse) %>%
    dplyr::mutate(prop = count/sum(count))

  # Proportion plot ---
  (prop_coll_shift <- prop_coll_shift_df %>%
     ggplot()+
     geom_col(aes(x=did_collapse, y=count, fill=traj_simpl),
              position="fill", width=0.5)+
     geom_text(aes(x=did_collapse, y=count, label=count, group=traj_simpl,
                   col=grepl("PAS", traj_simpl)),
               size = 4, position = position_fill(vjust = 0.5))+
     scale_y_continuous(labels = scales::percent, expand = c(0,0))+
     scale_fill_manual(
       values =
         c("stable"="#E8DE9C",
           "decrease gradual"="#F08080", "negative PAS"="#660000",
           "increase gradual"="#99CCFF", "positive PAS"="#000066"))+
     scale_color_manual(values = c("TRUE"="white", "FALSE"="black"))+
     labs(y="frequency", fill="Trajectory")+
     labs(x="Did stock collapsed?")+
     theme_classic()+
     labs(tag="A")+
     theme(axis.text = element_text(size=12),
           axis.text.x = element_text(size=15),
           axis.title = element_text(size=15),
           plot.tag=element_text(size=20, face="bold"))+
     guides(color="none", fill="none"))


  # Abruptness plot ---
  (abr_coll_shift <- coll %>%

     dplyr::filter(class=="abrupt") %>%
     ggplot()+
     geom_boxplot(aes(x=did_collapse, y=mag, fill=trend), col="grey30"
                  # , outlier.shape = NA
     )+
     scale_fill_manual(values =
                         c("decrease"="#660000", "increase"="#000066"))+
     labs(y="shift magnitude")+
     labs(x="Did stock collapsed?")+
     coord_cartesian(y=c(-0.8, 0.8))+
     geom_hline(yintercept = 0, lty=1)+
     theme_classic()+
     labs(tag="B")+
     theme(axis.text = element_text(size=12),
           axis.text.x = element_text(size=15),
           axis.title = element_text(size=15),
           plot.tag=element_text(size=20, face="bold"))+
     guides(fill="none"))


  # Timing of abrupt shifts relative to collapse ---
  (time_coll_shift <- coll %>%
     dplyr::mutate(diff = loc_brk_chg-first_coll_y,
                   diff = ifelse(class == "abrupt", diff, NA)
     ) %>%
     dplyr::filter(class=="abrupt") %>%
     tidyr::drop_na(diff) %>%
     ggplot()+
     geom_histogram(aes(x=diff, fill=trend), col="grey30", binwidth = 2.5, boundary=0)+
     geom_vline(xintercept = 0, lty="dotted")+
     scale_fill_manual(values =
                         c("decrease"="#660000", "increase"="#000066"))+

     labs(x="year since collapse", y="number of stocks")+
     facet_wrap(vars(traj), nrow = 2)+
     theme_bw()+
     scale_y_continuous(expand=c(0,0))+
     expand_limits(y=c(0,4))+
     scale_x_continuous(breaks = seq(-50, 50, 5))+
     theme(strip.text = element_text(colour = "black"))+
     guides(fill="none")+
     labs(tag="C")+
     theme(axis.text = element_text(size=12),
           axis.title = element_text(size=15),
           # panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           strip.text = element_text(size=15),
           plot.tag=element_text(size=20, face="bold")))



  # Correlation proportion vs. warming ---

  atemp_LME <- readr::read_csv("data/spatial_data/lme_yearly_had_sst.csv")

  asstchange_LME <- atemp_LME %>%
    dplyr::filter(year>=1950) %>%
    dplyr::group_by(lme) %>%
    dplyr::summarise(sst_change=summary(lm(sst_c~year))$coefficients["year","Estimate"])

  # Load main LME
  main_lme <- readr::read_csv("data/spatial_data/main_lme.csv")

  # Correlation negative PAS vs. warming ---
  (nPAS_warm <- prop_warm_plot(traj_SProd, main_lme, asstchange_LME,
                               expression(class=="abrupt" & trend=="decrease"),
                               weight="", 1)+
     labs(y="Proportion of negative PAS",
          tag="D")+
     theme(axis.text = element_text(size=12),
           axis.title = element_text(size=15),
           axis.title.x = element_text(size=13),
           plot.tag=element_text(size=20, face="bold")))


  # Correlation positive PAS vs. warming ---
  (pPAS_warm <- prop_warm_plot(traj_SProd, main_lme, asstchange_LME,
                               expression(class=="abrupt" & trend=="increase"),
                               weight="", 1)+
     labs(y="Proportion of positive PAS",
          tag="E")+
     theme(axis.text = element_text(size=12),
           axis.title = element_text(size=15),
           axis.title.x = element_text(size=13),
           plot.tag=element_text(size=20, face="bold")))

  # Correlation collapse vs. warming ---

  (coll_warm <- prop_warm_plot(coll, main_lme, asstchange_LME,
                               expression(did_collapse=="Yes"),
                               weight="", 1)+
     labs(y="Proportion of collapses",
          tag="F")+
     theme(axis.text = element_text(size=12),
           axis.title = element_text(size=15),
           axis.title.x = element_text(size=13),
           plot.tag=element_text(size=20, face="bold")))

  comp_fig <- (prop_coll_shift + abr_coll_shift + time_coll_shift) /
    (nPAS_warm + pPAS_warm + coll_warm)

  return(comp_fig)

}
