###-###-###-###-###-###-###-###-###
#
# Analyse time series classification
#
# 03_analyses.R
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


# Load classified trajectories --------------------------------------------

traj_SProd <-
  readRDS("res/classif/classif_v4.61_SProd_minlen25_normTBavg_looTRUE_aicasd_start1950.rds") %>%
  `[[`("traj_ts_full") %>%
  dplyr::rename(stockid = simu_id) %>%
  dplyr::left_join(as.data.frame(stock) %>%
                     dplyr::select(stockid, scientificname),
                   by="stockid") %>%
  dplyr::left_join(taxonomy %>%
                     dplyr::distinct(scientificname, .keep_all=TRUE) %>%
                     dplyr::select(scientificname, phylum, classname,
                                   ordername, family),
                   by="scientificname") %>%
  dplyr::left_join(metadata %>%
                     dplyr::select(stockid, primary_FAOarea),
                   by="stockid")


dir.create("res/figs", showWarnings=FALSE)

## Fig 1 - Overview of trajectory shapes -----------------------------------

dir.create("res/figs/fig1", showWarnings=FALSE)

st_examples <-
  c("POLLNEAR", # no change
    "STMARLINWCNPAC", # linear decrease
    "SKJCWPAC", # linear increase
    "MORWONGWSE", # quadratic increase accelerated
    "GOLDREDNEAR", # quadratic decrease decelerated
    "HADROCK", # quadratic convex
    "FLSOLEBSAI", # abrupt decrease
    "HAKENRTN" # abrupt increase,
  )

for (i in 1:length(st_examples)){

  print(paste0("Classification in progress for ", st_examples[i]))

  set <- extract_RAM(st_examples[i], "SProd") %>%
    prep_data_simpl()

  mean_TBbest <- extract_RAM(st_examples[i], "TBbest") %>%
    dplyr::filter(year>=min(set$ts[[1]]$X) & year<=max(set$ts[[1]]$X)) %>%
    dplyr::pull(TBbest) %>%
    mean()

  set$ts[[1]]$Y <- set$ts[[1]]$Y/mean_TBbest

  trajs <- traj_class(sets=set, str="aic_asd", abr_mtd=c("chg", "asd"),
                      run_loo=TRUE, two_bkps=TRUE,
                      makeplots=TRUE, ind_plot="best",
                      showlastplot=TRUE)

  pdf(file=paste0("res/figs/fig1/", st_examples[i], ".pdf"), width=4, height=3)
  print(trajs$class_plot$plots[[1]])
  dev.off()
}

# => Panels assembled with Inkscape adding silhouettes and tags



## Fig 2 - Distribution of trajectories and abrupt shifts ----------------

dir.create("res/figs/fig2", showWarnings=FALSE)

### 2a - Overall ---------------------------------------------------------

traj_SProd_sum_det <- traj_SProd %>%
  dplyr::mutate(traj = dplyr::case_when(
    class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
    class=="abrupt" & trend=="increase" ~"increase_abrupt",
    class!="abrupt" ~traj
  )) %>%
  dplyr::group_by(traj, trend, class) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(
    spe_class =
      factor(traj,
             levels = c("stable_constant",
                        "stable_concave","stable_convex",
                        "decrease_constant","decrease_decelerated",
                        "decrease_accelerated","decrease_abrupt",
                        "increase_constant","increase_decelerated",
                        "increase_accelerated","increase_abrupt"))) %>%
  dplyr::ungroup()

traj_SProd_sum_det

# Global distribution of trajectories:
SProd_sunburst <- sunburst_traj_plot_det(traj_SProd_sum_det, size=4,
                                         no_text=FALSE, no_caption=FALSE)+
  labs(tag="A")+
  theme(plot.tag=element_text(size=20, face = "bold"))

pdf(file = "res/figs/fig2/fig2a.pdf",
    width=8, height=6)
print(SProd_sunburst)
dev.off()


### 2bc - Abrupt shifts in time -------------------------------

ts_avail_distr <- timeseries_values_views %>%
  dplyr::filter(stockid %in% (traj_SProd %>% dplyr::pull(stockid)),
                year>=1950) %>%
  dplyr::select(stockid, year, SProd) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(n=n())

ts_avail_distr2 <-
  ts_avail_distr[rep(seq_along(ts_avail_distr$n), ts_avail_distr$n), ]

abr_events <- traj_SProd %>%
  dplyr::filter(class=="abrupt")

(abr_time_pos <-
    ggplot(ts_avail_distr2)+
    geom_histogram(aes(x=year, y=after_stat(count)), binwidth=1,
                   fill="grey70", alpha=0.5, col=NA)+
    geom_density(data = abr_events %>%
                   dplyr::rename(year=loc_brk_chg) %>%
                   dplyr::filter(trend=="increase"),
                 aes(x=year, y=100*after_stat(count)),
                 adjust=1, fill="#000066", alpha=0.3, col=NA)+
    geom_histogram(data = abr_events %>%
                     dplyr::rename(year=loc_brk_chg) %>%
                     dplyr::filter(trend=="increase"),
                   aes(x=year, y=50*after_stat(count)),
                   binwidth=1, fill="#000066", alpha=0.7, col=NA)+

    scale_y_continuous(
      name = "Time series coverage",
      sec.axis = sec_axis(transform=~./50,
                          name="Abrupt increases"),
      position = "right",
      breaks = seq(0, 300, 100),
      expand = c(0, 0))+
    scale_x_continuous(breaks = seq(1950, 2020, 10),
                       expand = c(0, 0))+
    theme_classic()+
    theme(legend.position="none",
          axis.text=element_text(size=10),
          axis.text.x=element_text(vjust=0),
          axis.title=element_text(size=15),
          axis.title.y.left=element_text(colour = "#2a2ab1ff"),
          axis.title.y.right=element_text(colour = "grey50"),
          axis.ticks=element_line(linewidth=.7),
          axis.ticks.length=unit(2, "mm"),
          plot.tag=element_text(size=20, face = "bold"))+
    coord_cartesian(xlim=c(1950, 2020))+
    labs(x="", tag="B"))


(abr_time_neg <-
    ggplot(ts_avail_distr2)+
    geom_histogram(aes(x=year, y=after_stat(count)), binwidth=1,
                   fill="grey70", alpha=0.5, col=NA)+
    geom_density(data = abr_events %>%
                   dplyr::rename(year=loc_brk_chg) %>%
                   dplyr::filter(trend=="decrease"),
                 aes(x=year, y=100*after_stat(count)),
                 adjust=1, fill="#660000", alpha=0.3, col=NA)+
    geom_histogram(data = abr_events %>%
                     dplyr::rename(year=loc_brk_chg) %>%
                     dplyr::filter(trend=="decrease"),
                   aes(x=year, y=50*after_stat(count)),
                   binwidth=1, fill="#660000", alpha=0.7, col=NA)+
    scale_y_continuous(
      name = "Time series coverage",
      sec.axis = sec_axis(transform=~./50,
                          name="Abrupt declines"),
      position = "right",
      breaks = seq(0, 300, 100),
      expand = c(0, 0))+
    scale_x_continuous(breaks = seq(1950, 2020, 10),
                       expand = c(0, 0))+
    theme_classic()+
    theme(legend.position="none",
          axis.text=element_text(size=10),
          axis.text.x=element_text(vjust=0),
          axis.title=element_text(size=15),
          axis.title.y.left=element_text(colour = "#ae2727ff"),
          axis.title.y.right=element_text(colour = "grey50"),
          axis.ticks=element_line(linewidth=.7),
          axis.ticks.length=unit(2, "mm"),
          plot.tag=element_text(size=20, face = "bold"))+
    coord_cartesian(xlim=c(1950, 2020))+
    labs(x="", tag="C"))

abr_time <- abr_time_pos / abr_time_neg

ggsave(plot = abr_time,
       "res/figs/fig2/fig2bc.pdf", width=5, height=5.625)


### 2d - Trajectories in space -------------------------------
areas <-
  sf::read_sf("data/spatial_data/FAO_AREAS_ERASE_LOWRES/FAO_AREAS_ERASE_LOWRES.shp")

FAOarea_major <- areas %>%
  dplyr::filter(F_LEVEL=="MAJOR") %>%
  dplyr::select(F_CODE, geometry) %>%
  ggplot()+
  geom_sf()+
  theme_minimal()+
  theme(panel.grid = element_blank())

sf::st_coordinates(areas %>%
                     dplyr::filter(F_LEVEL=="MAJOR") %>%
                     dplyr::select(F_CODE, geometry))

centroids <- areas %>%
  dplyr::filter(F_LEVEL=="MAJOR") %>%
  dplyr::select(F_CODE, geometry) %>%
  sf::st_centroid() %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  left_join(traj_SProd %>%
              dplyr::group_by(primary_FAOarea) %>%
              dplyr::summarise(sz=n()) %>%
              dplyr::mutate(sz=sqrt(sz+10)*6),
            by=c("F_CODE"="primary_FAOarea"))

traj_SProd_sum_area <- traj_SProd %>%
  split(.$primary_FAOarea) %>%
  mapply(function(x,y){
    x %>%
      dplyr::mutate(traj = dplyr::case_when(
        class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
        class=="abrupt" & trend=="increase" ~"increase_abrupt",
        class!="abrupt" ~traj)) %>%
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
      dplyr::ungroup() %>%
      sunburst_traj_plot_det(.,
                             title=NULL,
                             size_caption=13,
                             size=2, no_text=TRUE, no_caption=FALSE)

  }, x=., y=names(.), SIMPLIFY = FALSE)


traj_SProd_sum_area <- traj_SProd %>%
  split(.$primary_FAOarea) %>%
  mapply(function(x,y){
    x %>%
      dplyr::mutate(traj = dplyr::case_when(
        class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
        class=="abrupt" & trend=="increase" ~"increase_abrupt",
        class!="abrupt" ~traj)) %>%
      dplyr::mutate(
        traj = ifelse(traj=="decrease_constant", "decrease_linear", traj),
        traj = ifelse(traj=="increase_constant", "increase_linear", traj),
        traj = ifelse(traj %in% c("decrease_decelerated", "decrease_accelerated"),
                      "decrease_quadratic", traj),
        traj = ifelse(traj %in% c("increase_decelerated", "increase_accelerated"),
                      "increase_quadratic", traj),
        traj = ifelse(traj %in% c("stable_concave", "stable_convex"),
                      "stable_quadratic", traj)) %>%
      dplyr::group_by(traj, trend, class) %>%
      dplyr::summarise(n=n()) %>%
      dplyr::mutate(
        spe_class =
          factor(traj,
                 levels = c("stable_constant","stable_quadratic",
                            "decrease_linear","decrease_quadratic",
                            "decrease_abrupt","increase_linear",
                            "increase_quadratic","increase_abrupt"))) %>%
      dplyr::ungroup() %>%
      sunburst_traj_plot(.,
                         title=NULL,
                         size_caption=13,
                         size=2, no_text=TRUE, no_caption=FALSE)

  }, x=., y=names(.), SIMPLIFY = FALSE)


centroid_major <-
  tibble(x=centroids$lon,
         y=centroids$lat,
         area=centroids$F_CODE,
         width=centroids$sz) %>%
  dplyr::right_join(tibble(area = names(traj_SProd_sum_area),
                           sunburst=traj_SProd_sum_area),
                    by = "area")

space <- ggplot(areas %>%
                  dplyr::filter(F_LEVEL=="MAJOR") %>%
                  dplyr::select(F_CODE, geometry))+
  geom_sf()+
  theme_void()+
  coord_sf(ylim = c(90, -60), expand = FALSE)+
  ggimage::geom_subview(aes(x=x, y=y, subview=sunburst,
                            width=width, height=width),
                        data=centroid_major)+
  theme(plot.tag=element_text(size=20, face = "bold"))+
  labs(tag="D")

pdf(file = "res/figs/fig2/fig2d.pdf",
    width=12, height=9)
print(space)
dev.off()


### 2e - Trajectory classification by taxonomic order -----------------------

# Number of stocks in RAMLDB:
sp_RAM <- timeseries_values_views %>%
  dplyr::filter(!is.na(TBbest)) %>%
  dplyr::group_by(stockid) %>%
  dplyr::slice(1) %>%
  dplyr::left_join(as.data.frame(stock) %>%
                     dplyr::select(stockid, scientificname), by="stockid") %>%
  dplyr::left_join(taxonomy %>%
                     dplyr::group_by(scientificname) %>%
                     dplyr::slice(1), by="scientificname") %>%
  # Excluding non-fish and Salmonids
  dplyr::filter(phylum == "Chordata" & ordername != "Salmoniformes") %>%
  dplyr::group_by(ordername) %>%
  dplyr::summarise(ram_stk=n())

traj_RAM <- traj_SProd %>%
  dplyr::mutate(traj_simpl = dplyr::case_when(
    traj %in% c("stable_constant", "stable_convex", "stable_concave") ~ "stable",
    traj %in% c("increase_constant", "increase_decelerated", "increase_accelerated") ~ "increase gradual",
    traj == "1_breakpoint" & trend=="increase" ~ "positive PAS",
    traj %in% c("decrease_constant", "decrease_decelerated") ~ "decrease gradual",
    traj == "1_breakpoint" & trend=="decrease" ~ "negative PAS"),
    traj_simpl = factor(traj_simpl,
                        levels=c("stable", "increase gradual", "positive PAS",
                                 "decrease gradual", "negative PAS"))) %>%
  dplyr::group_by(traj_simpl, ordername) %>%
  dplyr::summarise(stk_traj=n()) %>%
  dplyr::group_by(ordername) %>%
  dplyr::mutate(n_tot = sum(stk_traj)) %>%
  dplyr::right_join(sp_RAM, by="ordername")

(p_ram_traj <- ggplot(traj_RAM)+
    geom_col(aes(x=reorder(ordername, ram_stk), y=ram_stk), fill="grey",
             position="identity", width=0.7)+
    geom_col(aes(x=reorder(ordername, ram_stk), y=stk_traj, fill=traj_simpl),
             position="stack", width=0.7)+
    scale_fill_manual(
      values =
        c("stable"="#E8DE9C",
          "decrease gradual"="#F08080", "negative PAS"="#660000",
          "increase gradual"="#99CCFF", "positive PAS"="#000066"))+
    labs(y="Number of stocks", fill="trajectory")+
    labs(x="")+
    theme_bw()+
    coord_flip()+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    theme(panel.grid.major.y = element_blank(),
          axis.text=element_text(size=15),
          axis.title = element_text(size=15),
          legend.position="bottom"))


# FAO catch 1950-2022
fish_file <- readr::read_csv("data/FAO_data/Capture_2024.1.0/Capture_Quantity.csv")
species <- readr::read_csv("data/FAO_data/Capture_2024.1.0/CL_FI_SPECIES_GROUPS.csv")
area <- readr::read_csv("data/FAO_data/Capture_2024.1.0/CL_FI_WATERAREA_GROUPS.csv")
taxo <- readr::read_delim("data/FAO_data/ASFIS_sp_2023/ASFIS_sp_2023.txt")

FAOcatch_fish <- fish_file %>%
  dplyr::left_join(species %>%
                     dplyr::select(`3A_Code`, Scientific_Name,
                                   ISSCAAP_Group_En, Major_Group),
                   by=c("SPECIES.ALPHA_3_CODE"="3A_Code")) %>%
  dplyr::left_join(area %>%
                     dplyr::select(Code, Name_En) %>%
                     dplyr::rename(area=Name_En),
                   by=c("AREA.CODE"="Code")) %>%
  dplyr::filter(!grepl("Inland", area) & Major_Group=="PISCES") %>%
  dplyr::group_by(ISSCAAP_Group_En, Scientific_Name, SPECIES.ALPHA_3_CODE) %>%
  dplyr::summarise(catch = sum(VALUE), .groups = "keep") %>%
  dplyr::arrange(dplyr::desc(catch)) %>%
  dplyr::left_join(taxo %>%
                     dplyr::mutate(Order =
                                     ifelse(Family %in%
                                              c("OSMERIDAE", "PLECOGLOSSIDAE",
                                                "RETROPINNIDAE", "SALANGIDAE"),
                                            "Osmeriformes", Order)) %>%
                     dplyr::select(Alpha3_Code, Order),
                   by=c("SPECIES.ALPHA_3_CODE" = "Alpha3_Code")) %>%
  dplyr::mutate(Order = stringr::str_to_title(Order),
                Order = dplyr::recode(Order,
                                      Scombroidei = "Perciformes",
                                      Percoidei = "Perciformes",
                                      Gobioidei = "Perciformes",
                                      Zoarcoidei = "Perciformes",
                                      "Stromateoidei, Anabantoidei" = "Perciformes",
                                      "Other Perciformes" = "Perciformes",
                                      Trachinoidei = "Trachiniformes",
                                      Acanthuroidei = "Acanthuriformes",
                                      .default = Order)) %>%
  dplyr::filter(Order != "Pisces Miscellanea") %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(totcatch = sum(catch)) %>%
  dplyr::rename(ordername = Order) %>%
  dplyr::filter(ordername %in% sp_RAM$ordername) %>%
  arrange(totcatch) %>%
  dplyr::right_join(sp_RAM, by="ordername")

(p_fao_catch <- ggplot(FAOcatch_fish)+
    geom_col(aes(x=reorder(ordername, ram_stk), y=totcatch/1000000),
             fill="grey40", position="identity", width=0.7)+
    labs(y="Total catch 1950-2022 (million tons)")+
    labs(x="")+
    theme_bw()+
    coord_flip()+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(size=15),
          axis.title = element_text(size=15)))

pdf(file = paste0("res/figs/fig2/fig_2e_taxord.pdf"),
    width=12, height=5)
print(p_ram_traj + p_fao_catch)
dev.off()

# => Panels assembled with Inkscape to arrange the legends


### Statistical tests ----------------------------------------------------

# Distribution in time
set.seed(1)
ks.test(ts_avail_distr2 %>%
          dplyr::pull(year),
        abr_events %>%
          dplyr::filter(trend=="decrease") %>%
          dplyr::pull(loc_brk_chg),
        simulate.p.value = TRUE, B = 1e5)
# p-value = 0.0382

set.seed(1)
ks.test(ts_avail_distr2 %>%
          dplyr::pull(year),
        abr_events %>%
          dplyr::filter(trend=="increase") %>%
          dplyr::pull(loc_brk_chg),
        simulate.p.value = TRUE, B = 1e5)
# p-value = 0.3179


# By class
class_mat_fao <- traj_SProd %>%
  dplyr::group_by(primary_FAOarea, class) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::pivot_wider(names_from=class, values_from=n) %>%
  tibble::column_to_rownames("primary_FAOarea") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(class_mat_fao, main="class by FAO area", color = TRUE)

set.seed(1)
(chi_clas_fao <- chisq.test(class_mat_fao, simulate.p.value = TRUE, B=1e5))
# p-value = 0.3227

# class_chi_bal <- apply(class_mat_fao, 1,
#       function(i) chisq.test(x=i, p=c(0.25,0.25,0.25,0.25),
#                              simulate.p.value=TRUE, B=1e5)$p.value)


# Taxonomic order by trajectory (abrupt decrease, abrupt increase, nonabrupt)
traj_SProd_sum_ord <- traj_SProd %>%
  dplyr::mutate(traj = dplyr::case_when(
    traj == "1_breakpoint" & trend=="increase" ~ "positive PAS",
    traj == "1_breakpoint" & trend=="decrease" ~ "negative PAS",
    .default = "non abrupt")) %>%
  dplyr::group_by(traj, ordername) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup()

# Keep only the five more numerous orders in number of stocks included
traj_SProd_sum_ord <- traj_SProd_sum_ord %>%
  dplyr::filter(ordername %in%
                  c("Perciformes", "Gadiformes", "Scorpaeniformes",
                    "Pleuronectiformes", "Clupeiformes"))

traj_mat_ord <- traj_SProd_sum_ord %>%
  tidyr::pivot_wider(names_from=traj, values_from=n) %>%
  tibble::column_to_rownames("ordername") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(traj_mat_ord, main="class by taxonomic order", color = TRUE)

set.seed(1)
(chi_traj_ord <- chisq.test(traj_mat_ord, simulate.p.value=TRUE, B=1e5))
# p-value = 0.01203

chi_traj_ord$stdres %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="taxclass") %>%
  tidyr::pivot_longer(cols = -taxclass, names_to = "traj", values_to="stdres") %>%
  dplyr::arrange(desc(abs(stdres)))



## Fig 3 - Explanatory factors of PAS ----------------------------

dir.create("res/figs/fig3", showWarnings=FALSE)

# Add traits to stocks:
traits_default <- traits_fun(traj_SProd)

readr::write_csv(traits_default, "res/classif/traj_SProd_classif_traits.csv")
traits_default <- readr::read_csv("res/classif/traj_SProd_classif_traits.csv")

mod_df01_default <-
  traj_SProd %>%
  dplyr::mutate(length = last-first+1,
                length_prshf = ifelse(class=="abrupt", loc_brk_chg-first+1, length)) %>%
  dplyr::select(stockid, scientificname, class, trend,
                length, length_prshf, abruptness, rel_chg, weight_aic_abrupt) %>%
  dplyr::mutate(shift = ifelse(class=="abrupt", 1, 0),
                shift_inc = ifelse(class=="abrupt" & trend=="increase", 1, 0),
                shift_dec = ifelse(class=="abrupt" & trend=="decrease", 1, 0),
                shift_both = dplyr::case_when(shift_inc==1 ~"abr_inc",
                                              shift_dec==1 ~"abr_dec",
                                              shift_inc==0 & shift_dec==0 ~"no_abr")) %>%
  dplyr::rename(Species=scientificname) %>%
  dplyr::left_join(traits_default,
                   by=c("stockid", "Species")) %>%
  dplyr::mutate(DemersPelag = factor(DemersPelag, levels=c("pelagic", "demersal")))


### 3a - Abrupt declines -----

# Prepare dataset and scale variables
df_scl <- mod_df01_default %>%
  tidyr::drop_na(tm, Lm, DemersPelag, ordername, family, Species, X, Y,
                 sst_avg, sst_change, mean_ER_prshf, ER_change_prshf) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric) &
                                !dplyr::contains(c("shift", "X", "Y")),
                              ~ c(scale(.)))) %>%
  dplyr::mutate(Species = as.factor(Species),
                family = as.factor(family),
                ordername = as.factor(ordername))

# Run GAMM
modgam_dec <- mgcv::gam(
  formula=shift_dec ~ tm + Lm +
    DemersPelag +
    sst_avg_prshf + sst_change_prshf +
    mean_ER_prshf + ER_change_prshf +
    sst_avg_prshf:mean_ER_prshf +
    sst_change_prshf:ER_change_prshf +
    s(X, Y, bs="tp") +
    s(ordername, bs="re") + s(family, bs="re"),
  family=binomial(link='logit'),
  na.action=na.omit,
  data=df_scl)

summary(modgam_dec)
table(df_scl$shift_dec)

# Plot model estimates
(m_dec <- ggstats::ggcoef_model(
  modgam_dec,
  colour = NULL, stripped_rows = FALSE,
  conf.level = 0.95,
  significance = 0.05,
  add_reference_rows=TRUE,
  point_size = 3,
  point_stroke = 0.5,
  x = "estimate"))

pdf(file = "res/figs/fig3/fig3a_dec.pdf",
    width=7, height=4.5)
print(m_dec)
dev.off()


### 3b - Abrupt increases -----

# Run GAMM
modgam_inc <- mgcv::gam(
  formula=shift_inc ~ tm + Lm +
    DemersPelag+
    sst_avg_prshf + sst_change_prshf +
    mean_ER_prshf + ER_change_prshf +
    sst_avg_prshf:mean_ER_prshf +
    sst_change_prshf:ER_change_prshf +
    s(X, Y, bs="tp") +
    s(ordername, bs="re") + s(family, bs="re"),
  family=binomial(link='logit'),
  na.action=na.omit,
  data=df_scl)

summary(modgam_inc)
table(df_scl$shift_inc)

# Plot model estimates
(m_inc <- ggstats::ggcoef_model(
  modgam_inc,
  colour = NULL, stripped_rows = FALSE,
  conf.level = 0.95,
  significance = 0.05,
  add_reference_rows=TRUE,
  point_size = 3,
  point_stroke = 0.5,
  x = "estimate"))

pdf(file = "res/figs/fig3/fig3b_inc.pdf",
    width=7, height=4.5)
print(m_inc)
dev.off()


# => Panels assembled with Inkscape to arrange the legends


## Fig 4 - PAS collapse ---------------------------------------

dir.create("res/figs/fig4", showWarnings=FALSE)

fig4 <- shift_coll_plot(traj_SProd, coll_def=c("bavg", 0.25), csec_coll = 2)
# csec_coll=2 to avoid artifacts

ggsave(filename="res/figs/fig4/fig4_bavg0.25_cseccoll2.pdf",
       width=14, height=8, plot=fig4)

# => Improved rendering in Inkscape


### Statistical tests ----------------------------------------------------

# Magnitude of abrupt declines
coll <- add_collapsed(traj_SProd, coll_def=c("bavg", 0.25), csec_coll=2)

(t_dec <-
    t.test(coll[coll$traj=="decrease abrupt" & coll$did_collapse=="Yes",]$mag,
           coll[coll$traj=="decrease abrupt" & coll$did_collapse=="No",]$mag))
# p-value = 0.05855

# Test for abrupt increase magnitude not possible since only one stock collapsed

# Stocks with abrupt decrease followed by collapse
coll %>% dplyr::filter(traj=="decrease abrupt" & did_collapse=="Yes") %>%
  dplyr::mutate(diff = loc_brk_chg-first_coll_y) %>%
  dplyr::arrange(diff)


