###-###-###-###-###-###-###-###-###
#
# Analyse time series classification
#
# 03_analyses.R
#
###-###-###-###-###-###-###-###-###

# Load packages, functions, and data --------------------------------------
source("00_Rpackages.R")
source("R/functions_classification.R")
source("R/functions_analyses.R")

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
  dplyr::mutate(spe_class =
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
SProd_sunburst
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
      dplyr::mutate(spe_class =
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


centroid_major <- tibble(x=centroids$lon,
                         y=centroids$lat,
                         area=centroids$F_CODE,
                         width=centroids$sz
) %>%
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
space
dev.off()

# => Panels assembled with Inkscape to arrange the legends


### Statistical tests ----------------------------------------------------

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

# By trajectory (abrupt decrease, abrupt increase, nonabrupt)
traj_SProd_sum_fao <- traj_SProd %>%
  dplyr::mutate(traj = dplyr::case_when(
    class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
    class=="abrupt" & trend=="increase" ~"increase_abrupt",
    class!="abrupt" ~"nonabrupt")) %>%
  dplyr::group_by(traj, primary_FAOarea) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup()

traj_mat_fao <- traj_SProd_sum_fao %>%
  tidyr::pivot_wider(names_from=traj, values_from=n) %>%
  tibble::column_to_rownames("primary_FAOarea") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(traj_mat_fao, main="class by FAO area", color = TRUE)

set.seed(1)
(chi_traj_fao <- chisq.test(traj_mat_fao, simulate.p.value=TRUE, B=1e5))
# p-value = 0.01647

chi_traj_fao$stdres %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="FAOarea") %>%
  tidyr::pivot_longer(cols = contains("abrupt"), names_to = "traj", values_to="stdres") %>%
  dplyr::arrange(desc(abs(stdres)))

# fisher.test(traj_mat_fao, simulate.p.value=TRUE, B=1e5)


# Species by class
class_mat_ord <- traj_SProd %>%
  dplyr::group_by(ordername, class) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::pivot_wider(names_from=class, values_from=n) %>%
  tibble::column_to_rownames("ordername") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(class_mat_ord, main="class by taxonomic order", color = TRUE)

set.seed(1)
(chi_clas_ord <- chisq.test(class_mat_ord, simulate.p.value = TRUE, B=1e5))
# p-value = 0.09844


# Species by trajectory (abrupt decrease, abrupt increase, nonabrupt)
traj_SProd_sum_ord <- traj_SProd %>%
  dplyr::mutate(traj = dplyr::case_when(
    class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
    class=="abrupt" & trend=="increase" ~"increase_abrupt",
    class!="abrupt" ~"nonabrupt")) %>%
  dplyr::group_by(traj, ordername) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup()

traj_mat_ord <- traj_SProd_sum_ord %>%
  tidyr::pivot_wider(names_from=traj, values_from=n) %>%
  tibble::column_to_rownames("ordername") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(traj_mat_ord, main="class by taxonomic order", color = TRUE)

set.seed(1)
(chi_traj_ord <- chisq.test(traj_mat_ord, simulate.p.value=TRUE, B=1e5))
# p-value = 0.1743


# Species by class
class_mat_fam <- traj_SProd %>%
  dplyr::group_by(family, class) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::pivot_wider(names_from=class, values_from=n) %>%
  tibble::column_to_rownames("family") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(class_mat_fam, main="class by taxonomic family", color = TRUE)

set.seed(1)
(chi_clas_fam <- chisq.test(class_mat_fam, simulate.p.value = TRUE, B=1e5))
# p-value = 0.09844


# Species by trajectory (abrupt decrease, abrupt increase, nonabrupt)
traj_SProd_sum_fam <- traj_SProd %>%
  dplyr::mutate(traj = dplyr::case_when(
    class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
    class=="abrupt" & trend=="increase" ~"increase_abrupt",
    class!="abrupt" ~"nonabrupt")) %>%
  dplyr::group_by(traj, family) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup()

traj_mat_fam <- traj_SProd_sum_fam %>%
  tidyr::pivot_wider(names_from=traj, values_from=n) %>%
  tibble::column_to_rownames("family") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(traj_mat_fam, main="class by taxonomic family", color = TRUE)

set.seed(1)
(chi_traj_fam <- chisq.test(traj_mat_fam, simulate.p.value=TRUE, B=1e5))
# p-value = 0.1743


## Fig 3 - Explanatory factors of PAS ----------------------------

dir.create("res/figs/fig3", showWarnings=FALSE)

# Add traits to stocks:
traits_default <- traits_fun(traj_SProd)

# readr::write_csv(traits_default, "res/classif/traj_SProd_classif_traits.csv")
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
modgam_dec <- mgcv::gamm(
  formula=shift_dec ~ tm + Lm +
    DemersPelag +
    sst_avg_prshf + sst_change_prshf +
    # sst_avg + sst_change +
    mean_ER_prshf + ER_change_prshf +
    s(X, Y) + s(length, bs="re") +
    s(ordername, bs="re") + s(family, bs="re") + s(Species, bs="re"),
  family=binomial(link='logit'),
  na.action=na.omit,
  data=df_scl)

summary(modgam_dec$gam)

# Plot model estimates
(m_dec <- ggstats::ggcoef_model(modgam_dec$gam,
                                colour = NULL, stripped_rows = FALSE,
                                conf.level = 0.95,
                                significance = 0.05,
                                add_reference_rows=TRUE,
                                point_size = 3,
                                point_stroke = 0.5,
                                x = "estimate"))

pdf(file = "res/figs/fig3/fig3a_dec.pdf",
    width=6, height=4)
m_dec
dev.off()

table(df_scl$shift_dec)


### 3b - Abrupt increases -----

# Run GAMM
modgam_inc <- mgcv::gamm(
  formula=shift_inc ~ tm + Lm +
    DemersPelag+
    sst_avg_prshf + sst_change_prshf +
    mean_ER_prshf + ER_change_prshf +
    s(X, Y) + s(length, bs="re") +
    s(ordername, bs="re") + s(family, bs="re") + s(Species, bs="re"),
  family=binomial(link='logit'),
  na.action=na.omit,
  data=df_scl)

summary(modgam_inc$gam)

# Plot model estimates
(m_inc <- ggstats::ggcoef_model(modgam_inc$gam,
                                colour = NULL, stripped_rows = FALSE,
                                conf.level = 0.95,
                                significance = 0.05,
                                add_reference_rows=TRUE,
                                point_size = 3,
                                point_stroke = 0.5,
                                x = "estimate"))

pdf(file = "res/figs/fig3/fig3b_inc.pdf",
    width=6, height=4)
m_inc
dev.off()

table(df_scl$shift_inc)

## Fig 4  ------------------------------------------------------------------


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

  # 4a - Proportion plot ----
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


  # 4b - Abruptness plot ----
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


  # 4c - Timing of abrupt shifts relative to collapse ----
  (time_coll_shift <- coll %>%
     dplyr::mutate(diff = loc_brk_chg-first_coll_y,
                   diff = ifelse(class == "abrupt", diff, NA)
     ) %>%
     dplyr::filter(class=="abrupt") %>%
     tidyr::drop_na(diff) %>%
     ggplot()+
     geom_histogram(aes(x=diff, fill=trend), col="grey30", binwidth = 2.5, origin=0)+
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



  # Correlation proportion vs. warming --------------------------

  atemp_LME <- readr::read_csv("data/spatial_data/lme_yearly_had_sst.csv")

  asstchange_LME <- atemp_LME %>%
    dplyr::filter(year>=1950) %>%
    dplyr::group_by(lme) %>%
    dplyr::summarise(sst_change=summary(lm(sst_c~year))$coefficients["year","Estimate"])

  # Define main LME
  frac_lme <- readr::read_csv("data/spatial_data/frac_lme.csv")
  main_lme <- frac_lme %>%
    dplyr::filter(!grepl("Ocean", LME_NAME)) %>%
    dplyr::group_by(stockid) %>%
    dplyr::slice_max(prop_area)

  # Specify high-sea "LME" manually
  {
    main_lme[main_lme$stockid=="ALBANATL","LME_NAME"] <- "North Atlantic Ocean"
    main_lme[main_lme$stockid=="BIGEYECWPAC","LME_NAME"] <- "North Pacific Ocean"
    main_lme[main_lme$stockid=="BLKMARLINIO","LME_NAME"] <- "Indian Ocean"
    main_lme[main_lme$stockid=="BLSHARIO","LME_NAME"] <- "Indian Ocean"
    main_lme[main_lme$stockid=="BMARLINATL","LME_NAME"] <- "South Atlantic Ocean"
    main_lme[main_lme$stockid=="BMARLINIO","LME_NAME"] <- "Indian Ocean"
    main_lme[main_lme$stockid=="PACBTUNA","LME_NAME"] <- "South Pacific Ocean"
    main_lme[main_lme$stockid=="SAILEATL","LME_NAME"] <- "South Atlantic Ocean"
    main_lme[main_lme$stockid=="SAILWATL","LME_NAME"] <- "South Atlantic Ocean"
    main_lme[main_lme$stockid=="SBT","LME_NAME"] <- "Indian Ocean"
    main_lme[main_lme$stockid=="SKJCIO","LME_NAME"] <- "Indian Ocean"
    main_lme[main_lme$stockid=="SKJCWPAC","LME_NAME"] <- "North Pacific Ocean"
    main_lme[main_lme$stockid=="SKJEATL","LME_NAME"] <- "South Atlantic Ocean"
    main_lme[main_lme$stockid=="SKJWATL","LME_NAME"] <- "North Atlantic Ocean"
    main_lme[main_lme$stockid=="STMARLINIO","LME_NAME"] <- "Indian Ocean"
    main_lme[main_lme$stockid=="STMARLINNEPAC","LME_NAME"] <- "North Pacific Ocean"
    main_lme[main_lme$stockid=="STMARLINSWPO","LME_NAME"] <- "South Pacific Ocean"
    main_lme[main_lme$stockid=="STMARLINWCNPAC","LME_NAME"] <- "North Pacific Ocean"
    main_lme[main_lme$stockid=="SWORDEPAC","LME_NAME"] <- "South Pacific Ocean"
    main_lme[main_lme$stockid=="SWORDNATL","LME_NAME"] <- "North Atlantic Ocean"
    main_lme[main_lme$stockid=="SWORDNPAC","LME_NAME"] <- "North Pacific Ocean"
    main_lme[main_lme$stockid=="WMARLINATL","LME_NAME"] <- "South Pacific Ocean"
    main_lme[main_lme$stockid=="YFINCWPAC","LME_NAME"] <- "North Pacific Ocean"
    }

  prop_warm_plot <- function(df, main_lme, filter, weight=NULL, min_tot=0){

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
      geom_point(aes(size=n_tot))+
      geom_abline(intercept=summ$coefficients[1,1],
                  slope=summ$coefficients[2,1],
                  col="grey30", lty=line)+
      geom_text(x = 0.007, y = 0.9, label = coeff, size=5, check_overlap = TRUE)+
      coord_cartesian(ylim=c(0,1))+
      scale_y_continuous(expand=c(0,0))+
      theme_classic()+
      guides(size="none")+
      labs(x="SST rate of change (1950-2020) [°C/y]")

    return(plot)
  }

  # 4d - Correlation negative PAS vs. warming --------------------------
  (fig4d <- prop_warm_plot(traj_SProd, main_lme,
                           expression(class=="abrupt" & trend=="decrease"),
                           weight="", 1)+
     labs(y="Proportion of negative PAS",
          tag="D")+
     theme(axis.text = element_text(size=12),
           axis.title = element_text(size=15),
           axis.title.x = element_text(size=13),
           plot.tag=element_text(size=20, face="bold")))


  # 4e - Correlation positive PAS vs. warming --------------------------
  (fig4e <- prop_warm_plot(traj_SProd, main_lme,
                           expression(class=="abrupt" & trend=="increase"),
                           weight="", 1)+
     labs(y="Proportion of positive PAS",
          tag="E")+
     theme(axis.text = element_text(size=12),
           axis.title = element_text(size=15),
           axis.title.x = element_text(size=13),
           plot.tag=element_text(size=20, face="bold")))

  # 4f - Correlation collapse vs. warming --------------------------

  (fig4f <- prop_warm_plot(coll, main_lme,
                           expression(did_collapse=="Yes"),
                           weight="", 1)+
      labs(y="Proportion of collapses",
           tag="F")+
     theme(axis.text = element_text(size=12),
           axis.title = element_text(size=15),
           axis.title.x = element_text(size=13),
           plot.tag=element_text(size=20, face="bold")))

  fig4 <- (prop_coll_shift + abr_coll_shift + time_coll_shift) /
    (fig4d + fig4e + fig4f)

  return(fig4)

}

fig4 <- shift_coll_plot(traj_SProd, coll_def=c("bavg", 0.25), csec_coll = 2)
# csec_coll=2 to avoid artifacts

ggsave(filename="res/figs/fig4_bavg0.25_cseccoll2.pdf",
       width=14, height=8, plot=fig4)

### Statistical tests ----------------------------------------------------

# Magnitude of abrupt declines
(t_dec <-
   t.test(coll[coll$traj=="decrease abrupt" & coll$did_collapse=="Yes",]$mag,
          coll[coll$traj=="decrease abrupt" & coll$did_collapse=="No",]$mag))

# Test for abrupt increase magnitude not possible since only one stock collapsed

