###-###-###-###-###-###-###-###-###
#
# Analyse time series classification (supplements)
#
# 04_supplementary.R
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


# Load classified trajectories and traits -------------------------------------

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

main_lme <- readr::read_csv("data/spatial_data/main_lme.csv")

dir.create("res/figs/supp", showWarnings=FALSE)


# Table S1 - Productivity trajectories summary ----------------------------

traj_SProd_summary <-
  traj_SProd %>%
  dplyr::left_join(main_lme %>% dplyr::select(stockid, LME_NAME),
                   by="stockid") %>%
  dplyr::mutate(traj = dplyr::case_when(
    class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
    class=="abrupt" & trend=="increase" ~"increase_abrupt",
    class!="abrupt" ~traj),
    class = sub("_", " ", class),
    traj = sub("_", " ", traj),
    mag = round(mag, 2),
    loc_brk_chg = ifelse(class=="abrupt", loc_brk_chg, NA),
    mag = ifelse(class=="abrupt", mag, NA)) %>%
  dplyr::rename(`shift year` = loc_brk_chg,
                `shift magnitude` = mag,
                LME = LME_NAME) %>%
  dplyr::select(stockid, scientificname, LME,
                class, traj, `shift year`, `shift magnitude`)

readr::write_csv(traj_SProd_summary,
                 "res/figs/supp/tab_s1_trajsummary.csv")


# Table S2 - Hier.part abrupt declines --------------------------------------

df <- mod_df01_default %>%
  tidyr::drop_na(tm, Lm, DemersPelag, ordername, family, Species,
                 sst_avg_prshf, sst_change_prshf,
                 mean_ER_prshf,
                 ER_change_prshf,
                 X, Y)

env <- df %>% dplyr::select(tm, Lm,
                            DemersPelag,
                            sst_avg_prshf, sst_change_prshf,
                            mean_ER_prshf,
                            ER_change_prshf, X, Y)

# hier.part for abrupt declines:
part_dec <- hier.part::hier.part(df$shift_dec, env, fam = "binomial",
                                 link="logit", gof = "logLik")
# randomization test (LONG STEP):
set.seed(1)
rand_dec <- hier.part::rand.hp(df$shift_dec, env, fam = "binomial",
                               link="logit", gof = "logLik",
                               num.reps = 999)$Iprobs

hier.part_dec <- dplyr::bind_cols(rand_dec, part_dec$IJ) %>%
  tibble::rownames_to_column(var="variable")

readr::write_csv(hier.part_dec, "res/figs/supp/tab_s2_hier.partdec.csv")


# Table S3 - Hier.part increasing ----------------------------------------------

# hier.part for abrupt increases:
part_inc <- hier.part::hier.part(df$shift_inc, env, fam = "binomial",
                                 link="logit", gof = "logLik")

# randomization test (LONG STEP):
set.seed(1)
rand_inc <- hier.part::rand.hp(df$shift_inc, env, fam = "binomial",
                               link="logit", gof = "logLik",
                               num.reps = 999)$Iprobs

hier.part_inc <- dplyr::bind_cols(rand_inc, part_inc$IJ) %>%
  tibble::rownames_to_column(var="variable")
readr::write_csv(hier.part_inc, "res/figs/supp/tab_s3_hier.partinc.csv")



# Fig S1 - Quality scores ------------------------------------------------

# wAICc score
set.seed(12)
traj_SProd_waicc_simple <- traj_SProd %>%
  dplyr::select(stockid, class, contains("weight_aic_")) %>%
  tidyr::pivot_longer(cols = contains("weight_aic_"),
                      names_to = "class_waicc", names_prefix = "weight_aic_",
                      values_to = "waicc") %>%
  dplyr::filter(class==class_waicc) %>%
  dplyr::arrange(class, waicc) %>%
  dplyr::mutate(class_waicc = sub("_", " ", class_waicc),
                class_waicc = factor(class_waicc,
                                     levels=c("no change", "linear", "quadratic", "abrupt"))) %>%
  tibble::rowid_to_column("id") %>%

  ggplot(aes(x=class_waicc, y=waicc, fill=class_waicc)) +
  geom_boxplot(width=0.6, outliers = FALSE)+
  geom_jitter(width=0.1, height=0, size=1, shape=21)+
  scale_fill_manual(values = c("no change" = "grey80", "linear" = "grey60",
                               "quadratic" = "grey40", "abrupt" = "grey30"))+
  ylim(0, 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=15),
    axis.title.y = element_text(size=20),
    plot.tag = element_text(size=25, face="bold")
  ) +
  labs(tag="A", y="wAICc", x="")

# LOO score
traj_SProd_loo_simple <- traj_SProd %>%
  dplyr::select(stockid, class, contains("loo_")) %>%
  tidyr::pivot_longer(cols = contains("loo_"),
                      names_to = "class_loo", names_prefix = "loo_",
                      values_to = "loo") %>%
  dplyr::filter(class==class_loo) %>%
  dplyr::arrange(class, loo) %>%
  dplyr::mutate(class_loo = sub("_", " ", class_loo),
                class_loo = factor(class_loo,
                                   levels=c("no change", "linear", "quadratic", "abrupt"))) %>%
  tibble::rowid_to_column("id") %>%

  ggplot(aes(x=class_loo, y=loo, fill=class_loo)) +
  geom_boxplot(width=0.6, outliers = FALSE)+
  geom_jitter(width=0.1, height=0, size=1, shape=21)+
  scale_fill_manual(values = c("no change" = "grey80", "linear" = "grey60",
                               "quadratic" = "grey40", "abrupt" = "grey30"))+
  ylim(0, 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=15),
    axis.title.y = element_text(size=20),
    plot.tag = element_text(size=25, face="bold")
  ) +
  labs(tag="B", y="LOO", x="")

# NRMSE score
traj_SProd_nrmse_simple <- traj_SProd %>%
  dplyr::select(stockid, class, contains("nrmse_") & !nrmse_asd) %>%
  tidyr::pivot_longer(cols = contains("nrmse_"),
                      names_to = "class_nrmse", names_prefix = "nrmse_",
                      values_to = "nrmse") %>%
  dplyr::filter(class==class_nrmse) %>%
  dplyr::arrange(class, nrmse) %>%
  dplyr::mutate(class_nrmse = sub("_", " ", class_nrmse),
                class_nrmse = factor(class_nrmse,
                                     levels=c("no change", "linear", "quadratic", "abrupt"))) %>%
  tibble::rowid_to_column("id") %>%

  ggplot(aes(x=class_nrmse, y=1-nrmse, fill=class_nrmse)) +
  geom_boxplot(width=0.6, outliers = FALSE)+
  geom_jitter(width=0.1, height=0, size=1, shape=21)+
  scale_fill_manual(values = c("no change" = "grey80", "linear" = "grey60",
                               "quadratic" = "grey40", "abrupt" = "grey30"))+
  ylim(0, 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=15),
    axis.title.y = element_text(size=20),
    plot.tag = element_text(size=25, face="bold")
  ) +
  labs(tag="C", y="1-NRMSE", x="")


(traj_SProd_waicc_simple + traj_SProd_loo_simple + traj_SProd_nrmse_simple)

pdf(file = "res/figs/supp/fig_s1.pdf",
    width=15, height=5)
traj_SProd_waicc_simple + traj_SProd_loo_simple + traj_SProd_nrmse_simple
dev.off()



# Fig S2 - Trajectories in space by LMEs (abr/not abr) --------------------------

dir.create("res/figs/supp/fig_s2", showWarnings=FALSE)

sf::sf_use_s2(FALSE) # To remove error for intersection

lme <- sf::read_sf(dsn = "data/spatial_data/LME66/", layer = "lmes66gcd")
oceans <- sf::read_sf(dsn = "data/spatial_data/ocean_simpl/",
                      layer = "ocean_simpl") %>%
  dplyr::filter(grepl("Ocean", name)) %>%
  dplyr::rename(LME_NAME = name) %>%
  dplyr::mutate(LME_NUMBER = 67:73)
# Group all LMEs
lme_merge <- lme %>% dplyr::summarise()
# Substract LMEs from oceans
high_seas <- sf::st_difference(oceans, lme_merge)
# Combine LMEs and high sea regions
lmes <- dplyr::bind_rows(lme, high_seas)

lme_cntrd <- readr::read_csv("data/spatial_data/lme_ocean_centroids.csv")
main_lme <- readr::read_csv("data/spatial_data/main_lme.csv")

centroids <- lme_cntrd %>%
  dplyr::left_join(
    dplyr::left_join(traj_SProd, main_lme, by="stockid") %>%
      dplyr::group_by(LME_NAME) %>%
      dplyr::summarise(sz=n()) %>%
      dplyr::mutate(sz=sqrt(sz+5)*6),
    by="LME_NAME")

traj_SProd_sum_area <- dplyr::left_join(traj_SProd, main_lme, by="stockid") %>%
  split(.$LME_NAME) %>%
  mapply(function(x,y){
    x %>%
      dplyr::mutate(traj = dplyr::case_when(
        class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
        class=="abrupt" & trend=="increase" ~"increase_abrupt",
        class!="abrupt" ~traj)) %>%
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
      dplyr::ungroup() %>%
      sunburst_traj_plot_det(.,
                             title=NULL,
                             size_caption=8,
                             size=2, no_text=TRUE, no_caption=FALSE)

  }, x=., y=names(.), SIMPLIFY = FALSE)

centroid_lme <- tibble::tibble(x=centroids$lon,
                               y=centroids$lat,
                               area=centroids$LME_NAME,
                               width=centroids$sz) %>%
  dplyr::right_join(tibble::tibble(area = names(traj_SProd_sum_area),
                                   sunburst=traj_SProd_sum_area),
                    by = "area")

space_lme <- ggplot(lmes)+
  geom_sf()+
  theme_void()+
  ggimage::geom_subview(aes(x=x, y=y, subview=sunburst,
                            width=width, height=width), data=centroid_lme)

pdf(file = "res/figs/supp/fig_s2/fig_s2_lmes.pdf",
    width=12, height=9)
space_lme
dev.off()



# Fig S3 - Trajectory classification by taxonomic order -----------------------

dir.create("res/figs/supp/fig_s3", showWarnings=FALSE)

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
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    theme(panel.grid.major.y = element_blank(),
          axis.text=element_text(size=20),
          axis.title = element_text(size=20),
          legend.position="bottom"))


# FAO catch 1950-2022
fish_file <- readr::read_csv("data/FAO_data/Capture_2024.1.0/Capture_Quantity.csv")
species <- readr::read_csv("data/FAO_data/Capture_2024.1.0/CL_FI_SPECIES_GROUPS.csv")
area <- readr::read_csv("data/FAO_data/Capture_2024.1.0/CL_FI_WATERAREA_GROUPS.csv")
taxo <- readr::read_delim("data/FAO_data/ASFIS_sp/ASFIS_sp_2023.txt")

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
  dplyr::summarise(catch = sum(VALUE)) %>%
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
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20)))

pdf(file = paste0("res/figs/supp/fig_s3/fig_s3_taxord.pdf"),
    width=12, height=8)
print(p_ram_traj + p_fao_catch)
dev.off()


## Statistical tests -------------------------------------------------------

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



# Fig S4 - LME collapse proportions vs. warming (alternative definitions) ----

dir.create("res/figs/supp/fig_s8", showWarnings=FALSE)

# Collapse <10% maximum stock biomass
bmax0.1 <- shift_coll_plot(traj_SProd, coll_def=c("bmax", 0.1), csec_coll = 2)

# Collapse <15% average stock biomass
bavg0.15 <- shift_coll_plot(traj_SProd, coll_def=c("bavg", 0.15), csec_coll = 2)

# Collapse <50% average stock biomass
bavg0.5 <- shift_coll_plot(traj_SProd, coll_def=c("bavg", 0.5), csec_coll = 2)

fig_s4 <- bmax0.1[[2]][[3]] + bavg0.15[[2]][[3]] + bavg0.5[[2]][[3]] +
  patchwork::plot_annotation(tag_levels = "A")

ggsave("res/figs/supp/fig_s8/fig_s4_collsst.pdf", width=15, height=5, plot=fig_s4)



# Fig S5 - Trajectory overview classified on AICc only --------------------

traj_SProd_aicc <-
  readRDS("res/classif/classif_v4.61_SProd_minlen25_normTBavg_looFALSE_aicasd_start1950_AICconly.rds") %>%
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

dir.create("res/figs/supp/fig_s5", showWarnings=FALSE)

### a - Overall ---------------------------------------------------------

traj_SProd_sum_det <- traj_SProd_aicc %>%
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

pdf(file = "res/figs/supp/fig_s5/fig_s5a.pdf",
    width=8, height=6)
SProd_sunburst
dev.off()

### bc - Abrupt shifts in time -------------------------------

ts_avail_distr <- timeseries_values_views %>%
  dplyr::filter(stockid %in% (traj_SProd_aicc %>% dplyr::pull(stockid)),
                year>=1950) %>%
  dplyr::select(stockid, year, SProd) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(n=n())

ts_avail_distr2 <-
  ts_avail_distr[rep(seq_along(ts_avail_distr$n), ts_avail_distr$n), ]

abr_events <- traj_SProd_aicc %>%
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
       "res/figs/supp/fig_s5/fig_s5bc.pdf", width=5, height=5.625)


### d - Trajectories in space -------------------------------
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
  left_join(traj_SProd_aicc %>%
              dplyr::group_by(primary_FAOarea) %>%
              dplyr::summarise(sz=n()) %>%
              dplyr::mutate(sz=sqrt(sz+10)*6),
            by=c("F_CODE"="primary_FAOarea"))

traj_SProd_aicc_sum_area <- traj_SProd_aicc %>%
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


traj_SProd_aicc_sum_area <- traj_SProd_aicc %>%
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
  dplyr::right_join(tibble(area = names(traj_SProd_aicc_sum_area),
                           sunburst=traj_SProd_aicc_sum_area),
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

pdf(file = "res/figs/supp/fig_s5/fig_5d.pdf",
    width=12, height=9)
space
dev.off()


### Statistical tests ----------------------------------------------------

# Distribution in time
set.seed(1)
ks.test(ts_avail_distr2 %>%
          dplyr::pull(year),
        abr_events %>%
          dplyr::filter(trend=="decrease") %>%
          dplyr::pull(loc_brk_chg),
        simulate.p.value = TRUE, B = 1e5)
# p-value = 0.01219

set.seed(1)
ks.test(ts_avail_distr2 %>%
          dplyr::pull(year),
        abr_events %>%
          dplyr::filter(trend=="increase") %>%
          dplyr::pull(loc_brk_chg),
        simulate.p.value = TRUE, B = 1e5)
# p-value = 0.01643


## Proportions by FAO areas

# By trajectory types
class_mat_fao <- traj_SProd_aicc %>%
  dplyr::group_by(primary_FAOarea, class) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::pivot_wider(names_from=class, values_from=n) %>%
  tibble::column_to_rownames("primary_FAOarea") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(class_mat_fao, main="class by FAO area", color = TRUE)

set.seed(1)
(chi_clas_fao <- chisq.test(class_mat_fao, simulate.p.value = TRUE, B=1e5))
# p-value = 0.04485

# By trajectories (abrupt decrease, abrupt increase, nonabrupt)
traj_SProd_sum_fao <- traj_SProd_aicc %>%
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
# p-value = 0.00053

chi_traj_fao$stdres %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="FAOarea") %>%
  tidyr::pivot_longer(cols = contains("abrupt"), names_to = "traj", values_to="stdres") %>%
  dplyr::arrange(desc(abs(stdres)))


## Proportions by taxonomic order

# By trajectory types
class_mat_ord <- traj_SProd_aicc %>%
  dplyr::group_by(ordername, class) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::pivot_wider(names_from=class, values_from=n) %>%
  tibble::column_to_rownames("ordername") %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.x), 0, .))) %>%
  as.matrix()

mosaicplot(class_mat_ord, main="class by taxonomic order", color = TRUE)

set.seed(1)
(chi_clas_ord <- chisq.test(class_mat_ord, simulate.p.value = TRUE, B=1e5))
# p-value = 0.07558


# By trajectories (abrupt decrease, abrupt increase, nonabrupt)
traj_SProd_sum_ord <- traj_SProd_aicc %>%
  dplyr::mutate(traj = dplyr::case_when(
    class=="abrupt" & trend=="decrease" ~"decrease_abrupt",
    class=="abrupt" & trend=="increase" ~"increase_abrupt",
    class!="abrupt" ~"nonabrupt")) %>%
  dplyr::group_by(traj, ordername) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
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
# p-value = 1e-05

chi_traj_ord$stdres %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="taxclass") %>%
  tidyr::pivot_longer(cols = -taxclass, names_to = "traj", values_to="stdres") %>%
  dplyr::arrange(desc(abs(stdres)))


# Fig S6 - Correlation matrices ----------------------------------------------

dir.create("res/figs/supp/fig_s6", showWarnings=FALSE)

col_fn <- function(data, mapping, method="p", use="pairwise", ...){

  # grab data
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  # calculate correlation
  corr <- cor(x, y, method=method, use=use)

  # calculate colour based on correlation value
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]

  GGally::ggally_cor(data = data, mapping = mapping, col="black", ...) +
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

## Correlation matrix FishLife traits ----
fishlife_multicorr <- GGally::ggpairs(
  mod_df01_default %>%
    dplyr::select(stockid, Species, Loo, K, Winfinity, tmax, tm, M, Lm) %>%
    dplyr::mutate(across(where(is.numeric), ~ c(scale(.)))) %>%
    tidyr::drop_na(stockid) %>%
    dplyr::select(-c(stockid, Species)),
  upper = list(continuous = col_fn),
  lower = list(continuous = GGally::wrap("points", alpha = 0.3, size=0.3),
               combo = GGally::wrap("dot", alpha = 0.4, size=0.3)))

pdf(file = "res/figs/supp/fig_s6/fig_s6a_fishlife.pdf",
    width=8.75, height=8.75)
fishlife_multicorr
dev.off()

# Correlation network visualization
corrr::correlate(
  mod_df01_default %>%
    dplyr::select(stockid, Species, Loo, K, Winfinity, tmax, tm, M, Lm) %>%
    dplyr::mutate(across(where(is.numeric), ~ c(scale(.)))) %>%
    tidyr::drop_na(stockid),
  method="pearson") %>%
  corrr::network_plot(min_cor = 0.7, curved=TRUE)


## Correlation matrix all traits retained ----
traits_multicorr_all <- GGally::ggpairs(
  mod_df01_default %>%
    dplyr::mutate(across(where(is.numeric) & !contains(c("X", "Y")),
                         ~ c(scale(.)))) %>%
    tidyr::drop_na(stockid) %>%
    dplyr::select(Lm, tm, sst_change_prshf, sst_avg_prshf,
                  mean_ER_prshf, ER_change_prshf, X, Y) %>%
    dplyr::rename(SST_mean = sst_avg_prshf,
                  SST_change = sst_change_prshf,
                  ER_mean = mean_ER_prshf,
                  ER_change = ER_change_prshf),
  upper = list(continuous = col_fn),
  lower = list(continuous = GGally::wrap("points", alpha = 0.3, size=0.3),
               combo = GGally::wrap("dot", alpha = 0.4, size=0.3)))

pdf(file = "res/figs/supp/fig_s6/fig_s6b_alltraits.pdf",
    width=10, height=10)
traits_multicorr_all
dev.off()

# Correlation network visualization
corrr::correlate(
  mod_df01_default %>%
    dplyr::mutate(across(where(is.numeric) & !contains(c("X", "Y")),
                         ~ c(scale(.)))) %>%
    tidyr::drop_na(stockid) %>%
    dplyr::select(Lm, tm, sst_change_prshf, sst_avg_prshf,
                  mean_ER_prshf, ER_change_prshf, X, Y),
  method="pearson") %>%
  corrr::network_plot(min_cor = 0.7, curved=TRUE)



# Fig S7 - Abrupt shifts and alternative collapse definitions -------

dir.create("res/figs/supp/fig_s7", showWarnings=FALSE)

fig_s7 <- bmax0.1[[1]] / bavg0.15[[1]] / bavg0.5[[1]] +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(filename="res/figs/supp/fig_s7/fig_s7_coll.pdf", width=14, height=12, plot=fig_s7)

## Statistical tests ----------------------------------------------------

# Collapse <10% maximum stock biomass
coll_bmax0.1 <- add_collapsed(traj_SProd, coll_def=c("bmax", 0.1), csec_coll=2)
(t_dec_bmax0.1 <-
    t.test(coll_bmax0.1[coll_bmax0.1$traj=="decrease abrupt" &
                          coll_bmax0.1$did_collapse=="Yes",]$mag,
           coll_bmax0.1[coll_bmax0.1$traj=="decrease abrupt" &
                          coll_bmax0.1$did_collapse=="No",]$mag))
# p-value = 0.0225

# Collapse <15% average stock biomass
coll_bavg0.15 <- add_collapsed(traj_SProd, coll_def=c("bavg", 0.15), csec_coll=2)
(t_dec_bavg0.15 <-
    t.test(coll_bavg0.15[coll_bavg0.15$traj=="decrease abrupt" &
                           coll_bavg0.15$did_collapse=="Yes",]$mag,
           coll_bavg0.15[coll_bavg0.15$traj=="decrease abrupt" &
                           coll_bavg0.15$did_collapse=="No",]$mag))
# p-value = 0.062

# Collapse <50% average stock biomass
coll_bavg0.5 <- add_collapsed(traj_SProd, coll_def=c("bavg", 0.5), csec_coll=2)
(t_dec_bavg0.5 <-
    t.test(coll_bavg0.5[coll_bavg0.5$traj=="decrease abrupt" &
                          coll_bavg0.5$did_collapse=="Yes",]$mag,
           coll_bavg0.5[coll_bavg0.5$traj=="decrease abrupt" &
                          coll_bavg0.5$did_collapse=="No",]$mag))
# p-value = 0.2347

(t_inc_bavg0.5 <-
    t.test(coll_bavg0.5[coll_bavg0.5$traj=="increase abrupt" &
                          coll_bavg0.5$did_collapse=="Yes",]$mag,
           coll_bavg0.5[coll_bavg0.5$traj=="increase abrupt" &
                          coll_bavg0.5$did_collapse=="No",]$mag))
# p-value = 0.8263







# Fig S8 - Explanatory factors of PAS (alternative models) ----------

dir.create("res/figs/supp/fig_s8", showWarnings=FALSE)

## Classification with AICc only ----------------------------------------

traits_aicc <- traits_fun(traj_SProd_aicc)

readr::write_csv(traits_aicc, "res/classif/traj_SProd_classif_traits_aicc.csv")
traits_aicc <- readr::read_csv("res/classif/traj_SProd_classif_traits_aicc.csv")

mod_df02_aicc <-
  traj_SProd_aicc %>%
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
  dplyr::left_join(traits_aicc,
                   by=c("stockid", "Species")) %>%
  dplyr::mutate(DemersPelag = factor(DemersPelag, levels=c("pelagic", "demersal")))


### a - Abrupt declines -----

# Prepare dataset and scale variables
df_scl_aicc <- mod_df02_aicc %>%
  tidyr::drop_na(tm, Lm, DemersPelag, ordername, family, Species, X, Y,
                 sst_avg, sst_change, mean_ER_prshf, ER_change_prshf) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric) &
                                !dplyr::contains(c("shift", "X", "Y")),
                              ~ c(scale(.)))) %>%
  dplyr::mutate(Species = as.factor(Species),
                family = as.factor(family),
                ordername = as.factor(ordername))

# Run GAMM
modgam_dec_aicc <- mgcv::gamm(
  formula=shift_dec ~ tm + Lm +
    DemersPelag +
    sst_avg_prshf + sst_change_prshf +
    mean_ER_prshf + ER_change_prshf +
    s(X, Y) + s(length, bs="re") +
    s(ordername, bs="re") + s(family, bs="re"),
  family=binomial(link='logit'),
  na.action=na.omit,
  data=df_scl_aicc)

summary(modgam_dec_aicc$gam)
table(df_scl_aicc$shift_dec)

# Plot model estimates
(m_dec_aicc <- ggstats::ggcoef_model(
  modgam_dec_aicc$gam,
  colour = NULL, stripped_rows = FALSE,
  conf.level = 0.95,
  significance = 0.05,
  add_reference_rows=TRUE,
  point_size = 3,
  point_stroke = 0.5,
  x = "estimate"))

pdf(file = "res/figs/supp/fig_s8/fig_s8a_decaicc.pdf",
    width=6, height=4)
m_dec_aicc
dev.off()


### b - Abrupt increases -----

# Run GAMM
modgam_inc_aicc <- mgcv::gamm(
  formula=shift_inc ~ tm + Lm +
    DemersPelag+
    sst_avg_prshf + sst_change_prshf +
    mean_ER_prshf + ER_change_prshf +
    s(X, Y) + s(length, bs="re") +
    s(ordername, bs="re") + s(family, bs="re"),
  family=binomial(link='logit'),
  na.action=na.omit,
  data=df_scl_aicc)

summary(modgam_inc_aicc$gam)
table(df_scl_aicc$shift_inc)

# Plot model estimates
(m_inc_aicc <- ggstats::ggcoef_model(
  modgam_inc_aicc$gam,
  colour = NULL, stripped_rows = FALSE,
  conf.level = 0.95,
  significance = 0.05,
  add_reference_rows=TRUE,
  point_size = 3,
  point_stroke = 0.5,
  x = "estimate"))

pdf(file = "res/figs/supp/fig_s8/fig_s8b_incaicc.pdf",
    width=6, height=4)
m_inc_aicc
dev.off()


## Fishing intensity estimated relative to MSY ----

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


### c - Abrupt declines -----

# Prepare dataset and scale variables
df_scl_Umsy <- mod_df01_default %>%
  tidyr::drop_na(tm, Lm, DemersPelag, ordername, family, Species, X, Y,
                 sst_avg, sst_change, mean_U_prshf, U_change_prshf) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric) &
                                !dplyr::contains(c("shift", "X", "Y")),
                              ~ c(scale(.)))) %>%
  dplyr::mutate(Species = as.factor(Species),
                family = as.factor(family),
                ordername = as.factor(ordername))

# Run GAMM
modgam_dec_Umsy <- mgcv::gamm(
  formula=shift_dec ~ tm + Lm +
    DemersPelag +
    sst_avg_prshf + sst_change_prshf +
    mean_U_prshf + U_change_prshf +
    s(X, Y) + s(length, bs="re") +
    s(ordername, bs="re") + s(family, bs="re"),
  family=binomial(link='logit'),
  na.action=na.omit,
  data=df_scl_Umsy)

summary(modgam_dec_Umsy$gam)
table(df_scl_Umsy$shift_dec)

# Plot model estimates
(m_dec_Umsy <- ggstats::ggcoef_model(
  modgam_dec_Umsy$gam,
  colour = NULL, stripped_rows = FALSE,
  conf.level = 0.95,
  significance = 0.05,
  add_reference_rows=TRUE,
  point_size = 3,
  point_stroke = 0.5,
  x = "estimate"))

pdf(file = "res/figs/supp/fig_s8/fig_s8c_decUmsy.pdf",
    width=6, height=4)
m_dec_Umsy
dev.off()

# Stocks removed
missed_negPAS <- mod_df01_default %>%
  dplyr::filter(is.na(mean_U_prshf) & shift_dec==1) %>%
  dplyr::pull(stockid)

coll_negPAS <- add_collapsed(traj_SProd, coll_def=c("bavg", 0.25), csec_coll=2) %>%
  dplyr::filter(did_collapse=="Yes" & traj=="decrease abrupt") %>%
  dplyr::pull(stockid)

missed_negPAS[missed_negPAS %in% coll_negPAS]

### d - Abrupt increases -----

# Run GAMM
modgam_inc_Umsy <- mgcv::gamm(
  formula=shift_inc ~ tm + Lm +
    DemersPelag +
    sst_avg_prshf + sst_change_prshf +
    mean_U_prshf + U_change_prshf +
    s(X, Y) + s(length, bs="re") +
    s(ordername, bs="re") + s(family, bs="re"),
  family=binomial(link='logit'),
  na.action=na.omit,
  data=df_scl_Umsy)

summary(modgam_inc_Umsy$gam)
table(df_scl_Umsy$shift_inc)

# Plot model estimates
(m_inc_Umsy <- ggstats::ggcoef_model(
  modgam_inc_Umsy$gam,
  colour = NULL, stripped_rows = FALSE,
  conf.level = 0.95,
  significance = 0.05,
  add_reference_rows=TRUE,
  point_size = 3,
  point_stroke = 0.5,
  x = "estimate"))

pdf(file = "res/figs/supp/fig_s8/fig_s8d_incUmsy.pdf",
    width=6, height=4)
m_inc_Umsy
dev.off()
