# Figures for Living Data Presentation Colombia


library(BioDivPoweR)
library(dplyr)
library(ggplot2)


pilot <- pilot_single_trt

extrafont::font_import()
boots <-bootstrap_pilot(pilot)
hist <- last_plot()
hist + ggview::canvas(4,3)
hist <- hist + theme(panel.grid = element_blank())
ggsave(hist, 
       filename = "figures/ldp/scenario1_hist.pdf",
       width = 4, 
       height = 3,
       device = cairo_pdf
       )

dev.off()
out <- subsample_boots(boots, pilot)
out[[3]]
out[[3]] + ggview::canvas(5.5,4)
ggsave(out[[3]], 
       filename = "figures/ldp/scenario1_percent_correct.pdf",
       width = 5.5, 
       height = 4,
       device = "pdf"
)

out[[4]] + ggview::canvas(5.5,4)
ggsave(out[[4]], 
       filename = "figures/ldp/scenario1_power.pdf",
       width = 5.5, 
       height = 4,
       device = cairo_pdf
)



log2(81/86)

log2(91/86)

86 - 81
91 - 86



# scenario 2 --------------------------------------------------------------
# here, we'll compare Fawn point and St. Ann's bank
rm(list = ls())
gc()

atl <- readr::read_csv("data/dfo_mpa/Atlantic 12S Data Filtered (1).csv")
sab <- atl %>%
  filter(stringr::str_detect(Site, "SAB"))

sab2 <- sab[,c(T,F,(colSums(sab[,-c(1,2)]) > 0))]


sab_boots <- bootstrap_pilot(sab2,
                             seed = 1,
                             n_eff_size_bins = 60)

sab_out <- subsample_boots(sab_boots, sab2,
                           target_eff_size = .176,
                           seed = 4)

sab_out[[4]] + ggview::canvas(5.5,4)
ggsave(sab_out[[4]],
       filename = "figures/ldp/scenario2_sab.pdf",
       device = cairo_pdf,
       width = 5.5,
       height = 4)

# repeat for pacific
pac <- readr::read_csv("data/dfo_mpa/Pacific 12S Data Filtered (2).csv")
pac %>% glimpse()
fawn <- pac %>%
  filter(stringr::str_detect(Site, "Fawn"))

fawn <- fawn[,c(T,F,(colSums(fawn[,-c(1,2)]) > 0))]


fawn_boots <- bootstrap_pilot(fawn, 
                              min_exp_n = 40,
                              seed = 1,
                              n_eff_size_bins = 40)

fawn_out <- subsample_boots(fawn_boots, 
                            fawn,
                            seed = 2)
fawn_out[[4]] +
  geom_hline(yintercept = .024,
             linetype = "dotdash")

ggsave(fawn_out[[4]] +
         geom_hline(yintercept = .024,
                    linetype = "dotdash"),
       filename = "figures/ldp/scenario2_fawn.pdf",
       device = cairo_pdf,
       width = 5.5,
       height = 4)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Scenario 3: comparing eDNA and Beach Sein
rm(list = ls())
gc()

# get seine data ===============
ben_seine <- readr::read_csv("data/bens_data/beach_seining_long.csv")
# reformat to matrix: site should be column 1, all species as additional columns
ben_seine <- ben_seine %>%
  select(dat_site, species, p_a) %>%
  tidyr::pivot_wider(names_from = species, 
                     values_from = p_a) 

# check to make sure no species are never seen
(colSums(ben_seine[,c(-1)]) > 0) %>% table()
# ok, one species is never seen. Remove that column.
which(colSums(ben_seine[,c(-1)]) == 0)
ben_seine[,65] <- NULL
# check again
(colSums(ben_seine[,c(-1)]) > 0) %>% table()

# get edna data ================
ben_edna <- readr::read_csv("data/bens_data/eDNA_long.csv")
ben_edna %>% glimpse()

# reformat to matrix: site should be column 1, all species as additional columns
ben_edna <- ben_edna %>%
  select(dat_site, LCT_shared, p_a) %>%
  tidyr::pivot_wider(names_from = LCT_shared, 
                     values_from = p_a) 
ben_edna %>% glimpse
ben_edna$dat_site

# check to make sure no species are never seen
(colSums(ben_edna[,c(-1)]) > 0) %>% table()


# run workflow ==============
seine_boots <- bootstrap_pilot(ben_seine, 
                               n_eff_size_bins = 35,
                               seed = 1)
p1_boots <- last_plot()
p1_boots <- p1_boots + scale_x_continuous(breaks = c(-.2, 0,.2))

edna_boots <- bootstrap_pilot(ben_edna, 
                              n_eff_size_bins = 32,
                              n_boots = 10000,
                              seed = 1
                              )
p1_edna <- last_plot()

# subsample the boots
seine_out <- subsample_boots(seine_boots, ben_seine,
                             cost_per_sample = 2230,
                             target_eff_size = .08,
                             seed = 2,
                             power = 70)
seine_out[[4]] + ggview::canvas(7,5)
seine_out[[5]]
ggsave(seine_out[[4]],
       filename = "figures/ldp/scenario3_seine.pdf",
       width = 7,
       height = 5,
       device = cairo_pdf)


edna_out <- subsample_boots(edna_boots, ben_edna,
                             cost_per_sample = 434.01,
                             target_eff_size = .08,
                            seed = 2)
edna_out[[4]] + coord_cartesian(ylim = c(0,.35)) + ggview::canvas(7,5)
edna_out[[5]] 

ggsave(edna_out[[4]] + coord_cartesian(ylim = c(0,.35)),
       filename = "figures/ldp/scenario3_edna.pdf",
       width = 7,
       height = 5,
       device = cairo_pdf)




# repeat with only 2018 data ----------------------------------------------
rm(list = ls())

ben_seine <- readr::read_csv("data/bens_data/beach_seining_long.csv")
ben_seine <- ben_seine %>%
  filter(stringr::str_detect(dat_site,"2018"))
# reformat to matrix: site should be column 1, all species as additional columns
ben_seine <- ben_seine %>%
  select(dat_site, species, p_a) %>%
  tidyr::pivot_wider(names_from = species, 
                     values_from = p_a) 

# check to make sure no species are never seen
(colSums(ben_seine[,c(-1)]) > 0) %>% table()
# ok, one species is never seen. Remove that column.
zero <- which(colSums(ben_seine[,c(-1)]) == 0) +1
ben_seine[,zero] <- NULL
# check again
(colSums(ben_seine[,c(-1)]) > 0) %>% table()


# get edna data ================
ben_edna <- readr::read_csv("data/bens_data/eDNA_long.csv")
ben_edna %>% glimpse()
ben_edna <-
  ben_edna %>%
  filter(stringr::str_detect(dat_site,"2018"))


# reformat to matrix: site should be column 1, all species as additional columns
ben_edna <- ben_edna %>%
  select(dat_site, LCT_shared, p_a) %>%
  tidyr::pivot_wider(names_from = LCT_shared, 
                     values_from = p_a) 
ben_edna %>% glimpse
ben_edna$dat_site

setequal(ben_edna$dat_site, ben_seine$dat_site)

# check to make sure no species are never seen
zero_edna <- which(colSums(ben_edna[,c(-1)]) == 0) +1
ben_edna[,zero_edna] <- NULL
# check again
(colSums(ben_edna[,c(-1)]) > 0) %>% table()

# find total richness
rich_seine <- ncol(ben_seine) - 1
rich_edna <- ncol(ben_edna) - 1

n <- 4
round(c(log2((rich_seine - n )/ rich_seine ), log2((rich_seine + n )/ rich_seine )),2)
round(c(log2((rich_edna - n )/ rich_edna ), log2((rich_edna + n )/ rich_edna )),2)

log2(57/62)

# run workflow ==============
seine_boots <- bootstrap_pilot(ben_seine, 
                               n_eff_size_bins = 35,
                               seed = 1)
p1_boots <- last_plot()
p1_boots <- p1_boots + scale_x_continuous(breaks = c(-.3, 0,.3))

edna_boots <- bootstrap_pilot(ben_edna, 
                              n_eff_size_bins = 32,
                              n_boots = 10000,
                              seed = 1
)
p1_edna <- last_plot()

# subsample the boots
seine_out <- subsample_boots(seine_boots, ben_seine,
                             cost_per_sample = 2230,
                             target_eff_size = .1,
                             seed = 4,
                             power = 75)
seine_out[[4]] + ggview::canvas(7,5)
seine_out[[5]]
ggsave(seine_out[[4]],
       filename = "figures/ldp/scenario3_seine.pdf",
       width = 7,
       height = 5,
       device = cairo_pdf)


edna_out <- subsample_boots(edna_boots, ben_edna,
                            cost_per_sample = 434.01,
                            target_eff_size = .1,
                            seed = 4,
                            power = 75)
edna_out[[4]] + coord_cartesian(ylim = c(0,.4)) + ggview::canvas(7,5)
edna_out[[5]] 

ggsave(edna_out[[4]] + coord_cartesian(ylim = c(0,.35)),
       filename = "figures/ldp/scenario3_edna.pdf",
       width = 7,
       height = 5,
       device = cairo_pdf)



