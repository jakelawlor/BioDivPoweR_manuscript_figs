# Fig 6: Minimum Biologically-Relevant Effect figure


library(BioDivPoweR)
library(dplyr)
library(ggplot2)


# two-treatment ===================================================
# two-treatment fish data
data("pilot_two_trts")
pilot_two_trts %>% glimpse()
pilot_two_trts <- pilot_two_trts %>%
  mutate(veg_score = recode(veg_score,
                            "high_veg" = "High Veg",
                            "low_veg" = "Low Veg"))


# get values --------------------------------------------------------------
# for the figure caption, we'll get the sample sizes, richnesses, 
# coverages, and rarified richnesses in both groups

## raw richnesses ---- 
pilot_two_trts %>%
  split(f = .$veg_score) %>%
  purrr::map(
    .x = .,
    .f = ~sum(colSums(.x[,-c(1:2)]) > 0)
  )
# raw richness: 
# 81 in High veg
# 72 in Low veg




# run simulation workflow -------------------------------------------------
## Create bootstraps
two_trt_boots <- bootstrap_pilot(pilot_two_trts,
                                 method = "two",
                                 "veg_score",
                                 seed = 2)

# get the histogram plot
p1.a <- last_plot() + 
  theme_bw() +
  theme(panel.grid = element_blank())

# see the range of bin values kept
two_trt_boots %>% distinct(eff_size) %>% pull()



## Subsample boots
two_trt_sub <- subsample_boots(two_trt_boots,
                               pilot_two_trts,
                               method = "two",
                               analysis_type = "minimum",
                               effect_minimum = .1,
                               seed = 1)
two_trt_sub[[5]] %>% glimpse()
p1.b <- two_trt_sub[[3]] 
p1.c <- two_trt_sub[[4]] 
two_trt_sub[[5]] %>% glimpse()

two_trt_sub[[4]] 


# single-treatment ===================================================
# single-treatment fish data
data("pilot_single_trt")
pilot_single_trt %>% glimpse()

# find coverage of initial community
occ <- colMeans(pilot_single_trt[,-1])
find_coverage(occ, 208)

# run simulation workflow -------------------------------------------------
## Create bootstraps
single_trt_boots <- bootstrap_pilot(pilot_single_trt,
                                    seed = 2,
                                    n_eff_size_bins = 40,
                                    n_boots = 5000)


# get the histogram plot
p1.a <- last_plot() + 
  theme_bw() +
  theme(panel.grid = element_blank())


## Subsample boots
single_trt_sub <- subsample_boots(single_trt_boots, 
                                  pilot_single_trt,
                                  seed = 1,
                                  analysis_type = "minimum",
                                  effect_minimum = .1)

single_trt_sub[[5]] %>% glimpse()
p1.b <- single_trt_sub[[3]]
p1.c <- single_trt_sub[[4]]

single_trt_sub[3]
single_trt_sub[4]

#rm(pilot_single_trt, single_trt_boots, single_trt_sub)


# align plots -------------------------------------------------------------
two_trt_sub_nolegend <- two_trt_sub[[4]] +
  theme(legend.position = "none") +
  labs(y = "Minimum Detectable Richness Difference\n(Log2 fold-difference, at 80% correct)")

two_trt_sub_legend <- get_plot_component(two_trt_sub[[4]] +
                               theme(legend.key.height = unit(50,"pt"),
                                     legend.key.width = unit(9,"pt")) +
                               guides(color = guide_none()),
                             "guide")

two_trt_threshold <- plot_grid(two_trt_sub_nolegend, ggdraw(two_trt_sub_legend), 
                               rel_widths = c(5,1))

two_trt_threshold +  ggview::canvas(8*.85,6*.85)


single_trt_sub_nolegend <- single_trt_sub[[4]] +
  theme(legend.position = "none") +
  labs(y = NULL)

single_trt_sub_legend <- get_plot_component(single_trt_sub[[4]] +
                                           theme(legend.key.height = unit(50,"pt"),
                                                 legend.key.width = unit(9,"pt")) +
                                           guides(color = guide_none()),
                                         "guide")

single_trt_threshold <- plot_grid(single_trt_sub_nolegend, ggdraw(single_trt_sub_legend), 
                               rel_widths = c(5,1))

single_trt_threshold +  ggview::canvas(8*.85,6*.85)

fig6 <- plot_grid(two_trt_threshold, single_trt_threshold, rel_widths = c(1,.9),
                  labels = c("a","b")) 

fig6 + ggview::canvas(16 * .85, 6 * .85)

ggsave(fig6,
       bg = "white",
       filename = "figures/raw_pdfs/fig6_minimum_threshold_panel.pdf",
       width = 16*.85,
       height = 6*.85,
       device = cairo_pdf)


rm(list = ls())
