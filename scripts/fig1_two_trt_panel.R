# Figure 1: two-treatment panel figure

# install package BioDivPoweR from attached
# (remotes::install_github() removed for anonymity)
library(BioDivPoweR)
library(dplyr)
library(ggplot2)


# Fig 1 -------------------------------------------------------------------
# two-treatment fish data
data("pilot_two_trts")
pilot_two_trts %>% glimpse()
pilot_two_trts <- pilot_two_trts %>%
  mutate(veg_score = recode(veg_score,
                            "high_veg" = "High Veg",
                            "low_veg" = "Low Veg"))


# get values --------------------------------------------------------------
# for the figure caption, we'll get the sample size, richness, 
# coverage, and rarefied richness in both groups

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
                               seed = 1)
two_trt_sub[[5]] %>% glimpse()
p1.b <- two_trt_sub[[3]] 
p1.c <- two_trt_sub[[4]] 
two_trt_sub[[5]] %>% glimpse()

two_trt_sub[[4]]  + ggview::canvas(9,5.5)

ggsave(two_trt_sub[[4]], 
       filename = "example_output.png",
       width = 9, 
       height = 5.5,
       unit = "in")

#rm(pilot_two_trts, two_trt_boots, two_trt_sub)

# align plots -------------------------------------------------------------
library(cowplot)
# remove legends before combining
p1.a_nolegend <- p1.a + 
  theme(legend.position = "none")  +
  labs(y = "# Simulations")
  
p1.b_nolegend <- p1.b + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -70,
                                   hjust = 0,
                                   vjust = .5)) +
  labs(y = "Proportion Correct")

p1.c_nolegend <- p1.c + 
  theme(legend.position = "none") +
  labs(y = "Minimum Detectable Richness Difference\n(Log2 fold-difference, at 80% correct)")

rm(p1.a, p1.b)

fig1_alist <- align_plots(p1.a_nolegend, 
                          #p1.b_nolegend, 
                          p1.c_nolegend, 
                          align = "hv", 
                          axis = "tblr",
                          greedy = F)


# Stack left column without extra space
col1 <- plot_grid(fig1_alist[[1]],
                  p1.b_nolegend + theme(axis.title.y = element_text(margin = margin(r=7))),
                  ncol = 1,
                  rel_heights = c(1.1, 1),  # equal square plots
                  align = "v",            # vertical alignment only inside column
                  axis = "lr",
                  labels = c("a","b"))            # align left/right axes only

# Combine with right column
full <- plot_grid(col1, fig1_alist[[2]],
                  labels = c(NA,"c"),
                  ncol = 2,
                  rel_widths = c(1, 2),
                  axis = "tb")           # align top/bottom axes

full + ggview::canvas(11*.85,6*.85)

legend <- get_plot_component(p1.c +
                               theme(legend.key.height = unit(50,"pt"),
                                     legend.key.width = unit(9,"pt")) +
                               guides(color = guide_none()),
                             "guide")

full2 <- plot_grid(full, ggdraw(legend), rel_widths = c(9.5,1))
full2 +  ggview::canvas(12*.9,6*.85)


ggsave(full2,
       filename = "figures/raw_pdfs/fig1_two_trt_panel.pdf",
       width = 12*.9,
       height = 6*.85,
       device = cairo_pdf)


two_trt_sub[[5]] %>% glimpse()

rm(list = ls())

