# Figure 1: two-treatment panel figure
# COI matrix input

devtools::install_github("jakelawlor/BioDivPoweR", force = T)
library(BioDivPoweR)
library(dplyr)
library(ggplot2)

# Fig S1 -------------------------------------------------------------------
# two-treatment CO1 data
pilot_two_trts <- readr::read_csv("data/species_by_sites_matrix_habitat_COI_orig_split.csv")
pilot_two_trts %>% glimpse()
pilot_two_trts$veg_score %>% table()

pilot_two_trts %>% glimpse()
pilot_two_trts <- pilot_two_trts %>%
  mutate(veg_score = recode(veg_score,
                            "high_veg" = "High Veg",
                            "low_veg" = "Low Veg"))
pilot_two_trts %>% count(veg_score)

# remove zero-occurrence species
pilot_two_trts <- pilot_two_trts[,c(TRUE, TRUE, colSums(pilot_two_trts[,-c(1:2)]) > 0 )]


# run simulation workflow -------------------------------------------------
## Create bootstraps
two_trt_boots <- bootstrap_pilot(pilot_two_trts,
                                 method = "two",
                                    "veg_score",
                                    seed = 1)
two_trt_boots %>% distinct(eff_size) %>% pull()

# get the histogram plot
p1.a <- last_plot() + 
  theme_bw() +
  theme(panel.grid = element_blank())


## Subsample boots
two_trt_sub <- subsample_boots(two_trt_boots,
                                  pilot_two_trts,
                                  method = "two",
                                  seed = 2)
two_trt_sub[[5]] %>% glimpse()
p1.b <- two_trt_sub[[3]]
p1.c <- two_trt_sub[[4]] +
  coord_cartesian(clip = "on",
                  xlim = c(-1,41),
                  ylim = c(.25,-.01),
                  expand = F) 
two_trt_sub[[5]] %>% glimpse()
rm(pilot_two_trts, two_trt_boots, two_trt_sub)

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
                             "guide-box-right")

full2 <- plot_grid(full, ggdraw(legend), rel_widths = c(9.1,1))
full2 +  ggview::canvas(12*.85,6*.85)

ggsave(full2,
       filename = "figures/raw_pdfs/figS1_two_trt_panel_CO1.pdf",
       width = 12*.85,
       height = 6*.85,
       device = cairo_pdf)

# ggsave(full2,
#        filename = "figures/raw_pdfs/figS1_two_trt_panel_CO1.png",
#        width = 12*.85,
#        height = 6*.85)

rm(list = ls())
