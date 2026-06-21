# Figure 2: single-treatment panel figure

library(BioDivPoweR)
library(dplyr)
library(ggplot2)



# Fig 2 -------------------------------------------------------------------
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
                                  )

single_trt_sub[[5]] %>% glimpse()
p1.b <- single_trt_sub[[3]]
p1.c <- single_trt_sub[[4]]

#rm(pilot_single_trt, single_trt_boots, single_trt_sub)

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

full2 <- plot_grid(full, ggdraw(legend), rel_widths = c(9.1,1))
full2 +  ggview::canvas(12*.85,6*.85)

ggsave(full2,
       filename = "figures/raw_pdfs/fig2_single_trt_panel.pdf",
       width = 12*.85,
       height = 6*.85,
       device = cairo_pdf)



single_trt_sub[[5]] %>% glimpse()
#Rows: 1
#Columns: 6
#$ group                 <chr> "achieved"
#$ power                 <chr> "80"
#$ min_detectable_effect <dbl> 0.08540994
#$ sample_size.total     <int> 208
#$ coverage              <dbl> 0.9250003
#$ raw_richness          <int> 86


rm(list = ls())
