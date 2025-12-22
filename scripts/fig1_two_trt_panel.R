# Figure 1: two-treatment panel figure

devtools::install_github("jakelawlor/BioDivPoweR", force = T)
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
                               seed = 1)
two_trt_sub[[5]] %>% glimpse()
p1.b <- two_trt_sub[[3]] 
p1.c <- two_trt_sub[[4]] 
two_trt_sub[[5]] %>% glimpse()

two_trt_sub[[4]] +
  geom_hline(yintercept = 0.2180057) +
  geom_vline(xintercept = 18)

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
#
#ggsave(full2,
#       filename = "figures/fig1_two_trt_panel.png",
#       width = 12*.9,
#       height = 6*.85)
#

two_trt_sub[[5]] %>% glimpse()

rm(list = ls())

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# ## raw coverage ---- 
# # find pilot coverage    
# split <- pilot_two_trts %>%
#   split(f = .$veg_score)
# dim(split[[1]])
# split[[1]] <- split[[1]][,c(T,T, colSums(split[[1]][,-c(1:2)]) > 0)]
# 
# dim(split[[2]])
# split[[2]] <- split[[2]][,c(T,T, colSums(split[[2]][,-c(1:2)]) > 0)]
# 
# split_occupancies <- purrr::map(
#   .x = split,
#   .f = ~colMeans(.x[,-c(1:2)])
# )
# 
# split_coverages <- purrr::map2(
#   .x = split_occupancies,
#   .y = split,
#   .f = ~find_coverage(.x)[1:nrow(.y)]
# )
# 
# 
# dplyr::first(which(split_coverages[[1]] >= min(max(split_coverages[[1]]), max(split_coverages[[2]]))))
# dplyr::first(which(split_coverages[[2]] >= min(max(split_coverages[[2]]), max(split_coverages[[1]]))))
# 
# n_sites_equal_coverage <- purrr::map2(
#   .x = split_coverages[[1]],
#   .y = split_coverages[[2]],
#   .f = ~c(
#     
#   )
# )
# 
# 
# # find 
# # find the first number of sites in 1:k to sample in each community
# n_sites_equal_coverage <- purrr::map2(
#   .x = eff_coverage[[1]],
#   .y = eff_coverage[[2]],
#   .f = ~c(
#     # find first number of comm1 samples where coverage is
#     # greater than the lower max of the two coverage values
#     "comm1" = dplyr::first(which(.x >= min(max(.x), max(.y)))),
#     # find first number of comm2 samples where coverage is
#     # greater than the lower max of the two coverage values
#     "comm2" = dplyr::first(which(.y >= min(max(.x), max(.y))))
#   )
# )
# 
# split_coverages_at_pilot_n <- purrr::map2(
#   .x = split_coverages,
#   .y = split,
#   .f = ~.x[nrow(.y)]
# )
# 
# 
# 
# pilot_two_trts %>%
#   split(f = .$veg_score) %>%
#   purrr::map(
#     .x = .,
#     .f = ~find_coverage(colMeans(.x[,-c(1:2)][colSums(.x[,-c(1:2) > 0])]))
#   )
# pilot_coverage <- find_coverage( colMeans(pilot[,-c(1)]))[nrow(pilot)]
# # this makes the entire coverage vector (k = 1-1000),
# # then draws the nth element, which is the number of sites in the pilot community
# 
# 

