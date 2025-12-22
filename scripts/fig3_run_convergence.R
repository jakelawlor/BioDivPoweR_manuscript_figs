# Figure 3: convergence of recommendations at different sample sizes


library(BioDivPoweR)
library(dplyr)
library(purrr)
library(ggplot2)


# upload original data ----------------------------------------------------
data("pilot_single_trt")
pilot_orig <- pilot_single_trt
rm(pilot_single_trt)

# iterate in a loop --------------------------------------------------
run_sims <- function(start_number, n_iter, seed, n_bins){
  
  # create a vector to store data
  out <- list()
  
  
  # run trials in a loop ----------------------------------------------------
  for(i in 1:n_iter){
    
    set.seed(i * start_number)
    
    # set pilot ---------------------------------------------------------------
    # if it's the first trial, use that many rows from the original pilot
    if(i == 1){
      pilot <- pilot_orig[sample(nrow(pilot_orig),
                                 start_number,
                                 replace = F),]
      pilot <- pilot[,c(TRUE, colSums(pilot[,-c(1)]) > 0 )]
    } else {
      # if it's not the first pilot, sample number of rows from last rec
      # get last recommendation
      last_rec <-  out[[i-1]] %>% pull(recommended.sample_size.total)
      # sample that many from pilot_orig
      pilot <- pilot_orig[sample(nrow(pilot_orig),
                                 last_rec,
                                 replace = F),]
      # filter to only species that are present at least once
      pilot <- pilot[,c(TRUE, colSums(pilot[,-c(1)]) > 0 )]
      rm(last_rec)
    }
    
    # find pilot richness
  #  pilot_richness <- sum(colSums(pilot[,-1]) > 0)

    # find pilot coverage    
   # pilot_coverage <- find_coverage( colMeans(pilot[,-c(1)]))[nrow(pilot)]
    # this makes the entire coverage vector (k = 1-1000),
    # then draws the nth element, which is the number of sites in the pilot community
    

    min_exp_n <- 40
    n_bins <- 40
  # min_exp_n <- if (nrow(pilot) <= 20) {
  #   80
  # } else if (nrow(pilot) >20 & nrow(pilot) <= 30) {
  #   60
  # } else {
  #   40
  # }
  # 
  # n_bins <- if (nrow(pilot) <= 20) {
  #   80
  # } else if (nrow(pilot) >20 & nrow(pilot) <= 30) {
  #   60
  # } else {
  #   40
  # }
    

    # run workflow --------------------------------------------------------
    print(paste0("run boots. seed = ",seed, " pilot_nrow = ",nrow(pilot)))
    boots <- bootstrap_pilot(pilot,
                             #n_eff_size_bins = n_bins,
                             seed = seed,
                             min_exp_n = min_exp_n,
                             n_eff_size_bins = n_bins)
    print(paste0("run subsample. seed = ",seed, " pilot_nrow = ",nrow(pilot)))
    
    sub <- subsample_boots(boots, 
                           pilot, 
                           target_eff_size = .1,
                           seed = seed)
    
    
    # save --------------------------------------------------------------
    out[[i]] <- sub[[5]] %>%
      mutate(group = recode(group,
                            "achieved" = "start",
                            "target" = "recommended")) %>%
      rename(eff_size = min_detectable_effect) %>%
      tidyr::pivot_wider(names_from = group,
                         values_from =  c(eff_size,
                                          sample_size.total,
                                          raw_richness,
                                          coverage),
                         names_glue = "{group}.{.value}") %>% 
      relocate(contains("start"),contains("recommended"), power)
    
    rm(boots, sub, pilot)
    print(i)
    
  } # end for loop
  

  outdf <- out %>%
    bind_rows(.id = "trial") %>% 
    mutate(trial = as.numeric(trial)) %>%
    mutate(starting_n = start_number)
  
#
  #outdf <- out %>%
  #  bind_rows(.id = "trial") %>%
  #  select(-power, -min_detectable_eff_size) %>%
  #  tidyr::pivot_wider(names_from = group,
  #                     values_from = n_samples) %>%
  #  rename(start = achieved,
  #         rec = target) %>%
  #  mutate(trial = as.numeric(stringr::str_remove(trial,"run"))) %>%
  #  mutate(starting_n = start_number)
  
  return(outdf)
  
}


run_sims_1 <- run_sims(start_number = 208 * (8/8),
                       n_iter = 11,
                       seed = 1)
run_sims_2 <- run_sims(start_number = 208 * (7/8),
                       n_iter = 11,
                       seed = 2)
run_sims_3 <- run_sims(start_number = 208 * (6/8),
                       n_iter = 11,
                       seed = 3)
run_sims_4 <- run_sims(start_number = 208 * (5/8),
                       n_iter = 11,
                       seed = 4)
run_sims_5 <- run_sims(start_number = 208 * (4/8),
                       n_iter = 11,
                       seed = 5)
run_sims_6 <- run_sims(start_number = 208 * (3/8),
                       n_iter = 11,
                       seed = 6)
run_sims_7 <- run_sims(start_number = 208 * (2/8),
                       n_iter = 11,
                       seed = 7)
run_sims_8 <- run_sims(start_number =  208 * (1/8),
                       n_iter = 11,
                       seed = 7)
run_sims_9 <- run_sims(start_number =  8,
                       n_iter = 11,
                       seed = 7)
gc()


 
sims_df <- rbind(run_sims_1,
      run_sims_2,
      run_sims_3,
      run_sims_4,
      run_sims_5,
      run_sims_6,
      run_sims_7,
      run_sims_8,
      run_sims_9
      )  %>%
  mutate(starting_percent = paste0(starting_n, " (",round(starting_n/208*100,0),"%)")) %>%
  mutate(starting_percent = forcats::fct_reorder(starting_percent, rev(starting_n))) %>%
  
  mutate(trial = trial)
  

sims_sum <- sims_df %>%
  group_by(trial) %>%
  summarize(p10 = quantile(start.sample_size.total,.1),
            p90 = quantile(start.sample_size.total,.9),
            se = sd(start.sample_size.total)/sqrt(n()),
            mean = mean(start.sample_size.total))

sims_sum %>%
  filter(trial > 2) %>%
  summarize(min_p10 = min(p10),
            max_p90 = max(p90))

sims_df %>%
  filter(trial > 2) %>% 
  summarize(avghigh = quantile(start.sample_size.total, .9),
            avglow = quantile(start.sample_size.total,.1))
  


p_time <- sims_df %>% 
  ggplot(aes(x = trial,
             y = start.sample_size.total)) +
  geom_ribbon(data = sims_sum ,
              inherit.aes = F,
              aes(ymax = p90,
                  ymin = p10,
                  x = trial),
              fill = "grey50",
              alpha = .3) +
  geom_line(data = sims_sum,
              inherit.aes = F,
              aes(y = mean,
                  x = trial),
              color = "darkblue",
            linewidth = 2,
              alpha = .7) +
  geom_line(linewidth = .25,
            aes(group = starting_percent)) +
  geom_point(size = 2,
             shape = 21,
             aes(fill = starting_percent,
                 group = starting_percent)) +
  scale_fill_viridis_d() +
  theme_bw()  +
  labs(x = "Iteration",
       y = "Recommended Sample Size",
       fill = "Starting\nSample\nSize\n(proportion of\n208 sites)") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(),
                     limits = c(0,208))


p_time + ggview::canvas(5,4)



p_time2 <- sims_df %>%
  ggplot(aes(x = trial,
             y = start.sample_size.total)) +
  geom_ribbon(data = sims_sum,
              inherit.aes = F,
              aes(ymax = p90,
                  ymin = p10,
                  x = trial),
              fill = "grey50",
              alpha = .3) +
  geom_line(data = sims_sum,
            inherit.aes = F,
            aes(y = mean,
                x = trial),
            color = "grey10",
            linewidth = 2,
            alpha = .7) +
  geom_line(linewidth = .5,
            aes(group = starting_percent,
                color = starting_percent)) +
  geom_point(aes(group = starting_percent,
                 color = starting_percent),
             size = .5) +
  scale_color_viridis_d() +
  theme_bw()  +
  labs(x = "Iteration",
       y = "Sample Size\n(n sites)",
       color = "Starting\nSample\nSize\n(proportion of\n208 sites)") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(),
                     limits = c(0,208))


p_time2 + ggview::canvas(5,4)



# achieved effect size plot -----------------------------------------------
sims_sum.eff_size <- sims_df %>%
  group_by(trial) %>%
  summarize(p10 = quantile(start.eff_size,.1),
            p90 = quantile(start.eff_size,.9),
            se = sd(start.eff_size)/sqrt(n()),
            mean = mean(start.eff_size))

p_eff_size <- sims_df %>%  ggplot(aes(x = trial,
             y = start.eff_size)) +
  geom_ribbon(data = sims_sum.eff_size,
              inherit.aes = F,
              aes(ymax = p90,
                  ymin = p10,
                  x = trial),
              fill = "grey50",
              alpha = .3) +
  geom_line(data = sims_sum.eff_size,
            inherit.aes = F,
            aes(y = mean,
                x = trial),
            color = "darkblue",
            linewidth = 2,
            alpha = .7) +
  geom_line(linewidth = .25,
            aes(group = starting_percent)) +
  geom_point(size = 2,
             shape = 21,
             aes(fill = starting_percent,
                 group = starting_percent)) +
  scale_fill_viridis_d() +
  theme_bw()  +
  labs(x = "Iteration",
       y = "Min Eff Size at 80% Power",
       fill = "Starting\nSample\nSize\n(proportion of\n208 sites)") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  scale_y_reverse() +
  geom_hline(yintercept = .1,
             color = "red",
             linetype = "dashed")

p_eff_size + ggview::canvas(5,4)


p_eff_size2 <- sims_df %>%  ggplot(aes(x = trial,
             y = start.eff_size)) +
  geom_ribbon(data = sims_sum.eff_size,
              inherit.aes = F,
              aes(ymax = p90,
                  ymin = p10,
                  x = trial),
              fill = "grey50",
              alpha = .3) +
  geom_line(data = sims_sum.eff_size,
            inherit.aes = F,
            aes(y = mean,
                x = trial),
            color = "grey10",
            linewidth = 2,
            alpha = .7) +
  geom_line(linewidth = .5,
            aes(group = starting_percent,
                color = starting_percent)) +
  geom_point(aes(group = starting_percent,
                 color = starting_percent),
             size = .5) +
  scale_color_viridis_d() +
  theme_bw()  +
  labs(x = "Iteration",
       y = "Min Eff Size at 80% Power",
       color = "Starting\nSample\nSize\n(proportion of\n208 sites)") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  scale_y_reverse()

p_eff_size2 + ggview::canvas(5,4)


# make plot of richness ---------------------------------------------------
sims_sum.richness <- sims_df %>%
  group_by(trial) %>%
  summarize(p10 = quantile(start.raw_richness,.1),
            p90 = quantile(start.raw_richness,.9),
            se = sd(start.raw_richness)/sqrt(n()),
            mean = mean(start.raw_richness))

p_rich <- sims_df %>% ggplot(aes(x = trial,
             y = start.raw_richness)) +
  geom_ribbon(data = sims_sum.richness,
              inherit.aes = F,
              aes(ymax = p90,
                  ymin = p10,
                  x = trial),
              fill = "grey50",
              alpha = .3) +
  geom_line(data = sims_sum.richness,
            inherit.aes = F,
            aes(y = mean,
                x = trial),
            color = "darkblue",
            linewidth = 2,
            alpha = .7) +
  geom_line(linewidth = .25,
            aes(group = starting_n)) +
  geom_point(size = 2,
             shape = 21,
             aes(fill = starting_percent,
                 group = starting_percent)) +
  scale_fill_viridis_d()  +
  theme_bw()  +
  labs(x = "Iteration",
       y = "Detected Richness",
       fill = "Starting\nSample\nSize\n(proportion of\n208 sites)") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) 


p_rich2 <- sims_df %>%  ggplot(aes(x = trial,
             y = start.raw_richness)) +
  geom_ribbon(data = sims_sum.richness,
              inherit.aes = F,
              aes(ymax = p90,
                  ymin = p10,
                  x = trial),
              fill = "grey50",
              alpha = .3) +
  geom_line(data = sims_sum.richness,
            inherit.aes = F,
            aes(y = mean,
                x = trial),
            color = "grey10",
            linewidth = 2,
            alpha = .7) +
  geom_line(linewidth = .5,
            aes(group = starting_percent,
                color = starting_percent)) +
  geom_point(aes(group = starting_percent,
                 color = starting_percent),
             size = .5) +
  scale_color_viridis_d()  +
  theme_bw()  +
  labs(x = "Iteration",
       y = "Detected Richness",
       color = "Starting\nSample\nSize\n(proportion of\n208 sites)") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) 


p_rich2 + ggview::canvas(5,4)


# convergence plot --------------------------------------------------------
sims_df %>%
  ggplot(aes(y = recommended.sample_size.total, 
           x = start.sample_size.total,
           group = starting_n)) +
  geom_path(aes(color = starting_n)) +
  geom_abline(slope = 1,
              intercept = 0) +
  coord_equal() +
  scale_color_viridis_c() +
  theme_minimal() +
  facet_wrap(~starting_n)


# add in 1:1 points for every recommendation
outdf2 <- sims_df %>%
  # make the extra rows
  mutate(start.sample_size.total = recommended.sample_size.total, 
         .after = recommended.sample_size.total,
         opt = "one_one") %>%
  # duplicate dataset, one copy is original, one is extra rows
  bind_rows(sims_df %>% mutate(opt = "actual"), .) %>%
  arrange(trial, opt) %>%
  # filter out last 1-to-1 point
  filter(!(trial == 11),
         !(trial == 10 &  opt == "one_one")) 

# add starting numbers
start_numbers_df <- 
  data.frame(trial = c(-1),
             starting_n = unique(outdf2$starting_n),
             start.sample_size.total = outdf2 %>% 
               filter(trial == 1) %>%
               pull(starting_n),
             recommended.sample_size.total= 0
               )

outdf3 <- outdf2 %>%
  dplyr::bind_rows(start_numbers_df) %>%
  arrange(desc(starting_n), trial) %>%
  mutate(starting_percent = paste0(round(starting_n/208*100,0),"% (",starting_n,")")) %>%
  mutate(starting_percent = forcats::fct_reorder(starting_percent, rev(starting_n))) 
  

p_spiral <- outdf3 %>%
  ggplot(aes(y = recommended.sample_size.total, 
             x = start.sample_size.total,
             group = starting_n)) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "longdash") +
  geom_path(linewidth = .25) +
  geom_point(data = . %>%
               filter(start.sample_size.total != recommended.sample_size.total),
               shape = 21,
             size = 2.5,
             alpha = .75,
             aes(fill = trial)) +
  coord_equal(ylim = c(15,215),
              xlim = c(15,215)) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Pilot Sample Size",
       y = "Recommended Sample Size",
       fill = "Iteration") +
  theme(legend.key.height = unit(40,"pt"),
        legend.key.width = unit(10,"pt"))

p_spiral + ggview::canvas(5,4)


p_spiral2 <- outdf3 %>%
  ggplot(aes(y = recommended.sample_size.total, 
             x = start.sample_size.total,
             group = starting_n,
             color= factor(starting_n))) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "longdash") +
  geom_path(linewidth = .25) +
  coord_equal(ylim = c(15,215),
              xlim = c(15,215)) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Pilot Sample Size",
       y = "Recommended Sample Size",
       fill = "Iteration") 

p_spiral2 + ggview::canvas(5,4)

p_spiral3 <- outdf3 %>%
  ggplot(aes(y = recommended.sample_size.total, 
             x = start.sample_size.total,
             group = starting_n,
             color= as.factor(starting_n))) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "longdash") +
  geom_path(linewidth = .25,
            alpha = .6,
            show.legend = F,
            color = "black") +
  geom_point(data = . %>% filter(trial == 15),
             shape = 21, 
             stroke = .3,
             aes(fill = as.factor(starting_n)),
             color = "black") +
  coord_equal(ylim = c(15,215),
              xlim = c(15,215)) +
  theme_bw() +
  labs(x = "Pilot Sample Size",
       y = "Recommended Sample Size",
       fill = "Starting Sample Size") 

p_spiral3 + ggview::canvas(5,4)

p_spiral4 <- outdf3 %>%
  ggplot(aes(y = recommended.sample_size.total, 
             x = start.sample_size.total,
             group = starting_n,
             color= trial)) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "longdash") +
  geom_path(linewidth = .5,
            alpha = .9) +
  geom_point(data = . %>% filter(trial == max(trial)), 
             stroke = .5,
             size = 2,
             shape = 21,
             fill = "green",
             color = "black") +
  geom_point(data = . %>% filter(trial == min(trial) &
                                  opt == "actual" ), 
             stroke = .5,
             size = 2,
             shape = 21,
             fill = "grey80",
             color = "black") +
  coord_equal(ylim = c(15,215),
              xlim = c(15,215)) +
  theme_bw() +
  scale_color_gradient(high = "black",
                       low = "grey80",
                       breaks = c(1,3,5,7,9,11),
                       labels = ~. + 1) +
  labs(x = "Pilot Sample Size",
       y = "Recommended Sample Size",
       fill = "Starting Sample Size",
       color = "Iteration") 

p_spiral4 + ggview::canvas(5,4)

p_spiral5 <- outdf3 %>%
  ggplot(aes(y = recommended.sample_size.total, 
             x = start.sample_size.total,
             group = starting_n,
             color= trial)) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "longdash") +
  geom_path(aes(color = starting_percent),
            linewidth = .5,
            alpha = .9) +
  geom_point(data = . %>% filter(trial == min(trial) &
                                   opt == "actual" ), 
             stroke = .5,
             size = 2,
             shape = 4,
             aes(color = starting_percent)) +
  geom_point(data = . %>% filter(trial == max(trial)), 
             stroke = .5,
             size = 2,
             shape = 23,
             aes(fill = starting_percent),
             color = "black") +
  coord_equal(ylim = c(15,215),
              xlim = c(15,215)) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x = "Pilot Sample Size",
       y = "Recommended Sample Size",
       fill = "Starting Sample Size",
       color = "Starting Sample Size") 

p_spiral5 + ggview::canvas(5,4)



p_spiral6 <- outdf3 %>%
  filter(starting_percent %in% c("4% (8)","50% (104)","100% (208)")) %>%
  ggplot(aes(y = recommended.sample_size.total, 
             x = start.sample_size.total,
             group = starting_n,
             color= trial)) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "longdash") +
  geom_path(aes(color = starting_percent),
            linewidth = .5,
            alpha = .9) +
  geom_point(data = . %>% filter(trial == min(trial) &
                                   opt == "actual" ), 
             stroke = .5,
             size = 2,
             shape = 4,
             aes(color = starting_percent)) +
  geom_point(data = . %>% filter(trial == max(trial)), 
             stroke = .5,
             size = 2,
             shape = 23,
             aes(fill = starting_percent),
             color = "black") +
  coord_equal(ylim = c(15,215),
              xlim = c(15,215)) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x = "Pilot Sample Size",
       y = "Recommended Sample Size",
       fill = "Starting Sample Size",
       color = "Starting Sample Size") 

p_spiral6 + ggview::canvas(5,4)




# align plots -------------------------------------------------------------
library(cowplot)
# aligned <- cowplot::align_plots(
#   plotlist = list(
#     p_rich + theme(legend.position = "none",
#                    axis.text.x = element_blank(),
#                    axis.title.x = element_blank(),
#                    axis.ticks.x = element_blank(),
#                    plot.background = element_blank()) +
#       labs(y = "Detected\nRichness"),
#     p_eff_size + theme(legend.position = "none",
#                        axis.text.x = element_blank(),
#                        axis.title.x = element_blank(),
#                        axis.ticks.x = element_blank(),
#                        plot.background = element_blank()) +
#       labs(y = "Detectable\nChange"), 
#     p_time + theme(legend.position = "bottom",
#                    legend.title.position = "top",
#                    plot.background = element_blank(),
#                    legend.key.spacing = unit(2,"pt")) +
#       labs(fill = "Sample size (proportion of 208 sites)") +
#       guides(fill = guide_legend(nrow = 1)), 
#     p_spiral6 + theme(legend.position = "bottom",
#                       legend.direction = "horizontal",
#                       legend.title.position = "top",
#                       legend.key.height = unit(8,"pt"),
#                       legend.key.width = unit(40,"pt"),
#                       plot.background = element_blank(),
#                       legend.margin = margin(t = 0)) +
#       coord_cartesian(xlim = c(0,210),
#                       ylim = c(0,210),
#                       expand = F)))
# 
# p1_v1 <- plot_grid(
#   plot_grid(aligned[[1]],
#           aligned[[2]],
#           aligned[[3]],
#           ncol = 1,
#           rel_heights = c(1,1,3),
#           align = "v"),
#   aligned[[4]],
#   axis = "b") +
#   theme(plot.background = element_rect(fill = "white",
#                                        color = "transparent")) 
# 
# p1_v1  + ggview::canvas(10,5.75)
# 
# ggsave(p1_v1, 
#        filename = "figures/fig3_time_and_spiral_panel.pdf",
#        device = cairo_pdf(),
#        width = 10,
#        height = 5.75)
# 

# try again ---------------------------------------------------------------

aligned2 <- cowplot::align_plots(
  plotlist = list(
    p_rich2 + theme(legend.position = "none",
                   axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   plot.background = element_blank()) +
      labs(y = "Detected\nRichness"),
    p_eff_size2 + theme(legend.position = "none",
                       axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       plot.background = element_blank()) +
      labs(y = "Detectable\nChange"), 
    p_time2 + theme(legend.position = "none",
                   legend.title.position = "top",
                   plot.background = element_blank(),
                   legend.key.spacing = unit(2,"pt")) +
      labs(color = "Sample size (proportion of 208 sites)") +
      guides(color = guide_legend(nrow = 1)), 
    p_spiral6 + theme(legend.position = "none",
                      legend.direction = "horizontal",
                      legend.title.position = "top",
                      legend.key.height = unit(8,"pt"),
                      legend.key.width = unit(40,"pt"),
                      plot.background = element_blank(),
                      legend.margin = margin(t = 0)) +
      coord_cartesian(xlim = c(0,210),
                      ylim = c(0,210),
                      expand = F)))

shared_legend <- get_plot_component(p_time2 +
                                      theme(legend.position = "bottom",
                                            legend.title.position = "top",
                                            plot.background = element_blank(),
                                            legend.background = element_blank(),
                                            legend.key.spacing = unit(2,"pt")) +
                                      labs(color = "Starting sample size\n(proportion of 208 sites)") +
                                      guides(color = guide_legend(nrow = 1)) +
                                      theme(legend.box.margin = margin(0,0,0,0)),
                                    "guide-box-bottom")

ggdraw(shared_legend)

p1_v2 <- plot_grid(
  plot_grid(aligned2[[1]],
            aligned2[[2]],
            aligned2[[3]],
            ncol = 1,
            rel_heights = c(1,1,2.5),
            align = "v",
            labels = c("a","b","c"),
            label_x = -.035,
            label_y = c(1.03,1.03,1)),
  aligned2[[4]],
  axis = "b",
  labels = c(NA,"d"))  +
  theme(plot.background = element_rect(fill = "white",
                                       color = "transparent"),
        plot.margin = margin(b = 40,
                             l = 12)) +
draw_grob(shared_legend,
            x = .5, 
            y = -.55,
            hjust = .5, 
            vjust = 0)

    
p1_v2  + ggview::canvas(8,4.5)

ggsave(p1_v2, 
       filename = "figures/raw_pdfs/fig3_time_and_spiral_panel_v2.pdf",
       device = cairo_pdf(),
       width = 8,
       height = 4.5)



# just the two plots ------------------------------------------------------

p1_v3 <- plot_grid(
  p_time2 + theme(legend.position = "none",
                  legend.direction = "horizontal",
                  legend.title.position = "top",
                  plot.background = element_blank(),
                  legend.key.spacing = unit(2,"pt")) +
    labs(color = "Starting Sample Size") +
    guides(color = guide_legend(nrow = 2)),
  p_spiral6 + theme(legend.position = "none",
                    legend.direction = "horizontal",
                    legend.title.position = "top"),
  align = "h",
  labels = c("a","b")
) +
  theme(plot.margin = margin(b = 40)) +
  draw_grob(shared_legend,
            x = .5, 
            y = -.55,
            hjust = .5, 
            vjust = 0) 

p1_v3 + ggview::canvas(8,4)

ggsave(p1_v3, 
       filename = "figures/raw_pdfs/fig3_time_and_spiral_panel_v3.pdf",
       device = cairo_pdf(),
       width = 8,
       height = 4)

ggsave(p1_v3, 
       filename = "figures/fig3_time_and_spiral_panel_v3.png",
       width = 8,
       height = 4)




# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------

test <- sub[[4]]


data2 <- test %>%
  ggplot2::layer_data() %>% 
  filter(x > 17)

test +
  geom_vline(xintercept = 17)

data2 %>%
  ggplot(aes(x = x, y=y)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = -.1)
  glimpse()

  
  
  # cut out everything below absolute sample number of 6 (arbitrary?) 
  # keep 6 - cut out below 6. 
  # fix downsampling step so that x values aren't repeated 
  # keep the highest coverage for each sample number 
  