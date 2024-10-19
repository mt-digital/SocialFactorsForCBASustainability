# source("../scripts/plot.R")
# source("../scripts/plot.R")

library(gtable)
library(grid)
library(gridExtra)

figure_1 <- function() {
  # Create directories to write processed data and figures.
  for (dd in c("data", "figures")) { 
    dir.exists(dd) || dir.create(dd) 
  }


  p_jitter <- success_over_groups_jitter(
    csv_dirs = c("data/main_8-12-2024/"), 
    write_dir = "figures", size = 2.5, 
  )

  tbl = load_from_parts(c("data/main_8-12-2024/"), 
                        "data/main_asymm_sync.tmp.csv")

  this_nagents = 1000
  this_mean_degree = 20
  this_min_group_frac = 0.05

  this_group_w_innovation = 1

  print(names(tbl))
  this_tbl = tbl %>%
            filter(group_w_innovation == this_group_w_innovation &
                   nagents == this_nagents &
                   mean_degree == this_mean_degree &
                   min_group_frac == this_min_group_frac)

  p_asym_1 = asymm_heatmap(this_tbl, this_group_w_innovation, 
                      file.path("figures", 
                                paste0("asymm_", 
                                       this_group_w_innovation, 
                                       ".pdf")), 
                      cmap_limits = c(0.0, 0.8), 
                      theme_extra = theme(legend.position = "none", 
                                          plot.margin = unit(c(0,0,0,0.5), "inches"))

  )

  this_group_w_innovation = 2

  this_tbl = tbl %>%
            filter(group_w_innovation == this_group_w_innovation &
                   nagents == this_nagents &
                   mean_degree == this_mean_degree &
                   min_group_frac == this_min_group_frac)


  p_asym_2 = asymm_heatmap(this_tbl, this_group_w_innovation, 
                      file.path("figures", 
                                paste0("asymm_", 
                                       this_group_w_innovation, 
                                       ".pdf")), 
                      cmap_limits = c(0.0, 0.8), 
                      theme_extra = 
                        theme(axis.text.y = element_blank(), 
                              axis.ticks.y = element_blank(), 
                              axis.title.y = element_blank())
  )

  p_grid <- grid.arrange(p_asym_1, p_asym_2, nrow=1)

  p_grid <- grid.arrange(p_grid, p_jitter, nrow=2)

  ggsave("figures/figure_1.pdf", p_grid, width=15, height=9)
}
