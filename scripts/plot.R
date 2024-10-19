require(ggplot2)
require(latex2exp)
require(dplyr)
require(readr)
require(purrr)
require(reshape2)
require(stringr)
require(tidyr)


mytheme = theme(axis.line = element_line(), legend.key=element_rect(fill = NA),
                text = element_text(size=22),# family = 'PT Sans'),
                legend.key.width = unit(2, 'cm'),
                legend.key.size = unit(1.5, 'lines'),
                panel.background = element_rect(fill = "white"))


group_start_remap <- function(group_start) {
  if (group_start == "1")
    return ("Minority")
  else if (group_start == "2")
    return ("Majority")
  else
    return ("Both")
}

success_over_groups_compare_maxsuccess_minfrac <- function(supp_dir = "data/supp_network_maxsteps/",
                                                           write_path = "figures/compare_max_m.pdf",
                                                           tbl_override = NULL) {
  if (is.null(tbl_override)) {
    full_df <- load_from_parts(csv_dirs, "data/supp_network_maxsteps_sync.csv")
  } else {
    full_df <- tbl_override
  }
  
  aggdf <- full_df %>%
    filter(mean_degree == 6) %>% filter(nagents == 1000) %>%
    group_by(nagents, min_homophily, maj_homophily, group_w_innovation, min_group_frac) %>%
    summarize(success_rate = mean(frac_a_curr_trait)) %>%
    group_by(min_group_frac, group_w_innovation) %>%
    summarize(max_success_rate = max(success_rate))
  
  aggdf <- rename(aggdf, start_group = group_w_innovation) 
  aggdf$start_group[aggdf$start_group == 1] <- "Minority"
  aggdf$start_group[aggdf$start_group == 2] <- "Majority"
  aggdf$start_group <- factor(aggdf$start_group, level=c("Majority", "Minority", "Both"))
  aggdf$min_group_frac <- factor(aggdf$min_group_frac)
  
  p <- ggplot(aggdf, aes(x = start_group,  #factor(start_group, level=c("Majority", "Minority", "Both")), 
                         y = max_success_rate,
                         group = min_group_frac)) + 
    geom_point(aes(color=min_group_frac), size=3.5) +
    geom_line(aes(color=min_group_frac), size=2) +
    scale_color_discrete(name=TeX("Minority frac., $m$")) +
    xlab("Start group") + 
    ylab("Maximum success rate") + 
    mytheme
  
  ggsave(write_path, width=8, height=3.75)
  
  return (p)
}


success_over_groups_compare_kbar <- function(main_dir = "data/main_network_maxsteps",
                                             supp_dir = "data/supp_network_maxsteps",
                                             write_filename = "figures/compare_max_k",
                                             tbl_override = NULL) {
  
  # Filter out files that are not k sensitivity files from supp_dir.
}


success_over_groups_compare_maxsuccess_N <- function(supp_dir = "data/supp_network_maxsteps/",
                                         write_path = "figures/compare_max_N.pdf",
                                         write_filename = "compareN.pdf",
                                         tbl_override = NULL) {
  if (is.null(tbl_override)) {
    full_df <- load_from_parts(csv_dirs, "data/supp_network_maxsteps_sync.csv")
  } else {
    full_df <- tbl_override
  }
  
  aggdf <- full_df %>%
    filter(mean_degree == 6) %>%
    group_by(nagents, min_homophily, maj_homophily, group_w_innovation, min_group_frac) %>%
    summarize(success_rate = mean(frac_a_curr_trait)) %>%
    group_by(nagents, group_w_innovation) %>%
    summarize(max_success_rate = max(success_rate))
  aggdf <- rename(aggdf, start_group = group_w_innovation) 
                  # minority_group_size = min_group_frac,
                  # h_min = min_homophily, h_maj = maj_homophily)

  aggdf$start_group <- factor(aggdf$start_group)#, level=c("Majority", "Minority", "Both"))
  aggdf$nagents <- factor(aggdf$nagents)
  
  p <- ggplot(aggdf, aes(x = start_group,  #factor(start_group, level=c("Majority", "Minority", "Both")), 
                         y = max_success_rate)) + 
    geom_point(aes(color=nagents)) +
    xlab("Start group") + 
    ylab("Success rate") + 
    mytheme
  
  # ggsave(write_path, width=6, height=3.75)
  
  return (p)
}


success_over_groups_jitter <- function(
    csv_dirs = c("data/main_parts"), 
    write_dir = "figures",
    write_filename = "success_over_groups_jitter.pdf",
    tbl_override = NULL,
    this_ylim = c(NA, NA),
    size=1.2,   # Begin geom_jitter arguments.
    alpha=0.1, 
    width=0.2, 
    color="#1b74bc", 
    stroke=NA, 
    ...  # provide additional geom_jitter arguments se quiser.
  ) {

  if (is.null(tbl_override)) {
    full_df <- load_from_parts(csv_dirs, "data/success_over_groups.csv")
  } else {
    full_df <- tbl_override
  }

  # full_df <- load_from_parts(csv_dirs, "data/success_over_groups_jitter.csv")
  write_path <- file.path(write_dir, write_filename)
  
  # Need each outcome for every homophily value separated.
  aggdf <- full_df %>%
    group_by(min_homophily, maj_homophily, group_w_innovation, min_group_frac) %>%
    summarize(success_rate = mean(frac_a_curr_trait))
  
  aggdf$group_w_innovation <- map_chr(aggdf$group_w_innovation, group_start_remap)
  aggdf$min_group_frac <- factor(aggdf$min_group_frac)
  
  aggdf <- rename(aggdf, start_group = group_w_innovation, 
                  minority_group_size = min_group_frac,
                  h_min = min_homophily, h_maj = maj_homophily)
  
  # Now create boxplot with start group on x-axis, success_rate on y-axis, 
  # and color-coded by minority group size.
  p_ret <- ggplot(aggdf, aes(x = factor(start_group, 
                                    level=c("Majority", "Minority", "Both")), 
                         y = success_rate)) + 
    # geom_violin() + 
    geom_jitter(size=size, alpha=alpha, width=width, 
                color=color, stroke=stroke, ...) +
    stat_summary(fun=mean, geom="line", color="black", size=1.5, aes(group=1)) +
    stat_summary(fun=mean, geom="point", shape="circle", #shape=23, 
                 size=1.5, alpha=0.8, fill="black", stroke=1.1) +
    xlab("Start group") + ylab("Success rate") + 
      # ylim(this_ylim) +
    mytheme
  

  ggsave(write_path, width=6, height=3.75)
  
  return (p_ret)
}

line_for_jitter_legend <- function(write_path = "figures/jitter_legend_line.pdf") {
  df <- data.frame(x=c(0, 1, 2), y=c(1, 1, 1))
  p <- ggplot(df, aes(x=x,y=y)) + 
    geom_line(color="#02bec3", size=1.5) + 
    geom_point(fill="#02bec3", shape=23, size=4) + 
    mytheme 
  
  ggsave(write_path, p)
  
  return (p)
}

point_for_jitter_legend <- function(write_path = "figures/jitter_legend_point.pdf") {
  df <- data.frame(x=c(0, 1, 2), y=c(1, 1, 1))
  p <- ggplot(df, aes(x=x,y=y)) + 
    geom_point(alpha=0.1, size=20) + 
    mytheme 
  
  ggsave(write_path, p)
  
  return (p)
}


steps_over_groups_success_failure <- function(csv_dirs = c("data/main_parts"), 
                                              write_dir = "figures/",
                                              write_filename = "steps_over_groups_success_failure.pdf", 
                                              minority_pop_sizes = c(0.05),
                                              tbl_override = NULL,
                                              this_ylim = c(0, 80)) {
  
  if (is.null(tbl_override)) {
    full_df <- load_from_parts(csv_dirs, "data/steps_over_groups.csv")
  } else {
    full_df <- tbl_override
  }
  
  write_path <- file.path(write_dir, write_filename)
  
  success_remap <- function(frac_a_curr_trait) {
    if (frac_a_curr_trait == 1.0)
      return ("Success")
    else
      return ("Failure")
  }
    full_df$success = factor(map_chr(full_df$frac_a_curr_trait, success_remap))
# full_df
    aggdf <- full_df %>%
      filter(min_group_frac %in% minority_pop_sizes) %>%
      group_by(group_w_innovation, success) %>%
      summarize(step = mean(step))
    
    aggdf$group_w_innovation <- map_chr(aggdf$group_w_innovation, group_start_remap)

      aggdf <- rename(aggdf, start_group = group_w_innovation)
    aggdf$Status <- factor(aggdf$success)
    
    color_values <- c("#1896ee", "#f0a0ad") 
    # color_values <- c("#56B4E9")
    p <- ggplot(aggdf, aes(x = factor(start_group,
                                level=c("Majority", "Minority", "Both")),
                      y=step, group=Status)) +
      geom_line(aes(color=Status), size=1.5) +
      geom_point(aes(fill=Status, shape=Status), color='black', size=5, stroke=1.1) +
      xlab("Start group") + ylab("Mean steps to\nsuccess/failure") + 
      ylim(this_ylim) +
      # scale_color_discrete("Status") +
      scale_shape_manual(values = c(21,24), breaks=c("Success", "Failure")) +
      scale_color_discrete(limits=c("Success", "Failure"), type=color_values) +
      scale_fill_discrete(limits=c("Success", "Failure"), type = color_values) +
      mytheme
    
    ggsave(write_path, width=7, height=3.25)
    
    return(p)
}


main_asymm_heatmaps <- function(csv_dir = "data/main_parts", 
                                write_dir = "figures/heatmaps/main", 
                                measure = "sustainability",
                                cmap_limits = c(0.0, 0.8))
{
  
  tbl <- load_from_parts(c(csv_dir), "data/main_asymm_sync.tmp.csv")
  
  group_w_vals <- c("1", "2", "Both")
  for (group_w_innovation in group_w_vals) {
    
    asymm_heatmap(tbl, group_w_innovation, 
                  file.path(write_dir, paste0(group_w_innovation, ".pdf")),
                  measure, cmap_limits)
  }
}


load_from_parts <- function(csv_dirs = c("data/main_parts", "data/supp_parts"),
                            sync_file = NULL, col_select = NULL) {
  
  group_w_vals <- c("1", "2", "Both")
  
  if (!is_null(sync_file)) {
    if (file.exists(sync_file)) {
      return (read_csv(sync_file))
    }
  }
  
  for (group_w_innovation in group_w_vals) {
    
    dir_idx = 1
    
    for (dd in csv_dirs) {
      files_part <- list.files(dd, 
                               pattern = paste0("group_w_innovation=", group_w_innovation),
                               full.names = TRUE)
      if (dir_idx == 1) {
        files <- files_part 
      }
      else {
        files <- rbind(files, files_part)
      }
    }
    
    tbl_part <- files %>% 
      map_df(~read_csv(., show_col_types = FALSE))
      # map_df(~read_csv(., col_select=col_select, show_col_types = FALSE))
      # ^Not working currently to have col_select = NA; will have to fix for supp I think.
    
    tbl_part$group_w_innovation = group_w_innovation
    
    if (group_w_innovation == "1") {
      tbl <- tbl_part
    }
    else {
      tbl <- rbind(tbl, tbl_part)
    }
  }
  
  # If we got here and the sync_file is not null, then there was no sync_file
  # to read and we create it now to sync the full output dataframe.
  if (!is_null(sync_file)) {
    write_csv(tbl, sync_file)
  }
  
  return (tbl)
}



supp_asymm_heatmaps <- function(csv_dir = "data/supp_parts", 
                                tbl = NULL,
                                write_dir = "figures/supp",
                                do_nagents = TRUE,
                                do_mean_degree = TRUE,
                                do_min_group_frac = TRUE,
                                col_select = c("nagents", "mean_degree", "min_group_frac",
                                               "min_homophily", "maj_homophily", 
                                               "frac_a_curr_trait", "step")) {
  
  
  if (is_null(tbl)) {
    tbl <- load_from_parts(csv_dir, "data/supp_asymm_heatmaps.tmp.csv", col_select) 
  }

  group_w_vals <- c("1", "2", "Both")
  
  for (this_group_w_innovation in group_w_vals) {

    # nagents sensitivity.
    if (do_nagents) {
      for (this_nagents in c(50, 100, 2000)) {
        
        if (this_nagents == 50) {
          this_mean_degree = 5
          this_min_group_frac = 0.25
          cmap_limits = c(0.0, 0.6)
        } else if (this_nagents == 100) {
          this_mean_degree = 10
          this_min_group_frac = 0.2
          cmap_limits = c(0.0, 0.6)
        } else {
          this_mean_degree = 20
          this_min_group_frac = 0.05
          cmap_limits = c(0.0, 0.9)
        }
        
        this_tbl <- tbl %>%
          filter(group_w_innovation == this_group_w_innovation &
                 nagents == this_nagents &
                 mean_degree == this_mean_degree &
                 min_group_frac == this_min_group_frac)
        
        this_write_dir <- file.path(write_dir, "nagents", this_nagents)
        if (!dir.exists(this_write_dir)) { 
          dir.create(this_write_dir, recursive=TRUE)
        }

        write_path <- file.path(this_write_dir, paste0(this_group_w_innovation, ".pdf"))
        
        asymm_heatmap(this_tbl, this_group_w_innovation, write_path, 
                      cmap_limits = cmap_limits)
      }
    }

    # mean_degree sensitivity.
    if (do_mean_degree) {
      for (this_mean_degree in c(10, 30, 50)) {
        this_min_group_frac = 0.05
        this_nagents = 1000
        cmap_limits = c(0.0, 0.8)
        
        this_tbl <- tbl %>%
          filter(group_w_innovation == this_group_w_innovation &
                 nagents == this_nagents &
                 mean_degree == this_mean_degree &
                 min_group_frac == this_min_group_frac)

        this_write_dir <- file.path(write_dir, "mean_degree", this_mean_degree)
        if (!dir.exists(this_write_dir)) { 
          dir.create(this_write_dir, recursive=TRUE)
        }

        write_path <- file.path(this_write_dir, paste0(this_group_w_innovation, ".pdf"))
        
        asymm_heatmap(this_tbl, this_group_w_innovation, write_path, 
                      cmap_limits = cmap_limits)

      }
    }

    # min_group_frac 
    if (do_min_group_frac) {
      for (this_min_group_frac in c(0.2, 0.4)) {
        this_nagents = 1000
        this_mean_degree = 20
        cmap_limits = c(0.0, 0.65)
        
        this_tbl <- tbl %>%
          filter(group_w_innovation == this_group_w_innovation &
                 nagents == this_nagents &
                 mean_degree == this_mean_degree &
                 min_group_frac == this_min_group_frac)

        this_write_dir <- file.path(write_dir, "min_group_frac", this_min_group_frac)
        if (!dir.exists(this_write_dir)) { 
          dir.create(this_write_dir, recursive=TRUE)
        }

        write_path <- file.path(this_write_dir, paste0(this_group_w_innovation, ".pdf"))
        
        asymm_heatmap(this_tbl, this_group_w_innovation, write_path, 
                      cmap_limits = cmap_limits)

      }
    }
  }
    # minority group size sensitivity.
    # for (this_min_group_frac in c(0.2, 0.4)) { #, 0.5)) {
    #   
    #   this_tbl <- tbl[tbl$min_group_frac == this_min_group_frac, ]
    #   this_write_dir <- file.path(write_dir, "m", this_min_group_frac)
    #   write_path <- file.path(this_write_dir, paste0(group_w_innovation, ".pdf"))
    #   
    #   asymm_heatmap(this_tbl, group_w_innovation, write_path, cmap_limits = c(0.0, 0.8))
    # }
    # f(a) sensitivity.
    # for (this_a_fitness in c(1.05, 1.4)) { #, 2.0)) {
    #   
    #   this_tbl <- tbl[tbl$a_fitness == this_a_fitness, ]
    #   this_write_dir <- file.path(write_dir, "a_fitness", this_a_fitness)
    #   write_path <- file.path(this_write_dir, paste0(group_w_innovation, ".pdf"))
    #   asymm_heatmap(this_tbl, group_w_innovation, write_path, cmap_limits = c(0.0, 1.0))
    #   # if (this_a_fitness %in% c(1.4, 2.0)) {
    #     # asymm_heatmap(this_tbl, group_w_innovation, write_path, cmap_limits = c(0.0, 1.0))
    #   # }
    #   # else {
    #   #   asymm_heatmap(this_tbl, group_w_innovation, write_path, cmap_limits = c(0.0, 1.0))
    #   # }
    
  return (tbl)
}



asymm_heatmap <- function(asymm_tbl, this_group_w_innovation, write_path, 
                          measure = "sustainability", 
                          cmap_limits = c(0.0, 0.7), theme_extra = NULL, ...) {
  
  asymm_agg <- asymm_tbl %>%
    filter(min_homophily != 0.99) %>% 
    filter(maj_homophily != 0.99) %>%
    group_by(min_homophily, maj_homophily, group_w_innovation) %>%
    summarize(sustainability = mean(frac_a_curr_trait),
              step = mean(step))
  # return (asymm_agg)
    # asymm_tbl %>% 
    #   subset(group_w_innovation == this_group_w_innovation)
  asymm_lim_agg <- asymm_agg[asymm_agg$group_w_innovation == this_group_w_innovation, ]
  
  if (measure == "sustainability") {
  
    asymm_max_line <-
      asymm_lim_agg %>% 
        group_by(min_homophily) %>% 
        filter(sustainability == max(sustainability))
    
    max_sustainability <- 
      asymm_max_line[asymm_max_line$sustainability == 
                       max(asymm_max_line$sustainability), ]
  } else if (measure == "step") {
    asymm_max_line <-
      asymm_lim_agg %>% 
      group_by(min_homophily) %>% 
      filter(step == max(step))
    
    max_sustainability <- 
      asymm_max_line[asymm_max_line$step == 
                       max(asymm_max_line$step), ]
  } else {
    stop ("measure not recognized")
  }
  
  h1max <- max_sustainability$min_homophily
  h2max <- max_sustainability$maj_homophily
  
  if (measure == "sustainability") {
    ggplotstart <- ggplot(asymm_lim_agg, aes(x = min_homophily, y = maj_homophily, fill = sustainability))
    measure_label <- "Success\nrate"
  } else if (measure == "step") {
    ggplotstart <- ggplot(asymm_lim_agg, aes(x = min_homophily, y = maj_homophily, fill = step))
    measure_label <- "Mean\nsteps"
    cmap_limits = c(0, 60)
  } else {
    stop("Measure not recognized.")
  }
  
  p_ret <- ggplotstart + 
    geom_tile() +
    scale_fill_gradient2(low = "#000000", mid = "#010101", high = "#FFFFFF", limits = cmap_limits) +
    geom_point(data = asymm_max_line, aes(x = min_homophily, y = maj_homophily)) +
    geom_smooth(data = asymm_max_line, aes(x = min_homophily, y = maj_homophily), se=FALSE, n = 5) +
    geom_point(data = max_sustainability, aes(x=min_homophily, y=maj_homophily), 
               shape='diamond', size=5, color='red') +
    labs(x = TeX("Minority group homophily, $h_{min}$"), 
         y = TeX("Majority group homophily, $h_{maj}$")) +
    coord_fixed() + labs(fill = measure_label) +
    mytheme
    
  # save_path <- file.path(write_dir, str_replace(basename(csv_loc), ".csv", ".pdf"))
  if (!is_null(theme_extra)) { # add theme_extra to current p_ret }
    p_ret = p_ret + theme_extra
  }
  
  ggsave(write_path, p_ret, width = 6.75, height = 5)

  return (p_ret)
}


plot_group_freq_series <- function(csv_loc, write_dir = "figures/group_prevalence") {

    df <- read.csv(csv_loc)
    names(df) <- c("step", "frac_a", "Minority", "Majority", "Trial")
    df <- df[c("step", "Majority", "Minority", "Trial")]
    enslim <- 10
    df <- filter(df, Trial <= enslim)
    df$Trial = factor(df$Trial, levels = 1:enslim)
  
  df <- melt(df, id=c("step", "Trial"), value.name = "Frequency", variable.name = "Group")

  ggplot(df, aes(x=step, y=Frequency)) + 
    geom_line(aes(color=Trial, linetype=Group), lwd=0.8) +
    # geom_line(aes(x=step, y=frac_a_max, color=Trial, linetype="Majority"), linetype=1, lwd=1.05) +
    xlab("Step") + ylab(TeX(r"(Adaptation prevalence)")) +
    # scale_linetype_manual(name = "Group", values=c("Majority", "Minority"), labels = c("Majority", "Minority")) +
    mytheme #+ guides(color=guide_legend(override.aes=list(fill=NA))) 
  
  save_path <- file.path(write_dir, str_replace(basename(csv_loc), ".csv", ".pdf"))
  
  ggsave(save_path, width = 7.5, height = 4.65)

}


supp_mean_plots <- function(csv_dir = "data/supp_parts", 
                            write_dir = "figures/supp") {

  tbl <- load_from_parts(csv_dirs = c("data/supp_parts"), 
                         sync_file = "data/supp_mean_steps.csv")

  success_write_dir = file.path(write_dir, "success_jitter")
  steps_write_dir = file.path(write_dir, "mean_steps")

  for (this_nagents in c(50, 100, 200)) {
    this_tbl <- tbl[tbl$nagents == this_nagents, ]

    success_over_groups_jitter(NULL, success_write_dir, 
                               paste0("nagents=", this_nagents, ".pdf"),
                               tbl_override = this_tbl, 
                               this_ylim = c(0.0, 1.0)) 

    steps_ylim = c(0, 60)

    steps_over_groups_success_failure(NULL, steps_write_dir, 
                                      paste0("nagents=", this_nagents, ".pdf"),
                                      tbl_override = this_tbl,
                                      this_ylim = steps_ylim) 
  }

  # Now do the same for minority fraction and...
  for (this_min_group_frac in c(0.2, 0.35, 0.5)) {
    this_tbl <- tbl[tbl$min_group_frac == this_min_group_frac, ]

    success_over_groups_jitter(NULL, success_write_dir, 
                               paste0("min_group_frac=", this_min_group_frac, ".pdf"),
                               tbl_override = this_tbl,
                               this_ylim = c(0.0, 1.0)) 

    steps_ylim = c(0, 80)

    steps_over_groups_success_failure(NULL, steps_write_dir, 
                                      paste0("min_group_frac=", this_min_group_frac, ".pdf"),
                                      minority_pop_sizes = c(this_min_group_frac),
                                      tbl_override = this_tbl,
                                      this_ylim = steps_ylim) 
  }

  # ...adaptation fitness, f(a).
  for (this_a_fitness in c(1.05, 1.4, 2.0)) {
    this_tbl <- tbl[tbl$a_fitness == this_a_fitness, ]

    success_over_groups_jitter(NULL, success_write_dir, 
                               paste0("a_fitness=", this_a_fitness, ".pdf"),
                               tbl_override = this_tbl, 
                               this_ylim = c(0.0, 1.0)) 

    steps_ylim = c(0, 200)

    steps_over_groups_success_failure(NULL, steps_write_dir, 
                                      paste0("a_fitness=", this_a_fitness, ".pdf"),
                                      tbl_override = this_tbl,
                                      this_ylim = steps_ylim)
  }
}


plot_extreme_success_rates <- function(results_tbl, write_path = "../Writing/minmajnet/Figures/extreme_success_rates.pdf") {
  
  p <- results_tbl %>%
    # Process results table for summarizing and plotting over h_min.
    group_by(group_w_innovation, min_homophily, maj_homophily) %>%
    summarize(success_rate = mean(frac_a_curr_trait)) %>%
    group_by(group_w_innovation, min_homophily) %>%
    summarize(min_success = min(success_rate), max_success = max(success_rate)) %>%
    pivot_longer(cols = c(min_success, max_success), names_to = "variable", values_to = "value") %>%
    # Plot.
    ggplot(aes(x=min_homophily, linetype=variable, y = value, color = group_w_innovation)) +
      geom_line(size=2) + xlab(TeX('Minority homophily, $h_{min}$')) +  ylab("Success rate") + ylim(c(0, 0.8)) + 
      scale_color_manual(name = "Start group", values = c("#AA00AA", "#BBBBFF", "#00AAAA")) + 
      scale_linetype_discrete(name = "Type", labels = c("Max", "Min")) + mytheme
  
  ggsave(write_path, width=8, height=4.5)
    
  return (p)
}


success_reciprocity <- function(results_tbl, adj_mat_dir, write_path = "../Writing/minmajnet/Figures/extreme_success_rates.pdf") {
  
  success_over_hmaj <- results_tbl %>%
    group_by(group_w_)
    
}

calculate_reciprocity <- function(adj_mat_dir, hmin = 0.0, sync_file = "data/reciprocity.csv", lim = 0) {
  
  if (!is_null(sync_file)) {
    if (file.exists(sync_file)) {
      return (read_csv(sync_file))
    }
  }
  
  files <- Sys.glob(file.path(adj_mat_dir, paste0("adjacency_hmin=", format(hmin, nsmall=1), "*.csv")))
  if (lim > 0) {
    files <- files[1:lim]
  }
  nfiles <- length(files)
  hmin <- rep(hmin, nfiles)
  # hmaj <- as.numeric( str_split( str_split(str_extract(files, "hmaj=.*"), "=")[[1]][2], "_")[[1]][1] )
  
  hmaj <- rep(0.0, nfiles)
  clustering <- rep(0.0, nfiles)
  reciprocity_vec <- rep(0.0, nfiles)
  
  for (f_idx in 1:nfiles) {
    adjacency_matrix <- load_adjacency(files[f_idx])
    hmaj[f_idx] <- as.numeric( str_split( str_split(str_extract(files[f_idx], "hmaj=.*"), "=")[[1]][2], "_")[[1]][1] )
    
    clustering[f_idx] <- ClustBCG(adjacency_matrix, type="directed")$GlobaltotalCC
    
    network <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
    reciprocity_vec[f_idx] <- reciprocity(network)
  }
  
  return (tibble(hmin, hmaj, clustering, reciprocity = reciprocity_vec))
  # return (tibble(hmin, hmaj, reciprocity = reciprocity_vec))
}
