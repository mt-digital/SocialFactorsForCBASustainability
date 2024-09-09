library(igraph)
library(purrrlyr)
##
# Create minority-majority asymmetric homophily networks 
#
BASE_EDGE_COLORS = c(adjustcolor("turquoise4", alpha.f = 0.8),
                     adjustcolor("darkgreen", alpha.f = 0.5), 
                     adjustcolor("goldenrod3", alpha.f = 0.5),
                     adjustcolor("pink", alpha.f = 0.215)
                    )
plot_network = function(net, min_frac, edge_color_options = BASE_EDGE_COLORS, 
                        coords = NULL) {
   

  metapop_size = length(net)
  N_min = round(min_frac * metapop_size)
  N_maj = metapop_size - N_min

  group_vec = c(rep(1, N_min), rep(2, N_maj))

  V(net)$color = c("turquoise4", "goldenrod3")[group_vec]


  N_min = sum(group_vec == 1)
  N_maj = sum(group_vec == 2)

  N = length(group_vec)

  if (is.null(coords)) {
    coords = matrix(NA, nrow = length(group_vec), ncol = 2)
    coords = matrix(0.0, nrow = N, ncol = 2)

    coords[group_vec == 1, 1] = -7.0
    coords[group_vec == 2, 1] = 7.0
    
    coords[group_vec == 1] = 
      coords[group_vec == 1] + 
      cbind(rnorm(N_min, sd = .75), rnorm(N_min, sd = 1.5))

    coords[group_vec == 2] = 
      coords[group_vec == 2] + 
      cbind(rnorm(N_maj, sd = 1.25), rnorm(N_maj, sd = 1.25))
  }

  edges_df = as_data_frame(net)
  
  edge_color_lookup = function(row) {
    color_code = NA
    if (group_vec[row$from] == 1 & group_vec[row$to] == 1) {
      color_code = 1
    } else if (group_vec[row$from] == 2 & group_vec[row$to] == 1) {
      color_code = 2
    } else if (group_vec[row$from] == 2 & group_vec[row$to] == 2) {
      color_code = 3
    } else {
      color_code = 4
    }

    return (edge_color_options[color_code])
  }

  edge_color = by_row(as_data_frame(net), edge_color_lookup, .collate = "rows")$.out

  p = plot(net, edge.arrow.size = 0.5, edge.width = 1.5, 
       edge.curved = 0.4, edge.color = edge_color, 
       xlim = c(-9, 11), ylim = c(-1.0, 1.0), vertex.label=NA, 
       vertex.size = 35, layout = coords, rescale = FALSE
      )

  return (coords)
}


# Make a minority-majority network using construction from Turner, Reynolds, and
# Jones.
make_min_maj_net = function(metapop_size = 200, min_frac = 0.05, 
                            mean_degree = 10,
                            homophily_min = 0.0,
                            homophily_maj = 0.0) {

  n_edges = mean_degree * metapop_size
  E = n_edges

  E_min = round(n_edges * min_frac)
  E_maj = E - E_min

  E_min_teach_min = round(E_min * ((1 + homophily_min) / 2.0))
  E_maj_teach_min = E_min - E_min_teach_min

  E_maj_teach_maj = round(E_maj * ((1 + homophily_maj) / 2.0))
  E_min_teach_maj = E_maj - E_maj_teach_maj
  
  min_size = round(min_frac * metapop_size)
  maj_size = metapop_size - min_size
  group = c(rep(1, min_size), rep(2, maj_size))
  
  network <- make_graph(c(), n = metapop_size, directed = TRUE)

  minority_ids <- 1:min_size
  majority_ids <- (min_size + 1):metapop_size

  for (xx in 1:E_min_teach_min) {

    # Generate a new edge that may already exist in the social network.
    resample = TRUE
    while (resample) {
      new_edge = sample(minority_ids, 2)
      resample = are_adjacent(network, new_edge[1], new_edge[2])
    }

    network = add_edges(network, new_edge)
  }

  for (xx in 1:E_maj_teach_maj) {
    
    resample = TRUE
    while (resample) {
      new_edge = sample(majority_ids, 2)
      resample = are_adjacent(network, new_edge[1], new_edge[2])
    }

    network = add_edges(network, new_edge)
  }

  for (xx in 1:E_min_teach_maj) {
    
    resample = TRUE
    while (resample) {
      new_edge = c(sample(minority_ids, 1), sample(majority_ids, 1))
      resample = are_adjacent(network, new_edge[1], new_edge[2])
    }

    network = add_edges(network, new_edge)
  }

  for (xx in 1:E_maj_teach_min) {
    
    resample = TRUE
    while (resample) {
      new_edge = c(sample(majority_ids, 1), sample(minority_ids, 1))
      resample = are_adjacent(network, new_edge[1], new_edge[2])
    }

    network = add_edges(network, new_edge)
  }

  return (network)
}


visualize_minmaj_collection = function(metapop_size = 200, min_frac = 0.05,
                                       mean_degree = 10,
                                       h_min_maj_pairs = 
                                         rbind(c(0.0, 0.0),
                                               c(-0.4, 0.7),
                                               c(0.7, 0.7)),
                                       savedir = 
                                         file.path("figures", "minmaj_networks"),
                                       seed = NULL
                                      ) {

  if (!dir.exists(savedir))
    dir.create(savedir, recursive = TRUE)

  for (pair_idx in 1:nrow(h_min_maj_pairs)) {

    h_min = h_min_maj_pairs[pair_idx, 1]
    h_maj = h_min_maj_pairs[pair_idx, 2]

    net = make_min_maj_net(metapop_size, min_frac, mean_degree, h_min, h_maj)

    savefile = file.path(
      savedir, paste0("h_min=", h_min, "_", "h_maj=", h_maj, ".pdf")
    )

    pdf(savefile, 7.5, 5)
      plot_network(net, min_frac)
    dev.off()
  }

}


visualize_network_creation = 
  function(metapop_size = 200, min_frac = 0.05, mean_degree = 10,
           h_min = 0.0, h_maj = 0.0,
           savedir = 
             file.path("figures", "minmaj_networks", "creation"),
           coords = NULL
          ) {
  
  if (!dir.exists(savedir))
    dir.create(savedir, recursive = TRUE)

  net = make_min_maj_net(metapop_size, min_frac, mean_degree, h_min, h_maj)

  for (condition in c("no-edge", "min-min", "maj-min", "maj-maj", "min-maj")) {

    if (condition == "no-edge") {
      edge_colors = rep(adjustcolor("white", alpha.f = 0.0), 4)
    } else if (condition == "min-min") {
      edge_colors = c(BASE_EDGE_COLORS[1], 
                      rep(adjustcolor("white", alpha.f = 0.0), 3))
    } else if (condition == "maj-min") {
      edge_colors = c(BASE_EDGE_COLORS[1], BASE_EDGE_COLORS[2], 
                      rep(adjustcolor("white", alpha.f = 0.0), 2))
    } else if (condition == "maj-maj") {
      edge_colors = c(BASE_EDGE_COLORS[1], BASE_EDGE_COLORS[2], 
                      BASE_EDGE_COLORS[3], adjustcolor("white", alpha.f = 0.0))
    } else {
      edge_colors = BASE_EDGE_COLORS
    }

    savefile = file.path(
        savedir, paste0("h_min=", h_min, "_", 
                        "h_maj=", h_maj, "_", 
                        condition, ".pdf")
    )

    pdf(savefile, 7.5, 5)

    coords = plot_network(net, min_frac, edge_color_options = edge_colors, coords)

    dev.off()
  }
}
