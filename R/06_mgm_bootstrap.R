# MGM bootstrap

library(here)
library(yaml)
library(data.table)
library(dplyr)
library(mgm)
library(igraph)

setwd(here())
config <- read_yaml("config/config.yaml")
set.seed(config$seeds$bootstrap)
dir.create("results/mgm", showWarnings = FALSE, recursive = TRUE)

benign_data <- fread("data/processed/benign.csv", na.strings = c("", "NA", "NaN", "null", "NULL"))
malignant_data <- fread("data/processed/malignant.csv", na.strings = c("", "NA", "NaN", "null", "NULL"))
benign_data <- as.data.frame(benign_data)
malignant_data <- as.data.frame(malignant_data)

lambda_star_df <- fread("results/mgm/lambda_star_selected.csv")
lambda_star_benign <- lambda_star_df[condition == "benign", lambda_star]
lambda_star_malignant <- lambda_star_df[condition == "malignant", lambda_star]

id_cols <- config$data$id_col
diagnosis_col <- config$data$diagnosis_col
exclude_cols <- c(id_cols, diagnosis_col)
analysis_vars <- setdiff(names(benign_data), exclude_cols)

prepare_mgm_data <- function(df, analysis_vars) {
  data_matrix <- df[, analysis_vars, drop = FALSE]
  
  type_vector <- character(length(analysis_vars))
  names(type_vector) <- analysis_vars
  
  for (var in analysis_vars) {
    var_class <- class(data_matrix[[var]])[1]
    if (var_class %in% c("numeric", "integer")) {
      n_unique <- length(unique(data_matrix[[var]][!is.na(data_matrix[[var]])]))
      type_vector[var] <- ifelse(n_unique <= 10 && n_unique > 2, "c", "g")
    } else if (var_class == "logical") {
      type_vector[var] <- "c"
    } else if (var_class == "factor") {
      type_vector[var] <- ifelse(is.ordered(data_matrix[[var]]), "g", "c")
    } else {
      type_vector[var] <- "c"
    }
  }
  
  data_cleaned <- data_matrix
  for (var in analysis_vars) {
    if (type_vector[var] == "c") {
      cat_counts <- table(data_cleaned[[var]], useNA = "no")
      rare_cats <- names(cat_counts)[cat_counts < 2]
      if (length(rare_cats) > 0) {
        most_common <- names(cat_counts)[which.max(cat_counts)]
        if (is.factor(data_cleaned[[var]])) {
          levels(data_cleaned[[var]]) <- c(levels(data_cleaned[[var]]), most_common)
        }
        data_cleaned[[var]][data_cleaned[[var]] %in% rare_cats] <- most_common
        if (is.factor(data_cleaned[[var]])) {
          data_cleaned[[var]] <- droplevels(data_cleaned[[var]])
        }
      }
    }
  }
  
  data_numeric <- data_cleaned[, analysis_vars, drop = FALSE]
  for (var in analysis_vars) {
    if (is.factor(data_numeric[[var]])) {
      data_numeric[[var]] <- as.numeric(data_numeric[[var]])
    } else if (is.logical(data_numeric[[var]])) {
      data_numeric[[var]] <- as.numeric(data_numeric[[var]])
    } else if (is.character(data_numeric[[var]])) {
      data_numeric[[var]] <- as.numeric(as.factor(data_numeric[[var]]))
    }
    data_numeric[[var]] <- as.numeric(data_numeric[[var]])
  }
  
  data_matrix_numeric <- as.matrix(data_numeric)
  type_vector <- type_vector[analysis_vars]
  
  level_vector <- integer(length(analysis_vars))
  names(level_vector) <- analysis_vars
  for (var in analysis_vars) {
    if (type_vector[var] == "c") {
      level_vector[var] <- length(unique(data_cleaned[[var]][!is.na(data_cleaned[[var]])]))
    } else {
      level_vector[var] <- 1
    }
  }
  
  var_variance <- apply(data_matrix_numeric, 2, function(x) stats::var(x, na.rm = TRUE))
  zero_var_vars <- names(var_variance)[var_variance == 0 | is.na(var_variance)]
  
  if (length(zero_var_vars) > 0) {
    valid_vars <- setdiff(analysis_vars, zero_var_vars)
    data_matrix_numeric <- data_matrix_numeric[, valid_vars, drop = FALSE]
    type_vector <- type_vector[valid_vars]
    level_vector <- level_vector[valid_vars]
  } else {
    valid_vars <- analysis_vars
  }
  
  return(list(data = data_matrix_numeric, types = type_vector, levels = level_vector, vars = valid_vars))
}

benign_prep <- prepare_mgm_data(benign_data, analysis_vars)
malignant_prep <- prepare_mgm_data(malignant_data, analysis_vars)

valid_vars <- intersect(benign_prep$vars, malignant_prep$vars)
benign_prep$data <- benign_prep$data[, valid_vars, drop = FALSE]
benign_prep$types <- benign_prep$types[valid_vars]
benign_prep$levels <- benign_prep$levels[valid_vars]
malignant_prep$data <- malignant_prep$data[, valid_vars, drop = FALSE]
malignant_prep$types <- malignant_prep$types[valid_vars]
malignant_prep$levels <- malignant_prep$levels[valid_vars]

run_bootstrap <- function(data_matrix, type_vector, level_vector, lambda_star, condition_name) {
  n <- nrow(data_matrix)
  p <- ncol(data_matrix)
  n_boot <- config$bootstrap$n_boot
  consensus_threshold <- config$bootstrap$consensus_threshold
  
  edge_selection_counts <- matrix(0, nrow = p, ncol = p)
  rownames(edge_selection_counts) <- colnames(edge_selection_counts) <- valid_vars
  
  successful_fits <- 0
  for (boot_idx in 1:n_boot) {
    bootstrap_indices <- sample(n, size = n, replace = TRUE)
    data_bootstrap <- data_matrix[bootstrap_indices, , drop = FALSE]
    
    mgm_result <- tryCatch({
      mgm(
        data = data_bootstrap,
        type = type_vector,
        level = level_vector,
        lambdaSeq = lambda_star,
        lambdaSel = "EBIC",
        lambdaGam = 0.25,
        ruleReg = "OR",
        pbar = FALSE,
        k = 2,
        threshold = "none"
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(mgm_result)) {
      adj_matrix <- mgm_result$pairwise$wadj
      edge_selection_counts <- edge_selection_counts + (adj_matrix != 0)
      successful_fits <- successful_fits + 1
    }
  }
  
  edge_frequencies <- edge_selection_counts / successful_fits
  
  consensus_adj <- (edge_frequencies >= consensus_threshold) * 1
  consensus_adj[lower.tri(consensus_adj)] <- t(consensus_adj)[lower.tri(consensus_adj)]
  diag(consensus_adj) <- 0
  
  consensus_weights <- edge_frequencies
  consensus_weights[consensus_adj == 0] <- 0
  
  return(list(
    edge_frequencies = edge_frequencies,
    consensus_adj = consensus_adj,
    consensus_weights = consensus_weights,
    successful_fits = successful_fits,
    total_boot = n_boot
  ))
}

bootstrap_benign <- run_bootstrap(
  benign_prep$data,
  benign_prep$types,
  benign_prep$levels,
  lambda_star_benign,
  "Benign"
)

bootstrap_malignant <- run_bootstrap(
  malignant_prep$data,
  malignant_prep$types,
  malignant_prep$levels,
  lambda_star_malignant,
  "Malignant"
)

extract_edge_list <- function(adj_matrix, condition_name) {
  adj_upper <- adj_matrix
  adj_upper[lower.tri(adj_upper, diag = TRUE)] <- 0
  edge_idx <- which(adj_upper > 0, arr.ind = TRUE)
  
  if (nrow(edge_idx) == 0) {
    return(data.frame(
      node1 = character(0),
      node2 = character(0),
      frequency = numeric(0),
      condition = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  edge_df <- data.frame(
    node1 = rownames(adj_matrix)[edge_idx[, 1]],
    node2 = colnames(adj_matrix)[edge_idx[, 2]],
    frequency = adj_upper[edge_idx],
    condition = condition_name,
    stringsAsFactors = FALSE
  )
  edge_df[order(edge_df$frequency, decreasing = TRUE), ]
}

edge_freq_benign <- extract_edge_list(bootstrap_benign$edge_frequencies, "benign")
edge_freq_malignant <- extract_edge_list(bootstrap_malignant$edge_frequencies, "malignant")
consensus_edges_benign <- extract_edge_list(bootstrap_benign$consensus_weights, "benign")
consensus_edges_malignant <- extract_edge_list(bootstrap_malignant$consensus_weights, "malignant")

write.csv(edge_freq_benign, "results/mgm/edge_frequencies_benign.csv", row.names = FALSE)
write.csv(edge_freq_malignant, "results/mgm/edge_frequencies_malignant.csv", row.names = FALSE)
write.csv(bootstrap_benign$consensus_adj, "results/mgm/adjacency_consensus_benign.csv", row.names = TRUE)
write.csv(bootstrap_malignant$consensus_adj, "results/mgm/adjacency_consensus_malignant.csv", row.names = TRUE)
write.csv(consensus_edges_benign, "results/mgm/edges_consensus_benign.csv", row.names = FALSE)
write.csv(consensus_edges_malignant, "results/mgm/edges_consensus_malignant.csv", row.names = FALSE)

g_benign <- graph_from_adjacency_matrix(bootstrap_benign$consensus_adj, mode = "undirected", weighted = TRUE)
g_malignant <- graph_from_adjacency_matrix(bootstrap_malignant$consensus_adj, mode = "undirected", weighted = TRUE)

all_nodes <- sort(unique(c(V(g_benign)$name, V(g_malignant)$name)))
combined_adj <- pmax(bootstrap_benign$consensus_adj, bootstrap_malignant$consensus_adj)
g_combined <- graph_from_adjacency_matrix(combined_adj, mode = "undirected")
V(g_combined)$name <- all_nodes

if (vcount(g_combined) > 0) {
  set.seed(42)
  layout_shared <- layout_with_kk(g_combined)
  layout_shared <- layout_shared - apply(layout_shared, 2, min)
  layout_range <- apply(layout_shared, 2, max) - apply(layout_shared, 2, min)
  layout_range[layout_range == 0] <- 1
  layout_shared <- (layout_shared / layout_range) * 2 - 1
  rownames(layout_shared) <- all_nodes
} else {
  layout_shared <- matrix(0, nrow = length(all_nodes), ncol = 2)
  rownames(layout_shared) <- all_nodes
}

assign_edge_weights <- function(g, consensus_weights, consensus_adj) {
  if (ecount(g) > 0) {
    edge_list <- as_edgelist(g)
    edge_weights <- numeric(ecount(g))
    weights_upper <- consensus_weights
    weights_upper[lower.tri(weights_upper, diag = TRUE)] <- 0
    for (i in 1:ecount(g)) {
      node1_idx <- which(rownames(consensus_adj) == edge_list[i, 1])
      node2_idx <- which(colnames(consensus_adj) == edge_list[i, 2])
      edge_weights[i] <- weights_upper[node1_idx, node2_idx]
    }
    E(g)$weight <- edge_weights
  }
  g
}

g_benign <- assign_edge_weights(g_benign, bootstrap_benign$consensus_weights, bootstrap_benign$consensus_adj)
g_malignant <- assign_edge_weights(g_malignant, bootstrap_malignant$consensus_weights, bootstrap_malignant$consensus_adj)

layout_benign <- layout_shared[V(g_benign)$name, , drop = FALSE]
layout_malignant <- layout_shared[V(g_malignant)$name, , drop = FALSE]

plot_consensus_network <- function(g, layout, condition, out_path) {
  png(out_path, width = 2400, height = 2400, res = 300)
  par(mar = c(2, 2, 4, 2), cex.main = 2.2)
  if (ecount(g) > 0) {
    plot(g,
         layout = layout,
         vertex.label = V(g)$name,
         vertex.size = 18,
         vertex.color = "lightblue",
         vertex.label.cex = 1.5,
         vertex.label.dist = 0.6,
         vertex.label.color = "black",
         edge.width = E(g)$weight * 6,
         edge.color = "gray40",
         main = paste("Consensus Network:", condition, "\n(Threshold:", config$bootstrap$consensus_threshold, ")"))
  } else {
    plot(g,
         layout = layout,
         vertex.label = V(g)$name,
         vertex.size = 18,
         vertex.color = "lightblue",
         vertex.label.cex = 1.5,
         vertex.label.dist = 0.6,
         vertex.label.color = "black",
         main = paste("Consensus Network:", condition, "\n(Threshold:", config$bootstrap$consensus_threshold, " - No edges)"))
  }
  dev.off()
  message("Saved: ", out_path)
}

plot_consensus_network(g_benign, layout_benign, "Benign", "results/mgm/network_consensus_benign.png")
plot_consensus_network(g_malignant, layout_malignant, "Malignant", "results/mgm/network_consensus_malignant.png")
