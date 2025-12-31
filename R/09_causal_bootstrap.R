# Causal bootstrap

library(here)
library(yaml)
library(data.table)
library(dplyr)
library(pcalg)
library(igraph)

setwd(here())
config <- read_yaml("config/config.yaml")
set.seed(config$seeds$bootstrap)
dir.create("results/causal_bootstrap", showWarnings = FALSE, recursive = TRUE)

benign_data <- fread("data/processed/benign.csv", na.strings = c("", "NA", "NaN", "null", "NULL"))
malignant_data <- fread("data/processed/malignant.csv", na.strings = c("", "NA", "NaN", "null", "NULL"))
benign_data <- as.data.frame(benign_data)
malignant_data <- as.data.frame(malignant_data)

id_cols <- config$data$id_col
diagnosis_col <- config$data$diagnosis_col
exclude_cols <- c(id_cols, diagnosis_col)
analysis_vars <- setdiff(names(benign_data), exclude_cols)

prepare_pc_data <- function(df, analysis_vars) {
  data_subset <- df[, analysis_vars, drop = FALSE]
  
  type_vector <- character(length(analysis_vars))
  names(type_vector) <- analysis_vars
  
  for (var in analysis_vars) {
    var_class <- class(data_subset[[var]])[1]
    if (var_class %in% c("numeric", "integer")) {
      n_unique <- length(unique(data_subset[[var]][!is.na(data_subset[[var]])]))
      type_vector[var] <- ifelse(n_unique <= 10 && n_unique > 2, "discrete", "continuous")
    } else if (var_class == "logical") {
      type_vector[var] <- "discrete"
    } else if (var_class == "factor") {
      type_vector[var] <- "discrete"
    } else {
      type_vector[var] <- "discrete"
    }
  }
  
  data_numeric <- data_subset[, analysis_vars, drop = FALSE]
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
  
  data_matrix <- as.matrix(data_numeric)
  
  var_variance <- apply(data_matrix, 2, function(x) stats::var(x, na.rm = TRUE))
  zero_var_vars <- names(var_variance)[var_variance == 0 | is.na(var_variance)]
  
  if (length(zero_var_vars) > 0) {
    valid_vars <- setdiff(analysis_vars, zero_var_vars)
    data_matrix <- data_matrix[, valid_vars, drop = FALSE]
    type_vector <- type_vector[valid_vars]
  } else {
    valid_vars <- analysis_vars
  }
  
  return(list(data = data_matrix, types = type_vector, vars = valid_vars))
}

benign_prep <- prepare_pc_data(benign_data, analysis_vars)
malignant_prep <- prepare_pc_data(malignant_data, analysis_vars)

valid_vars <- intersect(benign_prep$vars, malignant_prep$vars)
benign_prep$data <- benign_prep$data[, valid_vars, drop = FALSE]
benign_prep$types <- benign_prep$types[valid_vars]
malignant_prep$data <- malignant_prep$data[, valid_vars, drop = FALSE]
malignant_prep$types <- malignant_prep$types[valid_vars]

run_pcstable_single <- function(data_matrix, valid_vars, alpha) {
  n <- nrow(data_matrix)
  p <- ncol(data_matrix)
  
  suffStat <- list(C = cor(data_matrix, use = "complete.obs"), n = n)
  
  pc_result <- tryCatch({
    pc(suffStat = suffStat,
       indepTest = gaussCItest,
       alpha = alpha,
       labels = valid_vars,
       u2pd = "relaxed",
       skel.method = "stable",
       maj.rule = TRUE,
       solve.confl = TRUE,
       numCores = 1,
       verbose = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(pc_result)) {
    return(NULL)
  }
  
  skeleton_adj <- as(pc_result@graph, "matrix")
  skeleton_adj <- (skeleton_adj != 0) * 1
  
  cpdag_adj <- as(pc_result@graph, "matrix")
  
  v_structures <- list()
  for (y in 1:p) {
    parents <- which(cpdag_adj[, y] != 0 & cpdag_adj[y, ] == 0)
    if (length(parents) >= 2) {
      for (i in 1:(length(parents) - 1)) {
        for (j in (i + 1):length(parents)) {
          x <- parents[i]
          z <- parents[j]
          if (cpdag_adj[x, z] == 0 && cpdag_adj[z, x] == 0) {
            v_structures[[length(v_structures) + 1]] <- c(valid_vars[x], valid_vars[y], valid_vars[z])
          }
        }
      }
    }
  }
  
  return(list(
    skeleton_adj = skeleton_adj,
    cpdag_adj = cpdag_adj,
    v_structures = v_structures
  ))
}

run_causal_bootstrap <- function(data_matrix, valid_vars, condition_name) {
  n <- nrow(data_matrix)
  p <- ncol(data_matrix)
  n_boot <- config$causal$bootstrap$n_boot
  consensus_threshold <- config$causal$bootstrap$consensus_threshold
  alpha <- config$causal$alpha
  
  skeleton_edge_counts <- matrix(0, nrow = p, ncol = p)
  rownames(skeleton_edge_counts) <- colnames(skeleton_edge_counts) <- valid_vars
  
  directed_edge_counts <- matrix(0, nrow = p, ncol = p)
  rownames(directed_edge_counts) <- colnames(directed_edge_counts) <- valid_vars
  
  vstructure_counts <- list()
  successful_fits <- 0
  
  for (boot_idx in 1:n_boot) {
    bootstrap_indices <- sample(n, size = n, replace = TRUE)
    data_bootstrap <- data_matrix[bootstrap_indices, , drop = FALSE]
    
    pc_result <- run_pcstable_single(data_bootstrap, valid_vars, alpha)
    
    if (!is.null(pc_result)) {
      skeleton_edge_counts <- skeleton_edge_counts + (pc_result$skeleton_adj != 0)
      
      directed_edges <- (pc_result$cpdag_adj != 0) & (t(pc_result$cpdag_adj) == 0)
      directed_edge_counts <- directed_edge_counts + directed_edges
      
      for (vs in pc_result$v_structures) {
        vs_key <- paste(sort(c(vs[1], vs[3])), collapse = "|")
        vs_full <- paste(vs, collapse = "->")
        if (is.null(vstructure_counts[[vs_key]])) {
          vstructure_counts[[vs_key]] <- list()
        }
        if (is.null(vstructure_counts[[vs_key]][[vs_full]])) {
          vstructure_counts[[vs_key]][[vs_full]] <- 0
        }
        vstructure_counts[[vs_key]][[vs_full]] <- vstructure_counts[[vs_key]][[vs_full]] + 1
      }
      
      successful_fits <- successful_fits + 1
    }
  }
  
  skeleton_frequencies <- skeleton_edge_counts / successful_fits
  directed_frequencies <- directed_edge_counts / successful_fits
  
  consensus_cpdag <- matrix(0, nrow = p, ncol = p)
  rownames(consensus_cpdag) <- colnames(consensus_cpdag) <- valid_vars
  
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j) {
        if (skeleton_frequencies[i, j] >= consensus_threshold) {
          freq_i_to_j <- directed_frequencies[i, j]
          freq_j_to_i <- directed_frequencies[j, i]
          
          if (freq_i_to_j >= consensus_threshold) {
            consensus_cpdag[i, j] <- 1
          } else if (freq_j_to_i >= consensus_threshold) {
            consensus_cpdag[j, i] <- 1
          } else {
            consensus_cpdag[i, j] <- 1
            consensus_cpdag[j, i] <- 1
          }
        }
      }
    }
  }
  
  return(list(
    skeleton_frequencies = skeleton_frequencies,
    directed_frequencies = directed_frequencies,
    vstructure_counts = vstructure_counts,
    consensus_cpdag = consensus_cpdag,
    successful_fits = successful_fits,
    total_boot = n_boot
  ))
}

bootstrap_benign <- run_causal_bootstrap(benign_prep$data, valid_vars, "Benign")
bootstrap_malignant <- run_causal_bootstrap(malignant_prep$data, valid_vars, "Malignant")

freq_upper_benign <- bootstrap_benign$skeleton_frequencies
freq_upper_benign[lower.tri(freq_upper_benign, diag = TRUE)] <- 0
edges_benign_idx <- which(freq_upper_benign > 0, arr.ind = TRUE)

if (nrow(edges_benign_idx) == 0) {
  edge_freq_benign <- data.frame(
    node1 = character(0),
    node2 = character(0),
    frequency = numeric(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
} else {
  edge_freq_benign <- data.frame(
    node1 = rownames(bootstrap_benign$skeleton_frequencies)[edges_benign_idx[, 1]],
    node2 = colnames(bootstrap_benign$skeleton_frequencies)[edges_benign_idx[, 2]],
    frequency = freq_upper_benign[edges_benign_idx],
    condition = "benign",
    stringsAsFactors = FALSE
  )
  edge_freq_benign <- edge_freq_benign[order(edge_freq_benign$frequency, decreasing = TRUE), ]
}

freq_upper_malignant <- bootstrap_malignant$skeleton_frequencies
freq_upper_malignant[lower.tri(freq_upper_malignant, diag = TRUE)] <- 0
edges_malignant_idx <- which(freq_upper_malignant > 0, arr.ind = TRUE)

if (nrow(edges_malignant_idx) == 0) {
  edge_freq_malignant <- data.frame(
    node1 = character(0),
    node2 = character(0),
    frequency = numeric(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
} else {
  edge_freq_malignant <- data.frame(
    node1 = rownames(bootstrap_malignant$skeleton_frequencies)[edges_malignant_idx[, 1]],
    node2 = colnames(bootstrap_malignant$skeleton_frequencies)[edges_malignant_idx[, 2]],
    frequency = freq_upper_malignant[edges_malignant_idx],
    condition = "malignant",
    stringsAsFactors = FALSE
  )
  edge_freq_malignant <- edge_freq_malignant[order(edge_freq_malignant$frequency, decreasing = TRUE), ]
}

directed_edges_benign_idx <- which(bootstrap_benign$directed_frequencies > 0, arr.ind = TRUE)

if (nrow(directed_edges_benign_idx) == 0) {
  directed_freq_benign <- data.frame(
    node1 = character(0),
    node2 = character(0),
    frequency = numeric(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
} else {
  directed_freq_benign <- data.frame(
    node1 = rownames(bootstrap_benign$directed_frequencies)[directed_edges_benign_idx[, 1]],
    node2 = colnames(bootstrap_benign$directed_frequencies)[directed_edges_benign_idx[, 2]],
    frequency = bootstrap_benign$directed_frequencies[directed_edges_benign_idx],
    condition = "benign",
    stringsAsFactors = FALSE
  )
  directed_freq_benign <- directed_freq_benign[order(directed_freq_benign$frequency, decreasing = TRUE), ]
}

directed_edges_malignant_idx <- which(bootstrap_malignant$directed_frequencies > 0, arr.ind = TRUE)

if (nrow(directed_edges_malignant_idx) == 0) {
  directed_freq_malignant <- data.frame(
    node1 = character(0),
    node2 = character(0),
    frequency = numeric(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
} else {
  directed_freq_malignant <- data.frame(
    node1 = rownames(bootstrap_malignant$directed_frequencies)[directed_edges_malignant_idx[, 1]],
    node2 = colnames(bootstrap_malignant$directed_frequencies)[directed_edges_malignant_idx[, 2]],
    frequency = bootstrap_malignant$directed_frequencies[directed_edges_malignant_idx],
    condition = "malignant",
    stringsAsFactors = FALSE
  )
  directed_freq_malignant <- directed_freq_malignant[order(directed_freq_malignant$frequency, decreasing = TRUE), ]
}

vstruct_freq_benign <- NULL
for (vs_key in names(bootstrap_benign$vstructure_counts)) {
  for (vs_full in names(bootstrap_benign$vstructure_counts[[vs_key]])) {
    count <- bootstrap_benign$vstructure_counts[[vs_key]][[vs_full]]
    frequency <- count / bootstrap_benign$successful_fits
    
    parts <- strsplit(vs_full, "->")[[1]]
    if (length(parts) == 3) {
      new_row <- data.frame(
        X = parts[1],
        Y = parts[2],
        Z = parts[3],
        frequency = frequency,
        condition = "benign",
        stringsAsFactors = FALSE
      )
      if (is.null(vstruct_freq_benign)) {
        vstruct_freq_benign <- new_row
      } else {
        vstruct_freq_benign <- rbind(vstruct_freq_benign, new_row)
      }
    }
  }
}

if (is.null(vstruct_freq_benign)) {
  vstruct_freq_benign <- data.frame(
    X = character(0),
    Y = character(0),
    Z = character(0),
    frequency = numeric(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
} else {
  vstruct_freq_benign <- vstruct_freq_benign[order(vstruct_freq_benign$frequency, decreasing = TRUE), ]
}

vstruct_freq_malignant <- NULL
for (vs_key in names(bootstrap_malignant$vstructure_counts)) {
  for (vs_full in names(bootstrap_malignant$vstructure_counts[[vs_key]])) {
    count <- bootstrap_malignant$vstructure_counts[[vs_key]][[vs_full]]
    frequency <- count / bootstrap_malignant$successful_fits
    
    parts <- strsplit(vs_full, "->")[[1]]
    if (length(parts) == 3) {
      new_row <- data.frame(
        X = parts[1],
        Y = parts[2],
        Z = parts[3],
        frequency = frequency,
        condition = "malignant",
        stringsAsFactors = FALSE
      )
      if (is.null(vstruct_freq_malignant)) {
        vstruct_freq_malignant <- new_row
      } else {
        vstruct_freq_malignant <- rbind(vstruct_freq_malignant, new_row)
      }
    }
  }
}

if (is.null(vstruct_freq_malignant)) {
  vstruct_freq_malignant <- data.frame(
    X = character(0),
    Y = character(0),
    Z = character(0),
    frequency = numeric(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
} else {
  vstruct_freq_malignant <- vstruct_freq_malignant[order(vstruct_freq_malignant$frequency, decreasing = TRUE), ]
}

edges_directed_consensus_benign <- which(bootstrap_benign$consensus_cpdag != 0 & t(bootstrap_benign$consensus_cpdag) == 0, arr.ind = TRUE)
edges_undirected_consensus_benign <- which(bootstrap_benign$consensus_cpdag != 0 & t(bootstrap_benign$consensus_cpdag) != 0, arr.ind = TRUE)
edges_undirected_consensus_benign <- edges_undirected_consensus_benign[edges_undirected_consensus_benign[, 1] < edges_undirected_consensus_benign[, 2], , drop = FALSE]

consensus_edges_benign <- data.frame(
  node1 = character(0),
  node2 = character(0),
  edge_type = character(0),
  condition = character(0),
  stringsAsFactors = FALSE
)

if (nrow(edges_directed_consensus_benign) > 0) {
  directed_consensus_benign <- data.frame(
    node1 = rownames(bootstrap_benign$consensus_cpdag)[edges_directed_consensus_benign[, 1]],
    node2 = colnames(bootstrap_benign$consensus_cpdag)[edges_directed_consensus_benign[, 2]],
    edge_type = "directed",
    condition = "benign",
    stringsAsFactors = FALSE
  )
  consensus_edges_benign <- rbind(consensus_edges_benign, directed_consensus_benign)
}

if (nrow(edges_undirected_consensus_benign) > 0) {
  undirected_consensus_benign <- data.frame(
    node1 = rownames(bootstrap_benign$consensus_cpdag)[edges_undirected_consensus_benign[, 1]],
    node2 = colnames(bootstrap_benign$consensus_cpdag)[edges_undirected_consensus_benign[, 2]],
    edge_type = "undirected",
    condition = "benign",
    stringsAsFactors = FALSE
  )
  consensus_edges_benign <- rbind(consensus_edges_benign, undirected_consensus_benign)
}

edges_directed_consensus_malignant <- which(bootstrap_malignant$consensus_cpdag != 0 & t(bootstrap_malignant$consensus_cpdag) == 0, arr.ind = TRUE)
edges_undirected_consensus_malignant <- which(bootstrap_malignant$consensus_cpdag != 0 & t(bootstrap_malignant$consensus_cpdag) != 0, arr.ind = TRUE)
edges_undirected_consensus_malignant <- edges_undirected_consensus_malignant[edges_undirected_consensus_malignant[, 1] < edges_undirected_consensus_malignant[, 2], , drop = FALSE]

consensus_edges_malignant <- data.frame(
  node1 = character(0),
  node2 = character(0),
  edge_type = character(0),
  condition = character(0),
  stringsAsFactors = FALSE
)

if (nrow(edges_directed_consensus_malignant) > 0) {
  directed_consensus_malignant <- data.frame(
    node1 = rownames(bootstrap_malignant$consensus_cpdag)[edges_directed_consensus_malignant[, 1]],
    node2 = colnames(bootstrap_malignant$consensus_cpdag)[edges_directed_consensus_malignant[, 2]],
    edge_type = "directed",
    condition = "malignant",
    stringsAsFactors = FALSE
  )
  consensus_edges_malignant <- rbind(consensus_edges_malignant, directed_consensus_malignant)
}

if (nrow(edges_undirected_consensus_malignant) > 0) {
  undirected_consensus_malignant <- data.frame(
    node1 = rownames(bootstrap_malignant$consensus_cpdag)[edges_undirected_consensus_malignant[, 1]],
    node2 = colnames(bootstrap_malignant$consensus_cpdag)[edges_undirected_consensus_malignant[, 2]],
    edge_type = "undirected",
    condition = "malignant",
    stringsAsFactors = FALSE
  )
  consensus_edges_malignant <- rbind(consensus_edges_malignant, undirected_consensus_malignant)
}

write.csv(edge_freq_benign, "results/causal_bootstrap/edge_frequencies_benign.csv", row.names = FALSE)
write.csv(edge_freq_malignant, "results/causal_bootstrap/edge_frequencies_malignant.csv", row.names = FALSE)
write.csv(directed_freq_benign, "results/causal_bootstrap/directed_edge_frequencies_benign.csv", row.names = FALSE)
write.csv(directed_freq_malignant, "results/causal_bootstrap/directed_edge_frequencies_malignant.csv", row.names = FALSE)
write.csv(vstruct_freq_benign, "results/causal_bootstrap/vstructure_frequencies_benign.csv", row.names = FALSE)
write.csv(vstruct_freq_malignant, "results/causal_bootstrap/vstructure_frequencies_malignant.csv", row.names = FALSE)
write.csv(bootstrap_benign$consensus_cpdag, "results/causal_bootstrap/cpdag_consensus_benign.csv", row.names = TRUE)
write.csv(bootstrap_malignant$consensus_cpdag, "results/causal_bootstrap/cpdag_consensus_malignant.csv", row.names = TRUE)
write.csv(consensus_edges_benign, "results/causal_bootstrap/edges_cpdag_consensus_benign.csv", row.names = FALSE)
write.csv(consensus_edges_malignant, "results/causal_bootstrap/edges_cpdag_consensus_malignant.csv", row.names = FALSE)

plot_consensus_cpdag <- function(consensus_cpdag, condition, out_path) {
  g <- graph_from_adjacency_matrix(consensus_cpdag, mode = "directed", weighted = NULL)
  layout_pos <- layout_with_fr(g)
  png(out_path, width = 1200, height = 1200, res = 150)
  par(mar = c(0, 0, 2, 0))
  plot(g,
       layout = layout_pos,
       vertex.label = V(g)$name,
       vertex.size = 8,
       vertex.color = "lightblue",
       vertex.label.cex = 0.8,
       edge.arrow.size = 0.5,
       edge.color = "gray50",
       main = paste("Consensus CPDAG:", condition, "\n(Threshold:", config$causal$bootstrap$consensus_threshold, ")"))
  dev.off()
}

plot_consensus_cpdag(bootstrap_benign$consensus_cpdag, "Benign", "results/causal_bootstrap/cpdag_consensus_benign.png")
plot_consensus_cpdag(bootstrap_malignant$consensus_cpdag, "Malignant", "results/causal_bootstrap/cpdag_consensus_malignant.png")

baseline_edges_benign <- tryCatch({
  read.csv("results/causal/edges_cpdag_benign.csv")
}, error = function(e) {
  data.frame(node1 = character(0), node2 = character(0), edge_type = character(0))
})

baseline_edges_malignant <- tryCatch({
  read.csv("results/causal/edges_cpdag_malignant.csv")
}, error = function(e) {
  data.frame(node1 = character(0), node2 = character(0), edge_type = character(0))
})
