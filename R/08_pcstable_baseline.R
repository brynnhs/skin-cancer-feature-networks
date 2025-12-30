# PC-Stable baseline

library(here)
library(yaml)
library(data.table)
library(dplyr)
library(pcalg)
library(igraph)

setwd(here())
config <- read_yaml("config/config.yaml")
set.seed(config$seeds$main)
dir.create("results/causal", showWarnings = FALSE, recursive = TRUE)

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

run_pcstable <- function(data_matrix, type_vector, condition_name) {
  n <- nrow(data_matrix)
  p <- ncol(data_matrix)
  alpha <- config$causal$alpha
  max_cond_set_size <- config$causal$max_cond_set_size
  
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
  
  n_skeleton_edges <- sum(skeleton_adj != 0) / 2
  n_oriented_edges <- sum(cpdag_adj != 0 & t(cpdag_adj) == 0)
  
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
    pc_result = pc_result,
    skeleton_adj = skeleton_adj,
    cpdag_adj = cpdag_adj,
    n_skeleton_edges = n_skeleton_edges,
    n_oriented_edges = n_oriented_edges,
    v_structures = v_structures
  ))
}

extract_cpdag_edges <- function(pc_result, condition_name) {
  if (is.null(pc_result)) {
    return(list(
      edges = data.frame(
        node1 = character(0), node2 = character(0),
        edge_type = character(0), condition = character(0),
        stringsAsFactors = FALSE
      ),
      vstruct = data.frame(
        X = character(0), Y = character(0), Z = character(0),
        condition = character(0), stringsAsFactors = FALSE
      )
    ))
  }
  
  edges_directed <- which(pc_result$cpdag_adj != 0 & t(pc_result$cpdag_adj) == 0, arr.ind = TRUE)
  edges_undirected <- which(pc_result$cpdag_adj != 0 & t(pc_result$cpdag_adj) != 0, arr.ind = TRUE)
  edges_undirected <- edges_undirected[edges_undirected[, 1] < edges_undirected[, 2], , drop = FALSE]
  
  edges_df <- data.frame(
    node1 = character(0), node2 = character(0),
    edge_type = character(0), condition = character(0),
    stringsAsFactors = FALSE
  )
  
  if (nrow(edges_directed) > 0) {
    directed_df <- data.frame(
      node1 = rownames(pc_result$cpdag_adj)[edges_directed[, 1]],
      node2 = colnames(pc_result$cpdag_adj)[edges_directed[, 2]],
      edge_type = "directed",
      condition = condition_name,
      stringsAsFactors = FALSE
    )
    edges_df <- rbind(edges_df, directed_df)
  }
  
  if (nrow(edges_undirected) > 0) {
    undirected_df <- data.frame(
      node1 = rownames(pc_result$cpdag_adj)[edges_undirected[, 1]],
      node2 = colnames(pc_result$cpdag_adj)[edges_undirected[, 2]],
      edge_type = "undirected",
      condition = condition_name,
      stringsAsFactors = FALSE
    )
    edges_df <- rbind(edges_df, undirected_df)
  }
  
  if (length(pc_result$v_structures) == 0) {
    vstruct_df <- data.frame(
      X = character(0), Y = character(0), Z = character(0),
      condition = character(0), stringsAsFactors = FALSE
    )
  } else {
    vstruct_df <- data.frame(
      X = character(length(pc_result$v_structures)),
      Y = character(length(pc_result$v_structures)),
      Z = character(length(pc_result$v_structures)),
      condition = condition_name,
      stringsAsFactors = FALSE
    )
    for (i in seq_along(pc_result$v_structures)) {
      vs <- pc_result$v_structures[[i]]
      vstruct_df$X[i] <- vs[1]
      vstruct_df$Y[i] <- vs[2]
      vstruct_df$Z[i] <- vs[3]
    }
  }
  
  list(edges = edges_df, vstruct = vstruct_df)
}

pc_benign <- run_pcstable(benign_prep$data, benign_prep$types, "Benign")
pc_malignant <- run_pcstable(malignant_prep$data, malignant_prep$types, "Malignant")

benign_results <- extract_cpdag_edges(pc_benign, "benign")
edges_benign <- benign_results$edges
vstruct_benign <- benign_results$vstruct

malignant_results <- extract_cpdag_edges(pc_malignant, "malignant")
edges_malignant <- malignant_results$edges
vstruct_malignant <- malignant_results$vstruct

if (!is.null(pc_benign)) {
  write.csv(pc_benign$skeleton_adj, "results/causal/skeleton_benign.csv", row.names = TRUE)
  write.csv(pc_benign$cpdag_adj, "results/causal/cpdag_benign.csv", row.names = TRUE)
  write.csv(edges_benign, "results/causal/edges_cpdag_benign.csv", row.names = FALSE)
  write.csv(vstruct_benign, "results/causal/vstructures_benign.csv", row.names = FALSE)
}

if (!is.null(pc_malignant)) {
  write.csv(pc_malignant$skeleton_adj, "results/causal/skeleton_malignant.csv", row.names = TRUE)
  write.csv(pc_malignant$cpdag_adj, "results/causal/cpdag_malignant.csv", row.names = TRUE)
  write.csv(edges_malignant, "results/causal/edges_cpdag_malignant.csv", row.names = FALSE)
  write.csv(vstruct_malignant, "results/causal/vstructures_malignant.csv", row.names = FALSE)
}

plot_cpdag <- function(pc_result, condition, out_path) {
  if (is.null(pc_result)) return()
  g <- graph_from_adjacency_matrix(pc_result$cpdag_adj, mode = "directed", weighted = NULL)
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
       main = paste("PC-Stable CPDAG:", condition))
  dev.off()
}

plot_cpdag(pc_benign, "Benign", "results/causal/cpdag_benign.png")
plot_cpdag(pc_malignant, "Malignant", "results/causal/cpdag_malignant.png")
