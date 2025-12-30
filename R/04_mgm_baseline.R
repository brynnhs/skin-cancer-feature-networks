# MGM baseline

library(here)
library(yaml)
library(data.table)
library(dplyr)
library(mgm)
library(igraph)

setwd(here())
config <- read_yaml("config/config.yaml")
set.seed(config$seeds$main)
dir.create("results/mgm", showWarnings = FALSE, recursive = TRUE)

benign_data <- fread("data/processed/benign.csv", na.strings = c("", "NA", "NaN", "null", "NULL"))
malignant_data <- fread("data/processed/malignant.csv", na.strings = c("", "NA", "NaN", "null", "NULL"))
benign_data <- as.data.frame(benign_data)
malignant_data <- as.data.frame(malignant_data)

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
      if (n_unique <= 10 && n_unique > 2) {
        type_vector[var] <- "c"
      } else {
        type_vector[var] <- "g"
      }
    } else if (var_class == "logical") {
      type_vector[var] <- "c"
    } else if (var_class == "factor") {
      if (is.ordered(data_matrix[[var]])) {
        type_vector[var] <- "g"
      } else {
        type_vector[var] <- "c"
      }
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
      n_cats <- length(unique(data_cleaned[[var]][!is.na(data_cleaned[[var]])]))
      level_vector[var] <- n_cats
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

lambda_fixed <- (config$mgm$lambda_min + config$mgm$lambda_max) / 2

mgm_benign <- mgm(
  data = benign_prep$data,
  type = benign_prep$types,
  level = benign_prep$levels,
  lambdaSeq = lambda_fixed,
  lambdaSel = "EBIC",
  lambdaGam = 0.25,
  ruleReg = "OR",
  pbar = FALSE,
  k = 2,
  threshold = "none"
)

adj_benign <- mgm_benign$pairwise$wadj
rownames(adj_benign) <- colnames(adj_benign) <- valid_vars

mgm_malignant <- mgm(
  data = malignant_prep$data,
  type = malignant_prep$types,
  level = malignant_prep$levels,
  lambdaSeq = lambda_fixed,
  lambdaSel = "EBIC",
  lambdaGam = 0.25,
  ruleReg = "OR",
  pbar = FALSE,
  k = 2,
  threshold = "none"
)

adj_malignant <- mgm_malignant$pairwise$wadj
rownames(adj_malignant) <- colnames(adj_malignant) <- valid_vars

adj_benign_upper <- adj_benign
adj_benign_upper[lower.tri(adj_benign_upper, diag = TRUE)] <- 0
edges_benign_idx <- which(adj_benign_upper != 0, arr.ind = TRUE)

if (nrow(edges_benign_idx) == 0) {
  edges_benign <- data.frame(
    node1 = character(0),
    node2 = character(0),
    weight = numeric(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
} else {
  edges_benign <- data.frame(
    node1 = rownames(adj_benign)[edges_benign_idx[, 1]],
    node2 = colnames(adj_benign)[edges_benign_idx[, 2]],
    weight = adj_benign_upper[edges_benign_idx],
    condition = "benign",
    stringsAsFactors = FALSE
  )
  edges_benign <- edges_benign[order(abs(edges_benign$weight), decreasing = TRUE), ]
}

adj_malignant_upper <- adj_malignant
adj_malignant_upper[lower.tri(adj_malignant_upper, diag = TRUE)] <- 0
edges_malignant_idx <- which(adj_malignant_upper != 0, arr.ind = TRUE)

if (nrow(edges_malignant_idx) == 0) {
  edges_malignant <- data.frame(
    node1 = character(0),
    node2 = character(0),
    weight = numeric(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
} else {
  edges_malignant <- data.frame(
    node1 = rownames(adj_malignant)[edges_malignant_idx[, 1]],
    node2 = colnames(adj_malignant)[edges_malignant_idx[, 2]],
    weight = adj_malignant_upper[edges_malignant_idx],
    condition = "malignant",
    stringsAsFactors = FALSE
  )
  edges_malignant <- edges_malignant[order(abs(edges_malignant$weight), decreasing = TRUE), ]
}

write.csv(adj_benign, "results/mgm/adjacency_benign.csv", row.names = TRUE)
write.csv(adj_malignant, "results/mgm/adjacency_malignant.csv", row.names = TRUE)
write.csv(edges_benign, "results/mgm/edges_benign.csv", row.names = FALSE)
write.csv(edges_malignant, "results/mgm/edges_malignant.csv", row.names = FALSE)

g_benign <- graph_from_adjacency_matrix(adj_benign, mode = "undirected", weighted = TRUE, diag = FALSE)
g_benign <- delete.vertices(g_benign, degree(g_benign) == 0)

if (vcount(g_benign) > 0) {
  png("results/mgm/network_benign.png", width = 1200, height = 1200, res = 150)
  layout_benign <- layout_with_fr(g_benign)
  plot(g_benign,
       layout = layout_benign,
       vertex.label.cex = 0.7,
       vertex.size = 8,
       vertex.color = "lightblue",
       edge.width = abs(E(g_benign)$weight) * 3,
       edge.color = ifelse(E(g_benign)$weight > 0, "darkgreen", "red"),
       main = "MGM Network: Benign")
  dev.off()
}

g_malignant <- graph_from_adjacency_matrix(adj_malignant, mode = "undirected", weighted = TRUE, diag = FALSE)
g_malignant <- delete.vertices(g_malignant, degree(g_malignant) == 0)

if (vcount(g_malignant) > 0) {
  png("results/mgm/network_malignant.png", width = 1200, height = 1200, res = 150)
  layout_malignant <- layout_with_fr(g_malignant)
  plot(g_malignant,
       layout = layout_malignant,
       vertex.label.cex = 0.7,
       vertex.size = 8,
       vertex.color = "lightblue",
       edge.width = abs(E(g_malignant)$weight) * 3,
       edge.color = ifelse(E(g_malignant)$weight > 0, "darkgreen", "red"),
       main = "MGM Network: Malignant")
  dev.off()
}
