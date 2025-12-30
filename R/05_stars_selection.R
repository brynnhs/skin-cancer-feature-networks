# StARS selection

library(here)
library(yaml)
library(data.table)
library(dplyr)
library(mgm)
library(igraph)

setwd(here())
config <- read_yaml("config/config.yaml")
set.seed(config$seeds$stars)
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

lambda_min <- config$mgm$lambda_min
lambda_max <- config$mgm$lambda_max
n_lambda <- config$mgm$n_lambda
lambda_grid <- exp(seq(log(lambda_min), log(lambda_max), length.out = n_lambda))

run_stars <- function(data_matrix, type_vector, level_vector, condition_name) {
  n <- nrow(data_matrix)
  p <- ncol(data_matrix)
  n_subsamples <- config$stars$n_subsamples
  subsample_size <- round(n * config$stars$subsample_ratio)
  stability_threshold <- config$stars$stability_threshold
  
  instability_results <- data.frame(
    lambda = lambda_grid,
    instability = numeric(length(lambda_grid)),
    mean_edges = numeric(length(lambda_grid)),
    stringsAsFactors = FALSE
  )
  
  for (lambda_idx in seq_along(lambda_grid)) {
    lambda_val <- lambda_grid[lambda_idx]
    
    adj_list <- vector("list", n_subsamples)
    edge_counts <- numeric(n_subsamples)
    
    for (sub_idx in 1:n_subsamples) {
      subsample_indices <- sample(n, size = subsample_size, replace = FALSE)
      data_subsample <- data_matrix[subsample_indices, , drop = FALSE]
      
      adj_matrix <- matrix(0, nrow = p, ncol = p)
      
      mgm_result <- tryCatch({
        mgm(
          data = data_subsample,
          type = type_vector,
          level = level_vector,
          lambdaSeq = lambda_val,
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
        edge_counts[sub_idx] <- sum(adj_matrix != 0) / 2
      } else {
        edge_counts[sub_idx] <- 0
      }
      
      adj_list[[sub_idx]] <- adj_matrix
    }
    
    edge_selection_matrix <- matrix(0, nrow = p, ncol = p)
    for (sub_idx in 1:n_subsamples) {
      adj <- adj_list[[sub_idx]]
      edge_selection_matrix <- edge_selection_matrix + (adj != 0)
    }
    
    edge_frequency <- edge_selection_matrix / n_subsamples
    edge_variance <- edge_frequency * (1 - edge_frequency)
    edge_variance[lower.tri(edge_variance, diag = TRUE)] <- 0
    instability <- sum(edge_variance)
    
    max_edges <- p * (p - 1) / 2
    instability_normalized <- instability / max_edges
    
    instability_results$instability[lambda_idx] <- instability_normalized
    instability_results$mean_edges[lambda_idx] <- mean(edge_counts)
  }
  
  lambda_star_idx <- NA
  for (i in length(lambda_grid):1) {
    if (instability_results$instability[i] < stability_threshold) {
      lambda_star_idx <- i
      break
    }
  }
  
  if (is.na(lambda_star_idx)) {
    lambda_star_idx <- which.min(instability_results$instability)
  }
  
  lambda_star <- lambda_grid[lambda_star_idx]
  
  mgm_final <- tryCatch({
    mgm(
      data = data_matrix,
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
    if (lambda_star_idx < length(lambda_grid)) {
      lambda_star_alt <- lambda_grid[lambda_star_idx + 1]
      tryCatch({
        mgm(
          data = data_matrix,
          type = type_vector,
          level = level_vector,
          lambdaSeq = lambda_star_alt,
          lambdaSel = "EBIC",
          lambdaGam = 0.25,
          ruleReg = "OR",
          pbar = FALSE,
          k = 2,
          threshold = "none"
        )
      }, error = function(e2) {
        stop("Failed to refit MGM even with alternative lambda")
      })
    } else {
      stop("Failed to refit MGM with lambda*")
    }
  })
  
  adj_final <- mgm_final$pairwise$wadj
  rownames(adj_final) <- colnames(adj_final) <- valid_vars
  
  return(list(
    instability_results = instability_results,
    lambda_star = lambda_star,
    lambda_star_idx = lambda_star_idx,
    mgm_final = mgm_final,
    adj_final = adj_final
  ))
}

stars_benign <- run_stars(benign_prep$data, benign_prep$types, benign_prep$levels, "Benign")
stars_malignant <- run_stars(malignant_prep$data, malignant_prep$types, malignant_prep$levels, "Malignant")

adj_benign_upper <- stars_benign$adj_final
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
    node1 = rownames(stars_benign$adj_final)[edges_benign_idx[, 1]],
    node2 = colnames(stars_benign$adj_final)[edges_benign_idx[, 2]],
    weight = adj_benign_upper[edges_benign_idx],
    condition = "benign",
    stringsAsFactors = FALSE
  )
  edges_benign <- edges_benign[order(abs(edges_benign$weight), decreasing = TRUE), ]
}

adj_malignant_upper <- stars_malignant$adj_final
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
    node1 = rownames(stars_malignant$adj_final)[edges_malignant_idx[, 1]],
    node2 = colnames(stars_malignant$adj_final)[edges_malignant_idx[, 2]],
    weight = adj_malignant_upper[edges_malignant_idx],
    condition = "malignant",
    stringsAsFactors = FALSE
  )
  edges_malignant <- edges_malignant[order(abs(edges_malignant$weight), decreasing = TRUE), ]
}

write.csv(stars_benign$instability_results, "results/mgm/stars_curve_benign.csv", row.names = FALSE)
write.csv(stars_malignant$instability_results, "results/mgm/stars_curve_malignant.csv", row.names = FALSE)

lambda_star_summary <- data.frame(
  condition = c("benign", "malignant"),
  lambda_star = c(stars_benign$lambda_star, stars_malignant$lambda_star),
  instability = c(
    stars_benign$instability_results$instability[stars_benign$lambda_star_idx],
    stars_malignant$instability_results$instability[stars_malignant$lambda_star_idx]
  ),
  n_edges = c(nrow(edges_benign), nrow(edges_malignant))
)
write.csv(lambda_star_summary, "results/mgm/lambda_star_selected.csv", row.names = FALSE)

write.csv(stars_benign$adj_final, "results/mgm/adjacency_benign_stars.csv", row.names = TRUE)
write.csv(stars_malignant$adj_final, "results/mgm/adjacency_malignant_stars.csv", row.names = TRUE)
write.csv(edges_benign, "results/mgm/edges_benign_stars.csv", row.names = FALSE)
write.csv(edges_malignant, "results/mgm/edges_malignant_stars.csv", row.names = FALSE)

plot_stars_curve <- function(stars_result, condition, lambda_star, lambda_star_idx, out_path) {
  png(out_path, width = 2400, height = 1600, res = 300)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), cex.lab = 1.8, cex.axis = 1.6, cex.main = 2.0)
  plot(stars_result$instability_results$lambda, stars_result$instability_results$instability,
       type = "l", lwd = 3, col = "steelblue",
       xlab = "Lambda", ylab = "Instability",
       main = paste("StARS Curve:", condition),
       log = "x")
  abline(h = config$stars$stability_threshold, col = "red", lty = 2, lwd = 2)
  abline(v = lambda_star, col = "darkgreen", lty = 2, lwd = 2)
  points(lambda_star, 
         stars_result$instability_results$instability[lambda_star_idx],
         col = "darkgreen", pch = 19, cex = 2.5)
  legend("topright", 
         legend = c("Instability", "Threshold", "Lambda*"),
         col = c("steelblue", "red", "darkgreen"),
         lty = c(1, 2, 2), lwd = c(3, 2, 2), pch = c(NA, NA, 19),
         cex = 1.8)
  plot(stars_result$instability_results$lambda, stars_result$instability_results$mean_edges,
       type = "l", lwd = 3, col = "coral",
       xlab = "Lambda", ylab = "Mean Edges per Subsample",
       main = paste("Sparsity:", condition),
       log = "x")
  abline(v = lambda_star, col = "darkgreen", lty = 2, lwd = 2)
  dev.off()
  message("Saved: ", out_path)
}

plot_stars_curve(stars_benign, "Benign", stars_benign$lambda_star, stars_benign$lambda_star_idx, 
                 "results/mgm/stars_curve_benign.png")
plot_stars_curve(stars_malignant, "Malignant", stars_malignant$lambda_star, stars_malignant$lambda_star_idx,
                 "results/mgm/stars_curve_malignant.png")
