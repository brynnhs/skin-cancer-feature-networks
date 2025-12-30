# Network metrics

library(here)
library(yaml)
library(data.table)
library(igraph)

setwd(here())
config <- read_yaml("config/config.yaml")
dir.create("results/mgm", showWarnings = FALSE, recursive = TRUE)

adj_benign_df <- fread("results/mgm/adjacency_consensus_benign.csv", header = TRUE)
adj_benign <- as.matrix(adj_benign_df[, -1])
rownames(adj_benign) <- adj_benign_df[[1]]
colnames(adj_benign) <- names(adj_benign_df)[-1]

adj_malignant_df <- fread("results/mgm/adjacency_consensus_malignant.csv", header = TRUE)
adj_malignant <- as.matrix(adj_malignant_df[, -1])
rownames(adj_malignant) <- adj_malignant_df[[1]]
colnames(adj_malignant) <- names(adj_malignant_df)[-1]

n_edges_benign <- sum(adj_benign != 0) / 2
n_edges_malignant <- sum(adj_malignant != 0) / 2

compute_network_metrics <- function(adj_matrix, condition_name) {
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  
  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  
  metrics_df <- data.frame(
    node = V(g)$name,
    degree = degree(g),
    betweenness = betweenness(g, directed = FALSE, normalized = TRUE),
    closeness = if (n_edges > 0) closeness(g, normalized = TRUE) else rep(0, n_nodes),
    eigenvector = if (n_edges > 0) {
      ev <- eigen_centrality(g, directed = FALSE)$vector
      if (length(ev) == n_nodes) ev else rep(0, n_nodes)
    } else rep(0, n_nodes),
    stringsAsFactors = FALSE
  )
  
  if (n_edges > 0) {
    communities_walktrap <- tryCatch({
      walktrap.community(g)
    }, error = function(e) NULL)
    
    communities_louvain <- tryCatch({
      cluster_louvain(g)
    }, error = function(e) NULL)
    
    if (!is.null(communities_walktrap)) {
      metrics_df$community_walktrap <- membership(communities_walktrap)
    } else {
      metrics_df$community_walktrap <- rep(1, n_nodes)
    }
    
    if (!is.null(communities_louvain)) {
      metrics_df$community_louvain <- membership(communities_louvain)
    } else {
      metrics_df$community_louvain <- rep(1, n_nodes)
    }
  } else {
    metrics_df$community_walktrap <- rep(1, n_nodes)
    metrics_df$community_louvain <- rep(1, n_nodes)
  }
  
  metrics_df$condition <- condition_name
  
  return(list(metrics = metrics_df, graph = g))
}

metrics_benign <- compute_network_metrics(adj_benign, "benign")
metrics_malignant <- compute_network_metrics(adj_malignant, "malignant")

write.csv(metrics_benign$metrics, "results/mgm/centrality_benign.csv", row.names = FALSE)
write.csv(metrics_malignant$metrics, "results/mgm/centrality_malignant.csv", row.names = FALSE)

n_show_benign <- min(20, nrow(metrics_benign$metrics))
metrics_benign_sorted <- metrics_benign$metrics[order(metrics_benign$metrics$degree, decreasing = TRUE), ]
metrics_benign_plot <- metrics_benign_sorted[1:n_show_benign, ]

png("results/mgm/centrality_plot_benign.png", width = 1400, height = 800, res = 150)
par(mfrow = c(2, 2), mar = c(5, 4, 3, 2))
barplot(metrics_benign_plot$degree, names.arg = metrics_benign_plot$node, 
        las = 2, cex.names = 0.7, main = "Degree Centrality",
        ylab = "Degree")
barplot(metrics_benign_plot$betweenness, names.arg = metrics_benign_plot$node,
        las = 2, cex.names = 0.7, main = "Betweenness Centrality",
        ylab = "Betweenness")
barplot(metrics_benign_plot$closeness, names.arg = metrics_benign_plot$node,
        las = 2, cex.names = 0.7, main = "Closeness Centrality",
        ylab = "Closeness")
barplot(metrics_benign_plot$eigenvector, names.arg = metrics_benign_plot$node,
        las = 2, cex.names = 0.7, main = "Eigenvector Centrality",
        ylab = "Eigenvector")
mtext("Centrality Metrics: Benign", outer = TRUE, line = -1.5, cex = 1.2)
dev.off()

n_show_malignant <- min(20, nrow(metrics_malignant$metrics))
metrics_malignant_sorted <- metrics_malignant$metrics[order(metrics_malignant$metrics$degree, decreasing = TRUE), ]
metrics_malignant_plot <- metrics_malignant_sorted[1:n_show_malignant, ]

png("results/mgm/centrality_plot_malignant.png", width = 1400, height = 800, res = 150)
par(mfrow = c(2, 2), mar = c(5, 4, 3, 2))
barplot(metrics_malignant_plot$degree, names.arg = metrics_malignant_plot$node, 
        las = 2, cex.names = 0.7, main = "Degree Centrality",
        ylab = "Degree")
barplot(metrics_malignant_plot$betweenness, names.arg = metrics_malignant_plot$node,
        las = 2, cex.names = 0.7, main = "Betweenness Centrality",
        ylab = "Betweenness")
barplot(metrics_malignant_plot$closeness, names.arg = metrics_malignant_plot$node,
        las = 2, cex.names = 0.7, main = "Closeness Centrality",
        ylab = "Closeness")
barplot(metrics_malignant_plot$eigenvector, names.arg = metrics_malignant_plot$node,
        las = 2, cex.names = 0.7, main = "Eigenvector Centrality",
        ylab = "Eigenvector")
mtext("Centrality Metrics: Malignant", outer = TRUE, line = -1.5, cex = 1.2)
dev.off()
