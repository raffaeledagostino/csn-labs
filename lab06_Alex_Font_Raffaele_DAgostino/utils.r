# Useful Functions 

jaccard_table <- function(memb1, memb2) {
  # memb1, memb2: vectors of cluster labels, same length, one label per node
  
  if (length(memb1) != length(memb2)) {
    stop("memb1 and memb2 must have the same length")
  }
  
  # unique cluster labels in each clustering
  cl1 <- sort(unique(memb1))
  cl2 <- sort(unique(memb2))
  
  # matrix to store Jaccard indices
  J <- matrix(0, 
              nrow = length(cl1), 
              ncol = length(cl2),
              dimnames = list(as.character(cl1), as.character(cl2)))
  
  # compute Jaccard index for every pair of clusters
  for (i in seq_along(cl1)) {
    A <- memb1 == cl1[i]          # nodes in cluster i of labeling 1
    for (j in seq_along(cl2)) {
      B <- memb2 == cl2[j]        # nodes in cluster j of labeling 2
      inter <- sum(A & B)         # |A ∩ B|
      uni   <- sum(A | B)         # |A ∪ B|
      J[i, j] <- if (uni == 0) 0 else inter / uni
    }
  }
  
  return(J)
}

##########################################################################################################

match_clusters <- function(JS, name1 = "L1", name2 = "L2") {
  # JS: jaccard matrix, rows = cluster labeling 1, columns = labeling 2

  # for each row, finds the column with highest Jaccard Index
  best_col_index <- apply(JS, 1, which.max)
  best_score     <- apply(JS, 1, max)
  
  # clusters' names in the two labelings
  cl1_names <- rownames(JS)
  cl2_names <- colnames(JS)[best_col_index]
  
  # building table
  result <- data.frame(
    cluster1   = cl1_names,
    cluster2   = cl2_names,
    jaccard    = best_score,
    stringsAsFactors = FALSE
  )
  
  result$cluster1 <- paste0(name1, ".", result$cluster1)
  result$cluster2 <- paste0(name2, ".", result$cluster2)
  
  return(result)
}


##########################################################################################################

# Function to evaluate clustering multiple times for randomized algorithms
evaluate_multiple_runs <- function(graph, alg_list, gt_clustering = NULL, n_runs = 10) {
  
  # Measures we want (these are ROWS in the result)
  measures <- c("clustering coef", "expansion", "conductance", "modularity")
  
  # Get algorithm names (these are COLUMNS in the result)
  alg_names <- names(alg_list)
  if (!is.null(gt_clustering)) {
    alg_names <- c(alg_names, "GT")
  }
  
  # Store results from each run
  all_results <- list()
  
  for (i in 1:n_runs) {
    res <- evaluate_significance(
      graph,
      alg_list = alg_list,
      gt_clustering = gt_clustering
    )
    
    # Convert to data frame and select measures (ROWS)
    df <- as.data.frame(res)
    df_subset <- df[measures, , drop = FALSE]  # Select ROWS, keep all columns
    all_results[[i]] <- df_subset
  }
  
  # Initialize result data frames
  mean_df <- all_results[[1]]  # Template with correct structure
  sd_df   <- all_results[[1]]
  
  # Compute mean and sd for each measure (row) and algorithm (column)
  for (measure in measures) {
    for (alg in colnames(mean_df)) {
      # Extract this cell across all runs
      values <- sapply(all_results, function(df) df[measure, alg])
      mean_df[measure, alg] <- mean(values, na.rm = TRUE)
      sd_df[measure, alg]   <- sd(values, na.rm = TRUE)
    }
  }
  
  return(list(mean = mean_df, sd = sd_df, all_runs = all_results))
}


##########################################################################################################

# Function to create nice summary table: mean ± sd
create_summary_table <- function(results, network_name) {
  mean_df <- results$mean
  sd_df <- results$sd
  
  summary_df <- mean_df
  for (i in 1:nrow(mean_df)) {
    for (j in 1:ncol(mean_df)) {
      m <- round(mean_df[i, j], 3)
      s <- round(sd_df[i, j], 3)
      summary_df[i, j] <- paste0(m, " ± ", s)
    }
  }
  
  cat("\n", network_name, "Summary (Mean ± SD):\n")
  print(summary_df)
  latex_table <- xtable(summary_df, type = "latex")
  print(latex_table, file = paste0(network_name, "Summary (Mean ± SD).tex"))
  return(summary_df)
}


##########################################################################################################

Wmean <- function(MC, memb_GT) {
  # MC: output from match_clusters (data frame with columns: cluster1, cluster2, jaccard)
  # memb_GT: membership vector of ground truth clustering
  
  # Extract Jaccard values from the data frame
  if (is.data.frame(MC)) {
    jaccard_values <- MC$jaccard
  } else if (is.numeric(MC)) {
    jaccard_values <- as.numeric(MC)
  } else {
    stop("MC must be a data frame from match_clusters or a numeric vector")
  }
  
  # Get unique GT clusters (in order)
  gt_clusters <- sort(unique(memb_GT))
  
  # Compute size (number of nodes) for each GT cluster
  cluster_sizes <- sapply(gt_clusters, function(c) sum(memb_GT == c))
  
  # Weights = fraction of total nodes in each cluster
  weights <- cluster_sizes / sum(cluster_sizes)
  
  # Check dimensions match
  if (length(jaccard_values) != length(weights)) {
    stop("Number of Jaccard values doesn't match number of GT clusters")
  }
  
  # Weighted mean
  weighted_mean <- sum(weights * jaccard_values)
  
  return(weighted_mean)
}

##########################################################################################################

Wmin <- function(MC, memb_GT) {
  # MC: output from match_clusters (data frame with columns: cluster1, cluster2, jaccard)
  # memb_GT: membership vector of ground truth clustering
  
  # Extract Jaccard values from the data frame
  if (is.data.frame(MC)) {
    jaccard_values <- MC$jaccard
  } else if (is.numeric(MC)) {
    jaccard_values <- as.numeric(MC)
  } else {
    stop("MC must be a data frame from match_clusters or a numeric vector")
  }
  
  # Get unique GT clusters (in order)
  gt_clusters <- sort(unique(memb_GT))
  
  # Compute size (number of nodes) for each GT cluster
  cluster_sizes <- sapply(gt_clusters, function(c) sum(memb_GT == c))
  
  # Weights = fraction of total nodes in each cluster
  weights <- cluster_sizes / sum(cluster_sizes)
  
  # Check dimensions match
  if (length(jaccard_values) != length(weights)) {
    stop("Number of Jaccard values doesn't match number of GT clusters")
  }
  
  min_jaccard = min(jaccard_values)
  
  # Find which cluster(s) have the minimum Jaccard
  min_idx <- which(jaccard_values == min_jaccard)
  
  threshold <- 0.01
  relevant_indices <- which(weights >= threshold)
  
  # Find minimum Jaccard index
  min_jaccard <- min(jaccard_values[relevant_indices])
  
  
  return(list(
    min_jaccard = min_jaccard,
    worst_cluster_idx = min_idx
  ))
}


##########################################################################################################

Wharmonic <- function(MC, memb_GT) {
  # MC: output from match_clusters (data frame with columns: cluster1, cluster2, jaccard)
  # memb_GT: membership vector of ground truth clustering
  
  # Extract Jaccard values from the data frame
  if (is.data.frame(MC)) {
    jaccard_values <- MC$jaccard
  } else if (is.numeric(MC)) {
    jaccard_values <- as.numeric(MC)
  } else {
    stop("MC must be a data frame from match_clusters or a numeric vector")
  }
  
  # Get unique GT clusters (in order)
  gt_clusters <- sort(unique(memb_GT))
  
  # Compute size (number of nodes) for each GT cluster
  cluster_sizes <- sapply(gt_clusters, function(c) sum(memb_GT == c))
  
  # Weights = fraction of total nodes in each cluster
  weights <- cluster_sizes / sum(cluster_sizes)
  
  # Check dimensions match
  if (length(jaccard_values) != length(weights)) {
    stop("Number of Jaccard values doesn't match number of GT clusters")
  }
  
  # Check for zero Jaccard values (would cause division by zero)
  if (any(jaccard_values == 0)) {
    warning("Some Jaccard values are 0; harmonic mean will be 0")
    return(0)
  }
  
  # Weighted harmonic mean: 1 / sum(weights / values)
  weighted_harmonic <- 1 / sum(weights / jaccard_values)
  
  return(weighted_harmonic)
}


##########################################################################################################



##########################################################################################################



##########################################################################################################



##########################################################################################################



##########################################################################################################



##########################################################################################################
