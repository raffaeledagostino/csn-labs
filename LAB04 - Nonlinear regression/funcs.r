

###################################################################################################################
create.dependency.tree.metrics <- function() {
  base_dir <- "./pud"
  files <- list.files(path = base_dir,
                      pattern = "pud-ud-test\\.conllu$",
                      recursive = TRUE,
                      full.names = TRUE)
  
  if (!dir.exists("./data")) {
    dir.create("./data")
  }
  
  for (f in files) {
    lines <- readLines(f, encoding = "UTF-8")
    
    # split blocks by sentence
    sentence_blocks <- split(lines, cumsum(grepl("^# sent_id", lines)))
    
    # parse single sentence
    parse_sentence <- function(block) {
      token_lines <- grep("^[0-9]+\\t", block, value = TRUE)
      if (length(token_lines) == 0) return(NULL)
      
      df <- read.table(text = token_lines, sep = "\t", header = FALSE, 
                       quote = "",  
                       fill = TRUE, 
                       stringsAsFactors = FALSE,
                       comment.char = "")  
      
      if (ncol(df) < 8) return(NULL)
      
      colnames(df)[1:8] <- c("ID", "FORM", "LEMMA", "UPOS", "XPOS", 
                             "FEATS", "HEAD", "DEPREL")
      
      # composed token removal
      df <- df[!grepl("-", df$ID), ]
      
      df$ID <- as.numeric(df$ID)
      df$HEAD <- as.numeric(df$HEAD)
      df <- df[!is.na(df$HEAD), ]
      
      return(df)
    }
    
    sentences <- lapply(sentence_blocks, parse_sentence)
    sentences <- sentences[!sapply(sentences, is.null)]
    
    parent_folder <- basename(dirname(f))
    lang_name <- sub("^UD_(.*)-PUD$", "\\1", parent_folder)
    
    metrics_list <- lapply(sentences, function(df) {
      n <- nrow(df)
      
      degrees <- sapply(df$ID, function(id) {
        num_children <- sum(df$HEAD == id)
        has_parent <- ifelse(df$HEAD[df$ID == id] != 0, 1, 0)
        return(num_children + has_parent)
      })
      
      k2 <- sum(degrees^2) / n

      distances <- abs(df$HEAD - df$ID)
      distances <- distances[df$HEAD != 0]
      
      if (length(distances) > 0) {
        d_mean <- mean(distances)
      } else {
        d_mean <- 0.0  
      }
      
      return(c(as.integer(n), as.double(k2), as.double(d_mean)))
    })
    
    metrics_matrix <- do.call(rbind, metrics_list)
    
    out_file <- paste0("./data/", lang_name, "_dependency_tree_metrics.txt")
    
    write.table(metrics_matrix, file = out_file, 
                row.names = FALSE, col.names = FALSE, 
                sep = "\t", quote = FALSE)
    
    cat("Created:", out_file, "\n")
  }
}

###################################################################################################################

validate.dependency.metrics <- function() {
  
  data_dir <- "./data"
  metric_files <- list.files(path = data_dir,
                             pattern = "_dependency_tree_metrics\\.txt$",
                             full.names = TRUE)
  
  base_dir <- "./pud"
  pud_files <- list.files(path = base_dir,
                          pattern = "pud-ud-test\\.conllu$",
                          recursive = TRUE,
                          full.names = TRUE)
  
  # Parse a single sentence 
  parse_sentence <- function(block) {
    token_lines <- grep("^[0-9]+\\t", block, value = TRUE)
    if (length(token_lines) == 0) return(NULL)
    
    df <- read.table(text = token_lines, sep = "\t", header = FALSE, 
                     quote = "", fill = TRUE, stringsAsFactors = FALSE,
                     comment.char = "")
    
    colnames(df)[1:8] <- c("ID", "FORM", "LEMMA", "UPOS", "XPOS", 
                           "FEATS", "HEAD", "DEPREL")
    
    df <- df[!grepl("-", df$ID), ]
    df$ID <- as.numeric(df$ID)
    df$HEAD <- as.numeric(df$HEAD)
    
    return(df)
  }
  
  # Calculate q 
  calculate_q <- function(df) {
    degrees <- sapply(df$ID, function(id) {
      num_children <- sum(df$HEAD == id)
      has_parent <- ifelse(df$HEAD[df$ID == id] != 0, 1, 0)
      return(num_children + has_parent)
    })
    
    q <- sum(degrees %% 2)
    return(q)
  }
  
  summary_table <- data.frame()
  languages_with_violations <- c()
  
  for (pud_file in pud_files) {
    parent_folder <- basename(dirname(pud_file))
    lang_name <- sub("^UD_(.*)-PUD$", "\\1", parent_folder)
    
    metric_file <- file.path(data_dir, paste0(lang_name, "_dependency_tree_metrics.txt"))
    
    if (!file.exists(metric_file)) next
    
    metrics <- read.table(metric_file, header = FALSE, sep = "\t", 
                          colClasses = c("integer", "numeric", "numeric"))
    colnames(metrics) <- c("n", "k2", "d_mean") 
    
    lines <- readLines(pud_file, encoding = "UTF-8")
    sentence_blocks <- split(lines, cumsum(grepl("^# sent_id", lines)))
    sentences <- lapply(sentence_blocks, parse_sentence)
    sentences <- sentences[!sapply(sentences, is.null)]
    
    q_values <- sapply(sentences, calculate_q)
    metrics$q <- q_values[1:nrow(metrics)]
    
    metrics <- metrics %>%
      mutate(
        k2_lower_bound = 4 - 6/n,
        k2_upper_bound = n - 1,
        d_lower_bound = (1/4) * (n * k2 + q),
        d_upper_bound = (1/4) * (3 * (n - 1)^2 + 1 - (n %% 2)),
        
        # Round to 3 decimal places
        k2 = round(k2, digits = 3),
        k2_lower_bound = round(k2_lower_bound, digits = 3),
        k2_upper_bound = round(k2_upper_bound, digits = 3),
        d_lower_bound = round(d_lower_bound, digits = 3),
        d_upper_bound = round(d_upper_bound, digits = 3),
        
        # Check validity
        k2_valid = (k2 >= k2_lower_bound) & (k2 <= k2_upper_bound),
        d_valid = (d_mean * (n-1) >= d_lower_bound - 1e-6) & (d_mean * (n-1) <= d_upper_bound + 1e-6),
        all_valid = k2_valid & d_valid
      )
    
    # Print invalid rows for this language
    invalid_rows <- metrics[!metrics$all_valid, ]
    if (nrow(invalid_rows) > 0) {
      cat("\n========================================\n")
      cat(sprintf("VIOLATIONS IN: %s\n", lang_name))
      cat("========================================\n\n")
      
      for (i in 1:nrow(invalid_rows)) {
        row <- invalid_rows[i, ]
        sent_idx <- which(metrics$all_valid == row$all_valid & 
                            metrics$n == row$n & 
                            metrics$d_mean == row$d_mean)[1]
        
        cat(sprintf("Sentence %d (n=%d, q=%d):\n", sent_idx, row$n, row$q))
        
        if (!row$k2_valid) {
          cat(sprintf("  k² VIOLATION: %.3f not in [%.3f, %.3f]\n", 
                      row$k2, row$k2_lower_bound, row$k2_upper_bound))
        }
        
        if (!row$d_valid) {
          cat(sprintf("  <d> VIOLATION:\n"))
          cat(sprintf("    d_mean = %.3f\n", row$d_mean))
          cat(sprintf("    d_total = d_mean × (n-1) = %.3f × %d = %.3f\n", 
                      row$d_mean, row$n - 1, row$d_mean * (row$n - 1)))
          cat(sprintf("    Bounds: [%.3f, %.3f]\n", 
                      row$d_lower_bound, row$d_upper_bound))
        }
        cat("\n")
      }
      
      languages_with_violations <- c(languages_with_violations, lang_name)
    }
    
    n_sentences <- nrow(metrics)
    n_k2_valid <- sum(metrics$k2_valid, na.rm = TRUE)
    n_d_valid <- sum(metrics$d_valid, na.rm = TRUE)
    n_all_valid <- sum(metrics$all_valid, na.rm = TRUE)
    
    summary_table <- rbind(summary_table, data.frame(
      Language = lang_name,
      k2_perc = 100 * n_k2_valid / n_sentences,
      d_perc = 100 * n_d_valid / n_sentences,
      both_perc = 100 * n_all_valid / n_sentences,
      stringsAsFactors = FALSE
    ))
  }
  
  # Print results
  cat("\n========================================\n")
  cat("VALIDATION SUMMARY\n")
  cat("========================================\n\n")
  
  if (length(languages_with_violations) == 0) {
    cat("All metrics are valid for all languages.\n\n")
  } else {
    cat("Languages with violations:\n")
    for (lang in languages_with_violations) {
      lang_row <- summary_table[summary_table$Language == lang, ]
      cat(sprintf("  - %s: %.1f%% valid\n", lang, lang_row$both_perc))
    }
    cat("\n")
  }
  
  # Print summary table
  summary_table$k2_perc <- sprintf("%.1f%%", summary_table$k2_perc)
  summary_table$d_perc <- sprintf("%.1f%%", summary_table$d_perc)
  summary_table$both_perc <- sprintf("%.1f%%", summary_table$both_perc)
  colnames(summary_table) <- c("Language", "k² valid", "<d> valid", "Both valid")
  print(summary_table, row.names = FALSE)
  cat("\n")
}

####################################################################################################################
validate.dependency.metrics2 <- function() {
  
  data_dir <- "./data"
  metric_files <- list.files(path = data_dir,
                             pattern = "_dependency_tree_metrics\\.txt$",
                             full.names = TRUE)
  
  base_dir <- "./pud"
  pud_files <- list.files(path = base_dir,
                          pattern = "pud-ud-test\\.conllu$",
                          recursive = TRUE,
                          full.names = TRUE)
  
  parse_sentence <- function(block) {
    token_lines <- grep("^[0-9]+\\t", block, value = TRUE)
    if (length(token_lines) == 0) return(NULL)
    
    df <- read.table(text = token_lines, sep = "\t", header = FALSE, 
                     quote = "", fill = TRUE, stringsAsFactors = FALSE,
                     comment.char = "")
    
    colnames(df)[1:8] <- c("ID", "FORM", "LEMMA", "UPOS", "XPOS", 
                           "FEATS", "HEAD", "DEPREL")
    
    df <- df[!grepl("-", df$ID), ]
    df$ID <- as.numeric(df$ID)
    df$HEAD <- as.numeric(df$HEAD)
    
    return(df)
  }
  
  calculate_q <- function(df) {
    degrees <- sapply(df$ID, function(id) {
      num_children <- sum(df$HEAD == id)
      has_parent <- ifelse(df$HEAD[df$ID == id] != 0, 1, 0)
      return(num_children + has_parent)
    })
    
    q <- sum(degrees %% 2)
    return(q)
  }
  
  k2_table <- data.frame()
  d_table <- data.frame()
  
  for (pud_file in pud_files) {
    parent_folder <- basename(dirname(pud_file))
    lang_name <- sub("^UD_(.*)-PUD$", "\\1", parent_folder)
    
    metric_file <- file.path(data_dir, paste0(lang_name, "_dependency_tree_metrics.txt"))
    
    if (!file.exists(metric_file)) next
    
    metrics <- read.table(metric_file, header = FALSE, sep = "\t", 
                          colClasses = c("integer", "numeric", "numeric"))
    colnames(metrics) <- c("n", "k2", "d_mean") 
    
    lines <- readLines(pud_file, encoding = "UTF-8")
    sentence_blocks <- split(lines, cumsum(grepl("^# sent_id", lines)))
    sentences <- lapply(sentence_blocks, parse_sentence)
    sentences <- sentences[!sapply(sentences, is.null)]
    
    q_values <- sapply(sentences, calculate_q)
    metrics$q <- q_values[1:nrow(metrics)]
    
    metrics <- metrics %>%
      mutate(
        k2_lower_bound = 4 - 6/n,
        k2_upper_bound = n - 1,
        d_lower_bound = (1/4) * (n * k2 + q),
        d_upper_bound = (1/4) * (3 * (n - 1)^2 + 1 - (n %% 2)),
        d_lower_bound_normalized = d_lower_bound / (n - 1),
        d_upper_bound_normalized = d_upper_bound / (n - 1),
        
        k2 = round(k2, digits = 3),
        k2_lower_bound = round(k2_lower_bound, digits = 3),
        k2_upper_bound = round(k2_upper_bound, digits = 3),
        d_lower_bound_normalized = round(d_lower_bound_normalized, digits = 3),
        d_upper_bound_normalized = round(d_upper_bound_normalized, digits = 3),
        
        # Check validity
        k2_valid = (k2 >= k2_lower_bound) & (k2 <= k2_upper_bound),
        d_valid = (d_mean * (n-1) >= d_lower_bound - 1e-6) & (d_mean * (n-1) <= d_upper_bound + 1e-6)
      )
    
    k2_lang_table <- metrics %>%
      summarise(
        Language = lang_name,
        k2_mean = mean(k2, na.rm = TRUE),
        k2_lower_mean = mean(k2_lower_bound, na.rm = TRUE),
        k2_upper_mean = mean(k2_upper_bound, na.rm = TRUE),
        k2_valid_count = sum(k2_valid, na.rm = TRUE),
        k2_total = n()
      ) %>%
      mutate(
        k2_valid_status = ifelse(k2_valid_count == k2_total, "Yes", "No")
      )
    
    k2_table <- rbind(k2_table, k2_lang_table)
    
    d_lang_table <- metrics %>%
      summarise(
        Language = lang_name,
        d_mean_avg = mean(d_mean, na.rm = TRUE),
        d_lower_mean = mean(d_lower_bound_normalized, na.rm = TRUE),
        d_upper_mean = mean(d_upper_bound_normalized, na.rm = TRUE),
        d_valid_count = sum(d_valid, na.rm = TRUE),
        d_total = n()
      ) %>%
      mutate(
        d_valid_status = ifelse(d_valid_count == d_total, "Yes", "No")
      )
    
    d_table <- rbind(d_table, d_lang_table)
  }
  
  k2_table_formatted <- k2_table %>%
    mutate(
      k2_mean = round(k2_mean, 3),
      k2_lower_mean = round(k2_lower_mean, 3),
      k2_upper_mean = round(k2_upper_mean, 3)
    ) %>%
    select(Language, k2_mean, k2_lower_mean, k2_upper_mean, k2_valid_status)
  
  colnames(k2_table_formatted) <- c("Language", "k²", "Lower Bound", "Upper Bound", "Valid")
  
  cat("\n========================================\n")
  cat("k² VALIDATION TABLE\n")
  cat("========================================\n\n")
  print(k2_table_formatted, row.names = FALSE)
  cat("\n")
  
  d_table_formatted <- d_table %>%
    mutate(
      d_mean_avg = round(d_mean_avg, 3),
      d_lower_mean = round(d_lower_mean, 3),
      d_upper_mean = round(d_upper_mean, 3)
    ) %>%
    select(Language, d_mean_avg, d_lower_mean, d_upper_mean, d_valid_status)
  
  colnames(d_table_formatted) <- c("Language", "<d>", "Lower Bound", "Upper Bound", "Valid")
  
  cat("\n========================================\n")
  cat("<d> VALIDATION TABLE\n")
  cat("========================================\n\n")
  print(d_table_formatted, row.names = FALSE)
  cat("\n")
  
  # Return both tables invisibly for further use if needed
  invisible(list(k2_table = k2_table_formatted, d_table = d_table_formatted))
}



##############################################################################################################

library(dplyr)
compute.language.statistics <- function() {
  
  data_dir <- "./data"
  metric_files <- list.files(path = data_dir,
                             pattern = "_dependency_tree_metrics\\.txt$",
                             full.names = TRUE)
  
  results <- data.frame(
    Language = character(),
    N = integer(),
    mu_n = numeric(),
    sigma_n = numeric(),
    mu_d = numeric(),
    sigma_d = numeric(),
    stringsAsFactors = FALSE
  )
  
  language_tables <- list()
  
  for (metric_file in metric_files) {
    file_name <- basename(metric_file)
    lang_name <- sub("_dependency_tree_metrics\\.txt$", "", file_name)
    
    data <- read.table(metric_file, header = FALSE, sep = "\t",
                       colClasses = c("integer", "numeric", "numeric"))
    colnames(data) <- c("vertices", "degree_2nd_moment", "mean_length")
    
    # Order by number of vertices
    data <- data[order(data$vertices), ]
    
    # Store in the list
    language_tables[[lang_name]] <- data
    
    # Compute statistics (only for n and d_mean)
    N <- nrow(data)
    mu_n <- mean(data$vertices)
    sigma_n <- sd(data$vertices)
    mu_d <- mean(data$mean_length)
    sigma_d <- sd(data$mean_length)
    
    results <- rbind(results, data.frame(
      Language = lang_name,
      N = N,
      mu_n = mu_n,
      sigma_n = sigma_n,
      mu_d = mu_d,
      sigma_d = sigma_d
    ))
  }
  
  cat("\n========================================\n")
  cat("LANGUAGE STATISTICS SUMMARY\n")
  cat("========================================\n\n")
  
  results_display <- results %>%
    mutate(
      mu_n = round(mu_n, 2),
      sigma_n = round(sigma_n, 2),
      mu_d = round(mu_d, 2),
      sigma_d = round(sigma_d, 2)
    )
  
  print(results_display, row.names = FALSE)
  
  # Return both summary statistics and individual tables
  return(list(
    summary = results,
    tables = language_tables
  ))
}


#################################################################################################################################

logarithmic.binning.with.merge <- function(data, x_col, y_col, B = 1.2, 
                                           min_count = 5, alpha = NULL) {
  
  x <- data[[x_col]]
  y <- data[[y_col]]
  
  a0 <- min(x) * 0.95
  x_max <- max(x)
  num_bins <- ceiling(log(x_max / a0) / log(B))
  bin_edges <- a0 * (B ^ (0:num_bins))
  bin_indices <- findInterval(x, bin_edges, rightmost.closed = TRUE)
  
  binned_full <- data.frame(
    x = x,
    y = y,
    bin = bin_indices
  ) %>%
    filter(bin > 0, bin <= num_bins)
  
  bin_counts <- binned_full %>%
    group_by(bin) %>%
    summarise(count = n(), .groups = 'drop')
  
  bin_mapping <- data.frame(
    original_bin = 1:num_bins,
    merged_bin = 1:num_bins
  )
  
  for (i in 1:num_bins) {
    bin_count <- bin_counts$count[bin_counts$bin == i]
    
    if (length(bin_count) == 0 || bin_count < min_count) {
      merged <- FALSE
      for (j in (i+1):min(num_bins, i+10)) {
        next_count <- bin_counts$count[bin_counts$bin == j]
        if (length(next_count) > 0 && next_count >= min_count) {
          bin_mapping$merged_bin[i] <- j
          merged <- TRUE
          break
        }
      }
      if (!merged) {
        bin_mapping$merged_bin[i] <- i
      }
    }
  }

  first_valid_bin <- min(bin_counts$bin[bin_counts$count >= min_count])
  for (i in 1:num_bins) {
    if (i < first_valid_bin) {
      bin_count <- bin_counts$count[bin_counts$bin == i]
      if (length(bin_count) > 0) {
        bin_mapping$merged_bin[i] <- first_valid_bin
      }
    }
  }
  
  last_valid_bin <- max(bin_counts$bin[bin_counts$count >= min_count])
  for (i in num_bins:1) {
    if (i > last_valid_bin) {
      bin_count <- bin_counts$count[bin_counts$bin == i]
      if (length(bin_count) > 0) {
        bin_mapping$merged_bin[i] <- last_valid_bin
      }
    }
  }
  
  binned_full <- binned_full %>%
    left_join(bin_mapping %>% select(original_bin, merged_bin), 
              by = c("bin" = "original_bin")) %>%
    mutate(bin = merged_bin) %>%
    select(-merged_bin)
  
  binned_stats <- binned_full %>%
    group_by(bin) %>%
    summarise(
      count = n(),
      x_min = min(x),
      x_max = max(x),
      x_lower = a0 * (B ^ (min(bin) - 1)),
      x_upper = a0 * (B ^ max(bin)),
      x_geom = sqrt(x_lower * x_upper),
      y_mean = mean(y, na.rm = TRUE),
      y_sd = sd(y, na.rm = TRUE),
      y_se = sd(y, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  if (!is.null(alpha)) {
    binned_stats <- binned_stats %>%
      mutate(
        B_eff = x_upper / x_lower,
        correction_factor = ifelse(
          abs(alpha - 1) > 0.01,
          ((1 - B_eff^(2*(alpha-1))) / ((B_eff-1) * (B_eff^(alpha-1) - 1)))^(1/(alpha-1)),
          1
        ),
        x_star = x_geom * correction_factor
      )
  } else {
    binned_stats <- binned_stats %>%
      mutate(x_star = x_geom)
  }
  
  binned_full <- binned_full %>%
    left_join(binned_stats %>% select(bin, x_star), by = "bin")
  
  num_merged <- num_bins - nrow(binned_stats)
  
  return(list(
    binned_full = binned_full,
    binned_stats = binned_stats,
    bin_edges = bin_edges,
    B = B,
    a0 = a0,
    num_bins = num_bins,
    num_merged_bins = nrow(binned_stats),
    num_bins_merged = num_merged
  ))
}


#####################################################################################################

find_initial_parameters_binned <- function(lang_data, B = 1.2, min_count = 5, 
                                           alpha = NULL, use_binned_means = TRUE) {
  
  plot_dir <- ifelse(use_binned_means,
                     "./plots_initial_fit_binned_means",
                     "./plots_initial_fit_binned_raw")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  
  initial_params <- data.frame(
    Language  = character(),
    a_initial = numeric(),
    b_initial = numeric(),
    stringsAsFactors = FALSE
  )
  
  cat("\n========================================\n")
  if (use_binned_means) {
    cat("INITIAL PARAMETER ESTIMATION (BINNED MEANS)\n")
  } else {
    cat("INITIAL PARAMETER ESTIMATION (RAW DATA with x*)\n")
  }
  cat(sprintf("B = %.2f | Min count: %d\n", B, min_count))
  cat("========================================\n\n")
  
  plot_list <- list()
  
  for (lang_name in names(lang_data$tables)) {
    data <- lang_data$tables[[lang_name]]
    
    binned <- logarithmic.binning.with.merge(
      data      = data,
      x_col     = "vertices",
      y_col     = "mean_length",
      B         = B,
      min_count = min_count,
      alpha     = alpha
    )
    
    if (use_binned_means) {
      fit_data <- binned$binned_stats
      valid <- (fit_data$x_star > 0) & (fit_data$y_mean > 0)
      fit_data <- fit_data[valid, ]
      
      if (nrow(fit_data) < 2) next
      
      # Linear fit in log-log space
      fit <- lm(log(y_mean) ~ log(x_star), data = fit_data)
      
      coefs  <- coef(fit)
      a_init <- exp(coefs[1])
      b_init <- coefs[2]
      
      data_plot <- fit_data %>%
        mutate(
          log_x_star = log(x_star),
          log_y_mean = log(y_mean),
          fitted     = predict(fit)
        )
      
      p <- ggplot(data_plot, aes(x = log_x_star, y = log_y_mean)) +
        geom_point(aes(color = "Data"), size = 3) +
        geom_line(aes(y = fitted, color = "Linear fit"), linewidth = 1) +
        scale_color_manual(
          name   = NULL,
          values = c("Data" = "darkblue", "Linear fit" = "red")
        ) +
        labs(
          title    = paste(lang_name, "- Log-Log Linear Fit (Binned Means)"),
          subtitle = sprintf("B = %.2f | %d bins", B, nrow(fit_data)),
          x = "log(x*)",
          y = "log(mean dependency length)"
        ) +
        theme_minimal() +
        theme(
          plot.title   = element_text(size = 10, face = "bold"),
          plot.subtitle= element_text(size = 8),
          axis.title   = element_text(size = 8),
          axis.text    = element_text(size = 7),
          legend.position = "bottom",
          legend.text  = element_text(size = 6)
        )
      
    } else {
      fit_data <- binned$binned_full
      valid <- (fit_data$x_star > 0) & (fit_data$y > 0)
      fit_data <- fit_data[valid, ]
      
      if (nrow(fit_data) < 2) next
      
      # Linear fit in log-log space
      fit <- lm(log(y) ~ log(x_star), data = fit_data)
      
      coefs  <- coef(fit)
      a_init <- exp(coefs[1])
      b_init <- coefs[2]
      
      data_plot <- fit_data %>%
        mutate(
          log_x_star = log(x_star),
          log_y      = log(y),
          fitted     = predict(fit)
        )
      
      p <- ggplot(data_plot, aes(x = log_x_star, y = log_y)) +
        geom_point(aes(color = "Data"), size = 1.5, alpha = 0.3) +
        geom_line(aes(y = fitted, color = "Linear fit"), linewidth = 1) +
        scale_color_manual(
          name   = NULL,
          values = c("Data" = "darkblue", "Linear fit" = "red")
        ) +
        labs(
          title    = paste(lang_name, "- Log-Log Linear Fit (Raw Data)"),
          subtitle = sprintf("B = %.2f | %d points", B, nrow(fit_data)),
          x = "log(x*)",
          y = "log(mean dependency length)"
        ) +
        theme_minimal() +
        theme(
          plot.title   = element_text(size = 10, face = "bold"),
          plot.subtitle= element_text(size = 8),
          axis.title   = element_text(size = 8),
          axis.text    = element_text(size = 7),
          legend.position = "bottom",
          legend.text  = element_text(size = 6)
        )
    }
    
    initial_params <- rbind(initial_params, data.frame(
      Language  = lang_name,
      a_initial = a_init,
      b_initial = b_init,
      row.names = NULL
    ))
    
    cat(sprintf("%-15s  a: %.3f   b: %.3f\n", lang_name, a_init, b_init))
    
    plot_list[[lang_name]] <- p
  }
  
  cat("\n========================================\n\n")
  
  # Combine all subplots in a 7x3 grid
  title_text <- if (use_binned_means) {
    "Initial log-log fits (binned means)"
  } else {
    "Initial log-log fits (raw data)"
  }
  
  combined <- wrap_plots(plot_list, ncol = 3, nrow = 7, guides = "collect") +
    plot_annotation(
      title = title_text,
      subtitle = sprintf("B = %.2f | Min count = %d", B, min_count),
      theme = theme(
        plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)
      )
    )
  
  output_file <- if (use_binned_means) {
    file.path(plot_dir, "loglog_initial_fits_binned_means_grid.png")
  } else {
    file.path(plot_dir, "loglog_initial_fits_raw_grid.png")
  }
  
  ggsave(output_file, plot = combined, width = 18, height = 24, dpi = 150)
  cat("Grid plot saved to:", output_file, "\n")
  
  return(initial_params)
}


#####################################################################################################

library(dplyr)
library(tibble)

fit_all_models_binned <- function(lang_data, initial_params, B = 1.2, min_count = 5, 
                                  alpha = NULL, use_binned_means = TRUE) {
  
  plot_dir <- ifelse(use_binned_means, 
                     "./plots_nls_fits_binned_means", 
                     "./plots_nls_fits_binned_raw")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  
  results_df <- data.frame(
    Language = character(),
    Model = character(),
    Description = character(),
    Param_names = character(),
    Param_values = character(),
    Param_SEs = character(),
    RSS = numeric(),
    Residual_SE = numeric(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  cov_results <- tibble(
    Language   = character(),
    Model      = character(),
    Description= character(),
    Covariance = list()
  )
  
  cat("\n========================================\n")
  if (use_binned_means) {
    cat("FITTING MODELS ON BINNED MEANS\n")
  } else {
    cat("FITTING MODELS ON RAW DATA (with binned x*)\n")
  }
  cat(sprintf("B = %.2f | Min count: %d\n", B, min_count))
  cat("========================================\n\n")
  
  for (i in 1:nrow(initial_params)) {
    lang_name <- initial_params$Language[i]
    data <- lang_data$tables[[lang_name]]
    
    # Apply logarithmic binning 
    binned <- logarithmic.binning.with.merge(
      data = data,
      x_col = "vertices",
      y_col = "mean_length",
      B = B,
      min_count = min_count,
      alpha = alpha
    )
    
    if (use_binned_means) {
      fit_data <- binned$binned_stats
      valid <- (fit_data$x_star > 0) & (fit_data$y_mean > 0)
      fit_data <- fit_data[valid, ]
      x <- fit_data$x_star
      y <- fit_data$y_mean
    } else {
      fit_data <- binned$binned_full
      valid <- (fit_data$x_star > 0) & (fit_data$y > 0)
      fit_data <- fit_data[valid, ]
      x <- fit_data$x_star
      y <- fit_data$y
    }
    
    if (length(x) < 2) next
    
    n_obs <- length(x)
    
    models <- list()
    AICs <- numeric()
    
    ## Model 0: null model f(x*) = (x*+1)/3 (nessun parametro → niente vcov)
    model0_pred <- (x + 1) / 3
    RSS_0 <- sum((y - model0_pred)^2)
    p_0 <- 0
    AIC_0 <- n_obs * log(2 * pi) + n_obs * log(RSS_0 / n_obs) + n_obs + 2 * (p_0 + 1)
    s_0 <- sqrt(RSS_0 / max(1, n_obs - p_0))
    AICs <- c(AICs, AIC_0)
    models$model0 <- list(
      description = "(x*+1)/3", 
      param_names = NA, 
      param_values = NA, 
      param_ses = NA, 
      RSS = RSS_0, 
      s = s_0, 
      AIC = AIC_0, 
      fitted = model0_pred,
      vcov = NULL
    )
    
    ## Model 1: f(x*) = (x*/2)^b
    model1_nls <- nls(
      y ~ (x / 2)^b,
      start = list(b = initial_params$b_initial[i]),
      control = nls.control(warnOnly = TRUE, maxiter = 100)
    )
    summary_fit1 <- summary(model1_nls)
    se_params1  <- summary_fit1$parameters[, "Std. Error"]
    vcov1        <- vcov(model1_nls)
    
    b_1   <- coef(model1_nls)["b"]
    RSS_1 <- sum(resid(model1_nls)^2)
    p_1   <- 1
    AIC_1 <- n_obs * log(2 * pi) + n_obs * log(RSS_1 / n_obs) + n_obs + 2 * (p_1 + 1)
    s_1   <- sqrt(RSS_1 / (n_obs - p_1))
    AICs  <- c(AICs, AIC_1)
    models$model1 <- list(
      description = "(x*/2)^b", 
      param_names = "b", 
      param_values = as.numeric(b_1), 
      param_ses = as.numeric(se_params1), 
      RSS = RSS_1, 
      s = s_1, 
      AIC = AIC_1, 
      fitted = fitted(model1_nls),
      vcov = vcov1
    )
    
    ## Model 2: f(x*) = a * x*^b
    model2_nls <- nls(
      y ~ a * x ^ b,
      start = list(a = initial_params$a_initial[i],
                   b = initial_params$b_initial[i]),
      control = nls.control(warnOnly = TRUE, maxiter = 100)
    )
    summary_fit2 <- summary(model2_nls)
    se_params2   <- summary_fit2$parameters[, "Std. Error"]
    vcov2        <- vcov(model2_nls)
    
    a_2   <- coef(model2_nls)["a"]
    b_2   <- coef(model2_nls)["b"]
    RSS_2 <- sum(resid(model2_nls)^2)
    p_2   <- 2
    AIC_2 <- n_obs * log(2 * pi) + n_obs * log(RSS_2 / n_obs) + n_obs + 2 * (p_2 + 1)
    s_2   <- sqrt(RSS_2 / (n_obs - p_2))
    AICs  <- c(AICs, AIC_2)
    models$model2 <- list(
      description = "a * x*^b", 
      param_names = paste(c("a", "b"), collapse = ";"), 
      param_values = paste(signif(c(a_2, b_2), 4), collapse = ";"), 
      param_ses = paste(signif(se_params2, 4), collapse = ";"), 
      RSS = RSS_2, 
      s = s_2, 
      AIC = AIC_2, 
      fitted = fitted(model2_nls),
      vcov = vcov2
    )
    
    ## Model 3: f(x*) = a * exp(c * x*)
    model3_nls <- nls(
      y ~ a * exp(c * x),
      start = list(a = initial_params$a_initial[i], c = 0.01),
      control = nls.control(warnOnly = TRUE, maxiter = 100)
    )
    summary_fit3 <- summary(model3_nls)
    se_params3   <- summary_fit3$parameters[, "Std. Error"]
    vcov3        <- vcov(model3_nls)
    
    a_3   <- coef(model3_nls)["a"]
    c_3   <- coef(model3_nls)["c"]
    RSS_3 <- sum(resid(model3_nls)^2)
    p_3   <- 2
    AIC_3 <- n_obs * log(2 * pi) + n_obs * log(RSS_3 / n_obs) + n_obs + 2 * (p_3 + 1)
    s_3   <- sqrt(RSS_3 / (n_obs - p_3))
    AICs  <- c(AICs, AIC_3)
    models$model3 <- list(
      description = "a * exp(c*x*)", 
      param_names = paste(c("a", "c"), collapse = ";"), 
      param_values = paste(signif(c(a_3, c_3), 4), collapse = ";"), 
      param_ses = paste(signif(se_params3, 4), collapse = ";"), 
      RSS = RSS_3, 
      s = s_3, 
      AIC = AIC_3, 
      fitted = fitted(model3_nls),
      vcov = vcov3
    )
    
    ## Model 4: f(x*) = a * log(x*)
    model4_nls <- nls(
      y ~ a * log(x),
      start = list(a = initial_params$a_initial[i]),
      control = nls.control(warnOnly = TRUE, maxiter = 100)
    )
    summary_fit4 <- summary(model4_nls)
    se_params4   <- summary_fit4$parameters[, "Std. Error"]
    vcov4        <- vcov(model4_nls)
    
    a_4   <- coef(model4_nls)["a"]
    RSS_4 <- sum(resid(model4_nls)^2)
    p_4   <- 1
    AIC_4 <- n_obs * log(2 * pi) + n_obs * log(RSS_4 / n_obs) + n_obs + 2 * (p_4 + 1)
    s_4   <- sqrt(RSS_4 / (n_obs - p_4))
    AICs  <- c(AICs, AIC_4)
    models$model4 <- list(
      description = "a * log(x*)", 
      param_names = "a", 
      param_values = as.numeric(a_4), 
      param_ses = as.numeric(se_params4), 
      RSS = RSS_4, 
      s = s_4, 
      AIC = AIC_4, 
      fitted = fitted(model4_nls),
      vcov = vcov4
    )
    
    ## Model 5: f(x*) = a * x*^b * exp(c*x*)
    model5_nls <- nlsLM(
      y ~ a * x^b * exp(c * x),
      start = list(a = initial_params$a_initial[i],
                   b = initial_params$b_initial[i],
                   c = 0),
      control = nls.lm.control(maxiter = 500)
    )
    summary_fit5 <- summary(model5_nls)
    se_params5   <- summary_fit5$parameters[, "Std. Error"]
    vcov5        <- vcov(model5_nls)
    
    a_5   <- coef(model5_nls)["a"]
    b_5   <- coef(model5_nls)["b"]
    c_5   <- coef(model5_nls)["c"]
    RSS_5 <- sum(resid(model5_nls)^2)
    p_5   <- 3
    AIC_5 <- n_obs * log(2 * pi) + n_obs * log(RSS_5 / n_obs) + n_obs + 2 * (p_5 + 1)
    s_5   <- sqrt(RSS_5 / (n_obs - p_5))
    AICs  <- c(AICs, AIC_5)
    models$model5 <- list(
      description = "a * x*^b * exp(c*x*)", 
      param_names = paste(c("a", "b", "c"), collapse = ";"), 
      param_values = paste(signif(c(a_5, b_5, c_5), 4), collapse = ";"), 
      param_ses = paste(signif(se_params5, 4), collapse = ";"), 
      RSS = RSS_5, 
      s = s_5, 
      AIC = AIC_5, 
      fitted = fitted(model5_nls),
      vcov = vcov5
    )
    
    ## Delta AIC
    min_AIC <- min(AICs)
    delta_AICs <- AICs - min_AIC
    
    for (j in seq_along(models)) {
      model_name <- names(models)[j]
      m <- models[[model_name]]
      
      results_df <- rbind(results_df, data.frame(
        Language = lang_name,
        Model = model_name,
        Description = m$description,
        Param_names = ifelse(is.na(m$param_names), NA, m$param_names),
        Param_values = ifelse(is.na(m$param_values), NA, as.character(m$param_values)),
        Param_SEs = ifelse(is.na(m$param_ses), NA, as.character(m$param_ses)),
        RSS = m$RSS,
        Residual_SE = m$s,
        AIC = m$AIC,
        Delta_AIC = delta_AICs[j],
        stringsAsFactors = FALSE
      ))
      
      if (!is.null(m$vcov)) {
        cov_results <- add_row(
          cov_results,
          Language   = lang_name,
          Model      = model_name,
          Description= m$description,
          Covariance = list(m$vcov)
        )
      }
    }
    
    plot_data <- data.frame(
      x_star = x,
      observed = y,
      model0 = models$model0$fitted,
      model1 = models$model1$fitted,
      model2 = models$model2$fitted,
      model3 = models$model3$fitted,
      model4 = models$model4$fitted,
      model5 = models$model5$fitted
    )
    alpha_points <- if (use_binned_means) 1 else 0.3
    size_points  <- if (use_binned_means) 3 else 1
    
    y_max <- max(y) * 1.1
    y_min <- max(0, min(y) * 0.9)
    
    p <- ggplot(plot_data, aes(x = x_star)) +
      geom_point(aes(y = observed), size = size_points,
                 color = "black", alpha = alpha_points) +
      geom_line(aes(y = model0, color = "Model 0: (x*+1)/3"), linewidth = 0.8) +
      geom_line(aes(y = model1, color = "Model 1: (x*/2)^b"), linewidth = 0.8) +
      geom_line(aes(y = model2, color = "Model 2: a*x*^b"), linewidth = 0.8) +
      geom_line(aes(y = model3, color = "Model 3: a*exp(c*x*)"), linewidth = 0.8) +
      geom_line(aes(y = model4, color = "Model 4: a*log(x*)"), linewidth = 0.8) +
      geom_line(aes(y = model5, color = "Model 5: a*x*^b*exp(c*x*)"), linewidth = 0.8) +
      scale_x_log10() +
      scale_y_log10() +
      coord_cartesian(ylim = c(y_min, y_max)) +
      scale_color_manual(
        name = "Models",
        values = c(
          "Model 0: (x*+1)/3"              = "red",
          "Model 1: (x*/2)^b"              = "blue",
          "Model 2: a*x*^b"                = "green",
          "Model 3: a*exp(c*x*)"           = "orange",
          "Model 4: a*log(x*)"             = "purple",
          "Model 5: a*x*^b*exp(c*x*)"      = "brown"
        )
      ) +
      labs(
        title = paste(
          lang_name,
          ifelse(use_binned_means,
                 "- Nonlinear model fits (binned means)",
                 "- Nonlinear model fits (raw data with x*)")
        ),
        subtitle = sprintf("B = %.2f | N = %d", B, n_obs),
        x = expression(n^"*"~"(corrected bin position)"),
        y = "Mean dependency length"
      ) +
      theme_minimal() +
      theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title    = element_text(size = 10),
        axis.text     = element_text(size = 9),
        legend.position = "bottom",
        legend.text   = element_text(size = 8)
      )
    
    
    output_file <- paste0(plot_dir, "/", lang_name, "_nls_comparison.png")
    ggsave(output_file, plot = p, width = 10, height = 7, dpi = 120)
  }
  
  return(list(
    results_df  = results_df,
    cov_results = cov_results
  ))
}

###################################################################################################



library(dplyr)
library(ggplot2)
library(patchwork)

plot_covariance_grids_plus <- function(cov_results,
                                       out_dir = "./plots_corr_matrices_d") {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  models_plus <- c("model1", "model2", "model3", "model4", "model5")
  
  plots_by_model <- list()
  
  for (m in models_plus) {
    df_m <- cov_results %>%
      filter(Model == m)
    
    if (nrow(df_m) == 0) next
    
    lang_plots <- vector("list", length = nrow(df_m))
    
    for (i in seq_len(nrow(df_m))) {
      lang <- df_m$Language[i]
      desc <- df_m$Description[i]
      cov_mat <- df_m$Covariance[[i]]
      
      cor_mat <- cov2cor(cov_mat)
      
      cor_df <- as.data.frame(as.table(cor_mat))
      colnames(cor_df) <- c("Param1", "Param2", "Corr")
      
      cor_df$Param1 <- factor(cor_df$Param1, levels = colnames(cor_mat))
      cor_df$Param2 <- factor(cor_df$Param2, levels = rownames(cor_mat))
      
      p <- ggplot(cor_df, aes(x = Param1, y = Param2, fill = Corr)) +
        geom_tile(color = "grey80") +
        geom_text(aes(label = sprintf("%.2f", Corr)),
                  size = 2.5, color = "black") +
        scale_fill_gradient2(
          name   = "corr",
          low    = "#4575b4",
          mid    = "white",
          high   = "#d73027",
          limits = c(-1, 1),
          breaks = c(-1, -0.5, 0, 0.5, 1)
        ) +
        labs(
          title    = lang,
          subtitle = desc,
          x = NULL,
          y = NULL
        ) +
        theme_minimal() +
        theme(
          plot.title    = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6),
          axis.text.x   = element_text(size = 6, angle = 45, hjust = 1),
          axis.text.y   = element_text(size = 6)
        )
      
      lang_plots[[i]] <- p
    }
    
    combined <- wrap_plots(lang_plots, ncol = 3, nrow = 7, guides = "collect") +
      plot_annotation(
        title = paste("Parameter correlation matrix -", toupper(m))
      ) &
      theme(
        legend.position = "right",
        legend.title    = element_text(size = 7),
        legend.text     = element_text(size = 6),
        plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    
    file_name <- file.path(out_dir, paste0("corr_grid_", m, "_plus.png"))
    ggsave(
      filename = file_name,
      plot     = combined,
      width    = 18,
      height   = 24,
      dpi      = 150
    )
    
    message("Saved correlation grid: ", file_name)
    plots_by_model[[m]] <- combined
  }
  
  return(plots_by_model)
}
###################################################################################################

print_best_models <- function(results_df) {
  
  cat("\n========================================\n")
  cat("BEST MODEL BY LANGUAGE (minimum AIC)\n")
  cat("========================================\n\n")

  best_models <- results_df %>%
    group_by(Language) %>%
    slice_min(AIC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Language, Model, Description, AIC, Delta_AIC) %>%
    arrange(Language)
  

  for (i in 1:nrow(best_models)) {
    cat(sprintf("%-15s: %s (%s)\n", 
                best_models$Language[i], 
                best_models$Model[i],
                best_models$Description[i]))
  }
  
  cat("\n========================================\n\n")
  
  return(best_models)
}

#########################################################################################################################################

fit_all_models_binned_with_d <- function(lang_data, initial_params, B = 1.2, min_count = 5, 
                                         alpha = NULL, use_binned_means = TRUE, maxiter = 500) {
  
  plot_dir <- ifelse(use_binned_means,
                     "./plots_nls_fits_binned_means_d",
                     "./plots_nls_fits_binned_raw_d")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  
  results_df <- data.frame(
    Language = character(),
    Model = character(),
    Description = character(),
    Param_names = character(),
    Param_values = character(),
    Param_SEs = character(),
    RSS = numeric(),
    Residual_SE = numeric(),
    AIC = numeric(),
    Delta_AIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  cov_results <- tibble(
    Language   = character(),
    Model      = character(),
    Description= character(),
    Covariance = list()
  )
  
  cat("\n========================================\n")
  cat("FITTING MODELS WITH INTERCEPT d (using nlsLM)\n")
  if (use_binned_means) {
    cat("Data: Binned means\n")
  } else {
    cat("Data: Raw data with x*\n")
  }
  cat("Initial values: a, b from linear interpolation | c, d = 0\n")
  cat(sprintf("B = %.2f | Min count: %d\n", B, min_count))
  cat("========================================\n\n")
  
  for (i in 1:nrow(initial_params)) {
    lang_name <- initial_params$Language[i]
    data <- lang_data$tables[[lang_name]]
    
    binned <- logarithmic.binning.with.merge(
      data = data,
      x_col = "vertices",
      y_col = "mean_length",
      B = B,
      min_count = min_count,
      alpha = alpha
    )
    
    if (use_binned_means) {
      fit_data <- binned$binned_stats
      valid <- (fit_data$x_star > 0) & (fit_data$y_mean > 0)
      fit_data <- fit_data[valid, ]
      x <- fit_data$x_star
      y <- fit_data$y_mean
    } else {
      fit_data <- binned$binned_full
      valid <- (fit_data$x_star > 0) & (fit_data$y > 0)
      fit_data <- fit_data[valid, ]
      x <- fit_data$x_star
      y <- fit_data$y
    }
    
    if (length(x) < 3) next
    
    n_obs <- length(x)
    a_init <- initial_params$a_initial[i]
    b_init <- initial_params$b_initial[i]
    d_init <- 0
    c_init <- 0
    
    models <- list()
    AICs <- numeric()
    
    ## Model 0 (null) – come prima, ma senza vcov (niente parametri)
    model0_pred <- (x + 1) / 3
    RSS_0 <- sum((y - model0_pred)^2)
    p_0 <- 0
    AIC_0 <- n_obs * log(2 * pi) + n_obs * log(RSS_0 / n_obs) + n_obs + 2 * (p_0 + 1)
    s_0 <- sqrt(RSS_0 / max(1, n_obs - p_0))
    AICs <- c(AICs, AIC_0)
    models$model0 <- list(
      description = "(x*+1)/3",
      param_names = NA,
      param_values = NA,
      param_ses = NA,
      RSS = RSS_0,
      s = s_0,
      AIC = AIC_0,
      fitted = model0_pred,
      vcov  = NULL
    )
    
    ## Model 1: (x*/2)^b + d
    model1_nls <- nlsLM(y ~ (x / 2)^b + d,
                        start = list(b = b_init, d = d_init),
                        control = nls.lm.control(maxiter = maxiter))
    sum1 <- summary(model1_nls)
    se_params1 <- sum1$parameters[, "Std. Error"]
    vcov1 <- vcov(model1_nls)  # <- matrice di var-cov dei parametri
    
    b_1 <- coef(model1_nls)["b"]
    d_1 <- coef(model1_nls)["d"]
    RSS_1 <- sum(resid(model1_nls)^2)
    p_1 <- 2
    AIC_1 <- n_obs * log(2 * pi) + n_obs * log(RSS_1 / n_obs) + n_obs + 2 * (p_1 + 1)
    s_1 <- sqrt(RSS_1 / (n_obs - p_1))
    AICs <- c(AICs, AIC_1)
    models$model1 <- list(
      description = "(x*/2)^b + d",
      param_names = paste(c("b", "d"), collapse = ";"),
      param_values = paste(signif(c(b_1, d_1), 4), collapse = ";"),
      param_ses = paste(signif(se_params1, 4), collapse = ";"),
      RSS = RSS_1,
      s = s_1,
      AIC = AIC_1,
      fitted = fitted(model1_nls),
      vcov  = vcov1
    )
    
    ## Model 2: a * x*^b + d
    model2_nls <- nlsLM(y ~ a * x^b + d,
                        start = list(a = a_init, b = b_init, d = d_init),
                        control = nls.lm.control(maxiter = maxiter))
    sum2 <- summary(model2_nls)
    se_params2 <- sum2$parameters[, "Std. Error"]
    vcov2 <- vcov(model2_nls)
    
    a_2 <- coef(model2_nls)["a"]
    b_2 <- coef(model2_nls)["b"]
    d_2 <- coef(model2_nls)["d"]
    RSS_2 <- sum(resid(model2_nls)^2)
    p_2 <- 3
    AIC_2 <- n_obs * log(2 * pi) + n_obs * log(RSS_2 / n_obs) + n_obs + 2 * (p_2 + 1)
    s_2 <- sqrt(RSS_2 / (n_obs - p_2))
    AICs <- c(AICs, AIC_2)
    models$model2 <- list(
      description = "a * x*^b + d",
      param_names = paste(c("a", "b", "d"), collapse = ";"),
      param_values = paste(signif(c(a_2, b_2, d_2), 4), collapse = ";"),
      param_ses = paste(signif(se_params2, 4), collapse = ";"),
      RSS = RSS_2,
      s = s_2,
      AIC = AIC_2,
      fitted = fitted(model2_nls),
      vcov  = vcov2
    )
    
    ## Model 3: a * exp(c*x*) + d
    model3_nls <- nlsLM(y ~ a * exp(c * x) + d,
                        start = list(a = a_init, c = 0.01, d = d_init),
                        control = nls.lm.control(maxiter = maxiter))
    sum3 <- summary(model3_nls)
    se_params3 <- sum3$parameters[, "Std. Error"]
    vcov3 <- vcov(model3_nls)
    
    a_3 <- coef(model3_nls)["a"]
    c_3 <- coef(model3_nls)["c"]
    d_3 <- coef(model3_nls)["d"]
    RSS_3 <- sum(resid(model3_nls)^2)
    p_3 <- 3
    AIC_3 <- n_obs * log(2 * pi) + n_obs * log(RSS_3 / n_obs) + n_obs + 2 * (p_3 + 1)
    s_3 <- sqrt(RSS_3 / (n_obs - p_3))
    AICs <- c(AICs, AIC_3)
    models$model3 <- list(
      description = "a * exp(c*x*) + d",
      param_names = paste(c("a", "c", "d"), collapse = ";"),
      param_values = paste(signif(c(a_3, c_3, d_3), 4), collapse = ";"),
      param_ses = paste(signif(se_params3, 4), collapse = ";"),
      RSS = RSS_3,
      s = s_3,
      AIC = AIC_3,
      fitted = fitted(model3_nls),
      vcov  = vcov3
    )
    
    ## Model 4: a * log(x*) + d
    model4_nls <- nlsLM(y ~ a * log(x) + d,
                        start = list(a = a_init, d = d_init),
                        control = nls.lm.control(maxiter = maxiter))
    sum4 <- summary(model4_nls)
    se_params4 <- sum4$parameters[, "Std. Error"]
    vcov4 <- vcov(model4_nls)
    
    a_4 <- coef(model4_nls)["a"]
    d_4 <- coef(model4_nls)["d"]
    RSS_4 <- sum(resid(model4_nls)^2)
    p_4 <- 2
    AIC_4 <- n_obs * log(2 * pi) + n_obs * log(RSS_4 / n_obs) + n_obs + 2 * (p_4 + 1)
    s_4 <- sqrt(RSS_4 / (n_obs - p_4))
    AICs <- c(AICs, AIC_4)
    models$model4 <- list(
      description = "a * log(x*) + d",
      param_names = paste(c("a", "d"), collapse = ";"),
      param_values = paste(signif(c(a_4, d_4), 4), collapse = ";"),
      param_ses = paste(signif(se_params4, 4), collapse = ";"),
      RSS = RSS_4,
      s = s_4,
      AIC = AIC_4,
      fitted = fitted(model4_nls),
      vcov  = vcov4
    )
    
    ## Model 5: a * x*^b * exp(c*x*) + d
    model5_nls <- nlsLM(y ~ a * x^b * exp(c * x) + d,
                        start = list(a = a_init, b = b_init, c = c_init, d = d_init),
                        control = nls.lm.control(maxiter = maxiter))
    sum5 <- summary(model5_nls)
    se_params5 <- sum5$parameters[, "Std. Error"]
    vcov5 <- vcov(model5_nls)
    
    a_5 <- coef(model5_nls)["a"]
    b_5 <- coef(model5_nls)["b"]
    c_5 <- coef(model5_nls)["c"]
    d_5 <- coef(model5_nls)["d"]
    RSS_5 <- sum(resid(model5_nls)^2)
    p_5 <- 4
    AIC_5 <- n_obs * log(2 * pi) + n_obs * log(RSS_5 / n_obs) + n_obs + 2 * (p_5 + 1)
    s_5 <- sqrt(RSS_5 / (n_obs - p_5))
    AICs <- c(AICs, AIC_5)
    models$model5 <- list(
      description = "a * x*^b * exp(c*x*) + d",
      param_names = paste(c("a", "b", "c", "d"), collapse = ";"),
      param_values = paste(signif(c(a_5, b_5, c_5, d_5), 4), collapse = ";"),
      param_ses = paste(signif(se_params5, 4), collapse = ";"),
      RSS = RSS_5,
      s = s_5,
      AIC = AIC_5,
      fitted = fitted(model5_nls),
      vcov  = vcov5
    )
    
    # Delta AIC
    min_AIC <- min(AICs)
    delta_AICs <- AICs - min_AIC
    
    for (j in seq_along(models)) {
      model_name <- names(models)[j]
      m <- models[[model_name]]
      
      results_df <- rbind(results_df, data.frame(
        Language = lang_name,
        Model = model_name,
        Description = m$description,
        Param_names = ifelse(is.na(m$param_names), NA, m$param_names),
        Param_values = ifelse(is.na(m$param_values), NA, as.character(m$param_values)),
        Param_SEs = ifelse(is.na(m$param_ses), NA, as.character(m$param_ses)),
        RSS = m$RSS,
        Residual_SE = m$s,
        AIC = m$AIC,
        Delta_AIC = delta_AICs[j],
        stringsAsFactors = FALSE
      ))
      
      # aggiungi la vcov se esiste
      if (!is.null(m$vcov)) {
        cov_results <- add_row(
          cov_results,
          Language   = lang_name,
          Model      = model_name,
          Description= m$description,
          Covariance = list(m$vcov)
        )
      }
    }
    
    # Create plot
    plot_data <- data.frame(
      x_star = x,
      observed = y,
      model0 = models$model0$fitted,
      model1 = models$model1$fitted,
      model2 = models$model2$fitted,
      model3 = models$model3$fitted,
      model4 = models$model4$fitted,
      model5 = models$model5$fitted
    )
    
    y_max <- max(y) * 1.1
    y_min <- max(0, min(y) * 0.9)
    
    if (use_binned_means) {
      alpha_pts <- 1
      size_pts <- 3
    } else {
      alpha_pts <- 0.3
      size_pts <- 1
    }
    
    p <- ggplot(plot_data, aes(x = x_star)) +
      geom_point(aes(y = observed), alpha = alpha_pts, size = size_pts, color = "black") +
      geom_line(aes(y = model0, color = "Model 0: (x*+1)/3"), linewidth = 0.8) +
      geom_line(aes(y = model1, color = "Model 1: (x*/2)^b + d"), linewidth = 0.8) +
      geom_line(aes(y = model2, color = "Model 2: a*x*^b + d"), linewidth = 0.8) +
      geom_line(aes(y = model3, color = "Model 3: a*exp(c*x*) + d"), linewidth = 0.8) +
      geom_line(aes(y = model4, color = "Model 4: a*log(x*) + d"), linewidth = 0.8) +
      geom_line(aes(y = model5, color = "Model 5: a*x*^b*exp(c*x*) + d"), linewidth = 0.8) +
      scale_x_log10() +
      scale_y_log10() +
      coord_cartesian(ylim = c(y_min, y_max)) +
      scale_color_manual(
        name = "Models",
        values = c(
          "Model 0: (x*+1)/3" = "red",
          "Model 1: (x*/2)^b + d" = "blue",
          "Model 2: a*x*^b + d" = "green",
          "Model 3: a*exp(c*x*) + d" = "orange",
          "Model 4: a*log(x*) + d" = "purple",
          "Model 5: a*x*^b*exp(c*x*) + d" = "brown"
        )
      ) +
      labs(
        title = paste(lang_name, "- Model Fits with intercept d (nlsLM)"),
        subtitle = sprintf("B = %.2f | N = %d | c_init, d_init = 0", B, n_obs),
        x = "x* (corrected bin position)",
        y = "Mean dependency length"
      ) +
      labs(
        title = paste(
          lang_name,
          ifelse(use_binned_means,
                 "- Nonlinear model fits (binned means)",
                 "- Nonlinear model fits (raw data with x*)")
        ),
        subtitle = sprintf("B = %.2f | N = %d", B, n_obs),
        x = expression(n^"*"~"(corrected bin position)"),
        y = "Mean dependency length"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 8)
      )
    
    output_file <- paste0(plot_dir, "/", lang_name, "_nls_comparison_d.png")
    ggsave(output_file, plot = p, width = 10, height = 7, dpi = 120)
  }
  
  return(list(
    results_df = results_df,
    cov_results = cov_results
  ))
}


#########################################################################################################################################


extract_all_parameters <- function(results_no_d, results_with_d) {
  
  combined_results <- bind_rows(
    results_no_d %>% mutate(Model_Type = "no_d"),
    results_with_d %>% mutate(Model_Type = "with_d")
  )
  
  combined_results <- combined_results %>%
    mutate(
      Model_Label = case_when(
        Model_Type == "no_d" ~ gsub("model", "Model_", Model),
        Model_Type == "with_d" ~ paste0(gsub("model", "Model_", Model), "_plus")
      )
    )
  
  all_params <- combined_results %>%
    select(Language, Model_Label, Param_names, Param_values, Param_SEs) %>%
    arrange(Language, Model_Label)
  
  print(all_params)
  
  return(all_params)
}

################################################################################################################################
library(dplyr)
library(tidyr) 
create_aic_tables <- function(results_no_d, results_with_d) {
  
  combined_results <- bind_rows(
    results_no_d %>% mutate(Model_Type = "no_d"),
    results_with_d %>% mutate(Model_Type = "with_d")
  )
  
  combined_results <- combined_results %>%
    mutate(
      Model_Label = case_when(
        Model_Type == "no_d" ~ paste0("AIC_", gsub("model", "", Model)),
        Model_Type == "with_d" ~ paste0("AIC_", gsub("model", "", Model), "_plus")
      )
    )
  
  model_order <- c(
    "Language",
    "AIC_0", "AIC_1", "AIC_2", "AIC_3", "AIC_4", "AIC_5",
    "AIC_1_plus", "AIC_2_plus", "AIC_3_plus", "AIC_4_plus", "AIC_5_plus"
  )
  
  delta_aic <- combined_results %>%
    group_by(Language) %>%
    mutate(Delta_AIC = AIC - min(AIC)) %>%
    ungroup() %>%
    select(Language, Model_Label, Delta_AIC) %>%
    pivot_wider(names_from = Model_Label, values_from = Delta_AIC) %>%
    arrange(Language) %>%
    select(intersect(model_order, names(.)))
  
  
  absolute_aic <- combined_results %>%
    select(Language, Model_Label, AIC) %>%
    pivot_wider(names_from = Model_Label, values_from = AIC) %>%
    arrange(Language) %>%
    select(intersect(model_order, names(.)))
  
  print(delta_aic)
  print(absolute_aic)
  
  return(list(
    delta_aic = delta_aic,
    absolute_aic = absolute_aic
  ))
}

#########################################################################################################################################

create_residual_se_table <- function(results_no_d, results_with_d) {
  
  combined_results <- bind_rows(
    results_no_d %>% mutate(Model_Type = "no_d"),
    results_with_d %>% mutate(Model_Type = "with_d")
  )
  
  combined_results <- combined_results %>%
    mutate(
      Model_Label = case_when(
        Model_Type == "no_d" ~ paste0("SE_", gsub("model", "", Model)),
        Model_Type == "with_d" ~ paste0("SE_", gsub("model", "", Model), "_plus")
      )
    )
  
  
  residual_se_data <- combined_results %>%
    select(Language, Model_Label, Residual_SE)
  
  residual_se_wide <- residual_se_data %>%
    pivot_wider(
      names_from = Model_Label,
      values_from = Residual_SE
    ) %>%
    arrange(Language)
  
  model_order <- c(
    "Language",
    "SE_0",
    "SE_1", 
    "SE_2",
    "SE_3",
    "SE_4",
    "SE_5",
    "SE_1_plus",
    "SE_2_plus",
    "SE_3_plus",
    "SE_4_plus",
    "SE_5_plus"
  )
  
  existing_cols <- intersect(model_order, names(residual_se_wide))
  residual_se_wide <- residual_se_wide %>%
    select(all_of(existing_cols))
  
  print(residual_se_wide)
  
  return(residual_se_wide)
}


save_residual_se_table <- function(residual_se_table, filename = "residual_se_comparison.csv") {
  
  table_with_best <- highlight_best_se_models(residual_se_table)
  
  write.csv(table_with_best, filename, row.names = FALSE)
  cat(sprintf("\nResidual SE table saved to: %s\n", filename))
}
#########################################################################################################################################

extract_best_model_parameters <- function(results_no_d, results_with_d) {
  
  combined_results <- bind_rows(
    results_no_d %>% mutate(Model_Type = "no_d"),
    results_with_d %>% mutate(Model_Type = "with_d")
  )
  
  best_models <- combined_results %>%
    group_by(Language) %>%
    slice_min(AIC, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  params_table <- best_models %>%
    select(Language, Model, Description, Param_names, Param_values, Param_SEs, AIC, Residual_SE) %>%
    arrange(Language)
  
  params_expanded <- params_table %>%
    mutate(
      # Split param names, values, SEs
      param_list = strsplit(as.character(Param_values), ";"),
      se_list = strsplit(as.character(Param_SEs), ";"),
      name_list = strsplit(as.character(Param_names), ";")
    ) %>%
    rowwise() %>%
    mutate(
      params_formatted = list(
        if (!is.na(Param_names)) {
          paste0(name_list, " = ", param_list, " ± ", se_list)
        } else {
          NA_character_
        }
      )
    ) %>%
    ungroup() %>%
    select(Language, Model, Description, Param_names, Param_values, Param_SEs, AIC, Residual_SE)
  print(params_expanded)
  return(params_expanded)
}


#########################################################################################################################################

#########################################################################################################################################
#########################################################################################################################################