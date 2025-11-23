#plots
library(dplyr)
library(ggplot2)
library(patchwork)
########################################################################################

# LOG-LOG PLOTS GRID (raw data, non-aggregated)
plot.loglog.dependency.length.grid.raw <- function(lang_data) {
  
  languages <- names(lang_data$tables)
  plot_list <- list()
  
  for (lang_name in languages) {
    data <- lang_data$tables[[lang_name]]
    
    # Add random expectation to raw data
    data <- data %>%
      mutate(random_expectation = (vertices + 1) / 3)
    
    p <- ggplot(data, aes(x = vertices)) +
      geom_point(aes(y = mean_length),
                 color = "black", size = 0.5, alpha = 0.3) +
      geom_line(aes(y = random_expectation),
                color = "red", linewidth = 0.5, linetype = "dashed") +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = lang_name,
        x     = expression(n~"(number of vertices)"),
        y     = "Mean dependency length"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text  = element_text(size = 7)
      )
    
    plot_list[[lang_name]] <- p
  }
  
  combined <- wrap_plots(plot_list, ncol = 3, nrow = 7) +
    plot_annotation(
      title    = "Mean dependency length vs sentence length (raw data, log-log scale)",
      subtitle = "Black points: all observations | Red dashed: random expectation (n+1)/3"
    ) &
    labs(
      x = expression(n~"(number of vertices)"),
      y = "Mean dependency length"
    ) &
    theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title    = element_text(size = 11)
    )
  
  ggsave(
    "./loglog_grid_raw.png",
    plot   = combined,
    width  = 18,
    height = 24,
    dpi    = 150
  )
  cat("Grid plot saved to: ./loglog_grid_raw.png\n")
  
  return(combined)
}



plot.loglog.dependency.length.grid <- function(lang_data) {
  
  languages <- names(lang_data$tables)
  plot_list <- list()
  
  for (lang_name in languages) {
    data <- lang_data$tables[[lang_name]]
    
    mean_data <- data %>%
      group_by(vertices) %>%
      summarise(avg_mean_length = mean(mean_length), .groups = 'drop') %>%
      mutate(random_expectation = (vertices + 1) / 3)
    
    p <- ggplot(mean_data, aes(x = vertices)) +
      geom_point(aes(y = avg_mean_length), color = "black", size = 1.5) +
      geom_line(aes(y = avg_mean_length), color = "green", linewidth = 0.5) +
      geom_line(aes(y = random_expectation), color = "red", linewidth = 0.5, linetype = "dashed") +
      scale_x_log10() +
      scale_y_log10() +
      labs(title = lang_name) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_blank(),
        axis.text = element_text(size = 7)
      )
    
    plot_list[[lang_name]] <- p
  }
  
  combined <- wrap_plots(plot_list, ncol = 3, nrow = 7) +
    plot_annotation(
      title = "Mean Dependency Length vs Vertices (Aggregated)",
      subtitle = "Green: Observed | Red dashed: Random expectation (n+1)/3",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  ggsave("./loglog_grid.png", plot = combined, width = 18, height = 24, dpi = 150)
  cat("Grid plot saved to: ./loglog_grid.png\n")
  
  return(combined)
}

plot.binned.variance.grid <- function(lang_data, B = 1.2, min_count = 5, 
                                      alpha = NULL, plot_type = "variance") {
  
  languages <- names(lang_data$tables)
  plot_list <- list()
  
  for (lang_name in languages) {
    data <- lang_data$tables[[lang_name]]
    
    binned <- logarithmic.binning.with.merge(
      data = data,
      x_col = "vertices",
      y_col = "mean_length",
      B = B,
      min_count = min_count,
      alpha = alpha
    )
    
    binned_stats <- binned$binned_stats %>%
      mutate(y_var = y_sd^2)
    
    if (plot_type == "variance") {
      y_col <- "y_var"
      color_plot <- "darkred"
      y_label <- "Variance of mean dependency length"
    } else {
      y_col <- "y_sd"
      color_plot <- "darkblue"
      y_label <- "Standard deviation of mean dependency length"
    }
    
    p <- ggplot(binned_stats, aes(x = x_star, y = .data[[y_col]])) +
      geom_point(color = color_plot, size = 2) +
      geom_line(color = color_plot, alpha = 0.5, linewidth = 0.7) +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = lang_name,
        x     = expression(n^"*"~"(binned sentence length)"),
        y     = y_label
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text  = element_text(size = 7)
      )
    
    plot_list[[lang_name]] <- p
  }
  
  title_text <- ifelse(plot_type == "variance", "Variance", "Standard Deviation")
  
  combined <- wrap_plots(plot_list, ncol = 3, nrow = 7) +
    plot_annotation(
      title    = paste(title_text, "of binned mean dependency length"),
      subtitle = sprintf("B = %.2f | Min count = %d", B, min_count)
    ) &
    labs(
      x = expression(n^"*"~"(binned sentence length)"),
      y = y_label
    ) &
    theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title    = element_text(size = 11)
    )
  
  filename <- paste0("./", tolower(plot_type), "_grid.png")
  ggsave(filename, plot = combined, width = 18, height = 24, dpi = 150)
  cat(sprintf("Grid plot saved to: %s\n", filename))
  
  return(combined)
}


plot.binned.histogram.grid <- function(lang_data, B = 1.2, min_count = 5, alpha = NULL) {
  
  languages <- names(lang_data$tables)
  plot_list <- list()
  
  for (lang_name in languages) {
    data <- lang_data$tables[[lang_name]]
    
    binned <- logarithmic.binning.with.merge(
      data = data,
      x_col = "vertices",
      y_col = "mean_length",
      B = B,
      min_count = min_count,
      alpha = alpha
    )
    
    binned_stats <- binned$binned_stats
    
    p <- ggplot(binned_stats, aes(x = x_star, y = y_mean)) +
      geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se),
                    color = "darkblue", width = 0.02, linewidth = 0.5) +
      geom_point(color = "darkblue", size = 1) +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = lang_name,
        x     = expression(n^"*"~"(binned sentence length)"),
        y     = "Mean dependency length"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text  = element_text(size = 7)
      )
    
    plot_list[[lang_name]] <- p
  }
  
  combined <- wrap_plots(plot_list, ncol = 3, nrow = 7) +
    plot_annotation(
      title    = "Binned averages with standard error",
      subtitle = sprintf("B = %.2f | Min count = %d | Error bars: ± SE", B, min_count)
    ) &
    labs(
      x = expression(n^"*"~"(binned sentence length)"),
      y = "Mean dependency length"
    ) &
    theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title    = element_text(size = 11)
    )
  
  ggsave("./binned_histogram_grid.png", plot = combined, width = 18, height = 24, dpi = 150)
  cat("Grid plot saved to: ./binned_histogram_grid.png\n")
  
  return(combined)
}



############################################################################################################################
plot.loglog.raw.with.binned.grid <- function(lang_data,
                                             B = 1.18,
                                             min_count = 8,
                                             alpha = NULL) {
  
  languages <- names(lang_data$tables)
  plot_list <- list()
  
  for (lang_name in languages) {
    data <- lang_data$tables[[lang_name]]
    
    # Raw data + random expectation
    data <- data %>%
      mutate(random_expectation = (vertices + 1) / 3)
    
    # Logarithmic binning of mean_length vs vertices
    binned <- logarithmic.binning.with.merge(
      data      = data,
      x_col     = "vertices",
      y_col     = "mean_length",
      B         = B,
      min_count = min_count,
      alpha     = alpha
    )
    
    binned_stats <- binned$binned_stats
    
    p <- ggplot() +
      # Raw points (low opacity)
      geom_point(
        data = data,
        aes(x = vertices, y = mean_length, 
            color = "Raw data"),
        size  = 0.5,
        alpha = 0.1
      ) +
      # Random expectation on raw scale
      geom_line(
        data = data,
        aes(x = vertices, y = random_expectation, 
            color = "Null model"),
        linewidth = 0.4,
        linetype  = "dashed",
        alpha     = 0.7
      ) +
      # Binned means with error bars (± SE)
      geom_errorbar(
        data = binned_stats,
        aes(x = x_star,
            ymin = y_mean - y_se,
            ymax = y_mean + y_se,
            color = "Log-binned mean ± SE"),
        width     = 0.02,
        linewidth = 0.4
      ) +
      geom_point(
        data = binned_stats,
        aes(x = x_star, y = y_mean, 
            color = "Log-binned mean ± SE"),
        size = 1
      ) +
      scale_x_log10() +
      scale_y_log10() +
      scale_color_manual(
        name   = NULL,
        values = c(
          "Raw data"             = "black",
          "Log-binned mean ± SE" = "red",
          "Null model"           = "orange"
        )
      ) +
      labs(
        title = lang_name,
        x     = expression(n~"(number of vertices)"),
        y     = "Mean dependency length"
      ) +
      theme_minimal() +
      theme(
        plot.title    = element_text(size = 10, face = "bold"),
        axis.title    = element_text(size = 9),
        axis.text     = element_text(size = 7),
        legend.position = "bottom",
        legend.text   = element_text(size = 6)
      )
    
    plot_list[[lang_name]] <- p
  }
  
  combined <- wrap_plots(plot_list, ncol = 3, nrow = 7) +
    plot_annotation(
      title = "Mean dependency length vs sentence length (raw + log-binned, log-log scale)",
      subtitle = sprintf(
        "Black: raw data (low opacity) | Red: log-binned means ± SE | Orange dashed: random expectation (n+1)/3\nB = %.2f | Min count = %d",
        B, min_count
      )
    ) &
    labs(
      x = expression(n~"(number of vertices)"),
      y = "Mean dependency length"
    ) &
    theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title    = element_text(size = 11)
    )
  
  ggsave(
    "./loglog_raw_with_binned_grid.png",
    plot   = combined,
    width  = 18,
    height = 24,
    dpi    = 150
  )
  cat("Grid plot saved to: ./loglog_raw_with_binned_grid.png\n")
  
  return(combined)
}




############################################################################################################################
plot_all_best_models_grid <- function(lang_data, results_no_d, results_with_d,
                                      B = 1.2, min_count = 5, alpha = NULL) {
  
  combined_results <- bind_rows(
    results_no_d  %>% mutate(Model_Type = "no_d"),
    results_with_d %>% mutate(Model_Type = "with_d")
  )
  
  best_models <- combined_results %>%
    group_by(Language) %>%
    slice_min(AIC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Language)
  
  plot_list <- list()
  
  for (i in seq_len(nrow(best_models))) {
    lang_name  <- best_models$Language[i]
    model_desc <- best_models$Description[i]
    
    param_names  <- strsplit(as.character(best_models$Param_names[i]),  ";")[[1]]
    param_values <- as.numeric(strsplit(as.character(best_models$Param_values[i]), ";")[[1]])
    param_ses    <- as.numeric(strsplit(as.character(best_models$Param_SEs[i]),  ";")[[1]])
    
    if (!is.na(best_models$Param_names[i])) {
      param_lines <- paste0(param_names, " = ",
                            round(param_values, 4), " \u00b1 ",
                            round(param_ses, 4))
      param_text <- paste(param_lines, collapse = "\n")
    } else {
      param_text <- "No parameters"
    }
    
    data <- lang_data$tables[[lang_name]]
    binned <- logarithmic.binning.with.merge(
      data      = data,
      x_col     = "vertices",
      y_col     = "mean_length",
      B         = B,
      min_count = min_count,
      alpha     = alpha
    )
    
    binned_stats <- binned$binned_stats
    x    <- binned_stats$x_star
    y    <- binned_stats$y_mean
    y_se <- binned_stats$y_se  
    
    model_name <- best_models$Model[i]
    
    if (model_name == "model0") {
      y_fitted <- (x + 1) / 3
    } else if (model_name == "model1") {
      b <- param_values[1]
      d <- if (length(param_values) == 2) param_values[2] else 0
      y_fitted <- (x / 2)^b + d
    } else if (model_name == "model2") {
      a <- param_values[1]; b <- param_values[2]
      d <- if (length(param_values) == 3) param_values[3] else 0
      y_fitted <- a * x^b + d
    } else if (model_name == "model3") {
      a <- param_values[1]; c <- param_values[2]
      d <- if (length(param_values) == 3) param_values[3] else 0
      y_fitted <- a * exp(c * x) + d
    } else if (model_name == "model4") {
      a <- param_values[1]
      d <- if (length(param_values) == 2) param_values[2] else 0
      y_fitted <- a * log(x) + d
    } else if (model_name == "model5") {
      a <- param_values[1]; b <- param_values[2]; c <- param_values[3]
      d <- if (length(param_values) == 4) param_values[4] else 0
      y_fitted <- a * x^b * exp(c * x) + d
    }
    
    plot_data <- data.frame(
      x_star   = x, 
      observed = y, 
      fitted   = y_fitted,
      y_se     = y_se
    )
    
    x_range <- range(log10(x))
    y_range <- range(log10(y))
    
    p <- ggplot(plot_data, aes(x = x_star, y = observed)) +
      geom_errorbar(aes(ymin = observed - y_se, ymax = observed + y_se), 
                    width = 0.02, color = "black", linewidth = 0.3) +
      geom_point(size = 1.5, color = "black") +
      geom_line(aes(y = fitted), color = "red", linewidth = 0.5) +
      geom_text(
        data = data.frame(
          x     = 10^(x_range[1] + 0.05 * diff(x_range)),
          y     = 10^(y_range[2] - 0.05 * diff(y_range)),
          label = param_text
        ),
        aes(x = x, y = y, label = label),
        hjust = 0, vjust = 1,
        size  = 3, color = "black",
        lineheight = 0.9
      ) +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title    = lang_name,
        subtitle = model_desc,
        x        = expression(n^"*"~"(binned sentence length)"),
        y        = "Mean dependency length"
      ) +
      theme_minimal() +
      theme(
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8),
        axis.title    = element_text(size = 9),
        axis.text     = element_text(size = 7)
      )
    
    plot_list[[i]] <- p
  }
  
  combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 7) +
    plot_annotation(
      title    = "Best nonlinear model fits by language",
      subtitle = sprintf(
        "Logarithmic binning: B = %.2f, Min count = %d | Error bars: ± SE (SD/√n)",
        B, min_count
      )
    ) &
    labs(
      x = expression(n^"*"~"(binned sentence length)"),
      y = "Mean dependency length"
    ) &
    theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title    = element_text(size = 11)
    )
  
  ggsave("./best_models_grid.png",
         plot   = combined_plot,
         width  = 18,
         height = 24,
         dpi    = 150)
  
  cat("Grid plot saved to: ./best_models_grid.png\n")
  
  return(combined_plot)
}



library(ggrepel)

plot_mu_d_vs_mu_N <- function(lang_summary) {
  df <- lang_summary %>%
    rename(
      mu_N = mu_n,
      sd_N = sigma_n,
      mu_d = mu_d,
      sd_d = sigma_d
    )
  
  ggplot(df, aes(x = mu_N, y = mu_d, label = Language)) +
    geom_point(size = 1.5, color = "darkred") +
    geom_text_repel(
      size         = 1.5,
      max.overlaps = Inf,      
      box.padding  = 0.3,
      point.padding= 0.2,
    ) +
    labs(
      title = expression(paste("Mean dependency length ", mu[d],
                               " vs mean sentence length ", mu[N])),
      x     = expression(mu[N]~"(mean sentence length)"),
      y     = expression(mu[d]~"(mean dependency length)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 9)
    )
}


