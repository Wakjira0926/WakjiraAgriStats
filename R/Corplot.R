#' Advanced Interactive Triangular Correlation Heatmap (Robust Version)
#'
#' Generates a triangular correlation heatmap with numeric values, significance stars,
#' hierarchical clustering (optional), custom colors, and interactive Plotly if available.
#'
#' @param data A data frame with numeric columns.
#' @param method Correlation method: "pearson", "spearman", "kendall".
#' @param sig.level Significance thresholds for stars (default: c(0.05,0.01,0.001)).
#' @param title Plot title.
#' @param colors Optional vector of colors for gradient (default c("red","white","green")).
#' @param ncolors Number of color steps (default 100).
#' @param cluster Logical. If TRUE, reorder variables by hierarchical clustering (default FALSE).
#' @param interactive Logical. If TRUE, generate interactive Plotly plot if available (default TRUE).
#' @param theme_custom A ggplot2 theme to override default theme (default minimal + clean).
#' @param font.size Font size for axis labels (default 12).
#'
#' @return A ggplot2 object (if interactive = FALSE or plotly unavailable) or a plotly object (if interactive = TRUE and plotly installed).
#' @export
advancedInteractiveCorr <- function(data, method = c("pearson","spearman","kendall"),
                                    sig.level = c(0.05,0.01,0.001),
                                    title = "Triangular Correlation Heatmap",
                                    colors = NULL, ncolors = 100,
                                    cluster = FALSE, interactive = TRUE,
                                    theme_custom = NULL, font.size = 12) {

  library(ggplot2)
  library(reshape2)

  # Check plotly availability
  plotly_available <- requireNamespace("plotly", quietly = TRUE)

  method <- match.arg(method)
  num_data <- data[, sapply(data, is.numeric), drop = FALSE]
  vars <- names(num_data)
  n <- length(vars)

  # Compute correlation and p-values
  cor_mat <- matrix(NA, n, n, dimnames = list(vars, vars))
  p_mat <- matrix(NA, n, n, dimnames = list(vars, vars))

  for (i in 1:n) {
    for (j in i:n) {
      test <- cor.test(num_data[[i]], num_data[[j]], method = method)
      cor_mat[i,j] <- cor_mat[j,i] <- test$estimate
      p_mat[i,j] <- p_mat[j,i] <- test$p.value
    }
  }

  # Hierarchical clustering
  if (cluster) {
    dist_mat <- as.dist(1 - cor_mat)
    hc <- hclust(dist_mat)
    cor_mat <- cor_mat[hc$order, hc$order]
    p_mat <- p_mat[hc$order, hc$order]
    vars <- vars[hc$order]
  }

  # Melt matrices
  cor_df <- melt(cor_mat, na.rm = TRUE)
  names(cor_df) <- c("Var1","Var2","Correlation")
  p_df <- melt(p_mat, na.rm = TRUE)
  cor_df$p <- p_df$value

  # Lower triangle only
  cor_df <- cor_df[as.numeric(factor(cor_df$Var1, levels=vars)) >
                     as.numeric(factor(cor_df$Var2, levels=vars)), ]

  # Significance stars
  cor_df$Signif <- ""
  cor_df$Signif[cor_df$p <= sig.level[1]] <- "*"
  cor_df$Signif[cor_df$p <= sig.level[2]] <- "**"
  cor_df$Signif[cor_df$p <= sig.level[3]] <- "***"

  # Color gradient
  if (is.null(colors)) colors <- c("red","white","green")
  fill_scale <- scale_fill_gradientn(colors = colorRampPalette(colors)(ncolors),
                                     limits = c(min(cor_df$Correlation), max(cor_df$Correlation)),
                                     name = "Correlation")

  # Default theme
  if (is.null(theme_custom)) {
    theme_custom <- theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = font.size),
            axis.text.y = element_text(size = font.size),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(face = "bold", size = font.size+2))
  }

  # Base ggplot
  p <- ggplot(cor_df, aes(x = Var2, y = Var1, fill = Correlation,
                          text = paste0("Var1: ", Var1,
                                        "<br>Var2: ", Var2,
                                        "<br>Correlation: ", round(Correlation,2),
                                        "<br>p-value: ", signif(p,3),
                                        "<br>Signif: ", Signif))) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(round(Correlation,2), Signif)), size = 4) +
    fill_scale +
    scale_y_discrete(limits = rev(vars)) +
    ggtitle(title) +
    coord_fixed() +
    theme_custom

  # Interactive Plotly if available
  if (interactive) {
    if (plotly_available) {
      return(plotly::ggplotly(p, tooltip = "text"))
    } else {
      message("Plotly is not installed. Returning static ggplot2 plot instead.")
      return(p)
    }
  } else {
    return(p)
  }
}
