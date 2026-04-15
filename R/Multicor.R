#' Multi-type Correlation Function
#'
#' Computes correlation matrices (Pearson, Spearman, Kendall, or Polychoric/Polyserial)
#' for numeric columns of a data frame, with p-values.
#'
#' @param data A data frame containing numeric (or ordinal for polycor) columns.
#' @param method Correlation method: "pearson", "spearman", "kendall", or "polycor".
#' @param use How to handle missing data: "everything", "complete.obs", "pairwise.complete.obs".
#'
#' @return A list containing:
#' \item{cor}{Correlation matrix.}
#' \item{p}{Matrix of p-values for each correlation.}
#'
#' @examples
#' df <- data.frame(
#'   x = rnorm(10),
#'   y = rnorm(10),
#'   z = rnorm(10)
#' )
#' multiCor(df, method = "pearson")
#'
#' @importFrom stats cor cor.test
#' @export
multiCor <- function(data, method = c("pearson","spearman","kendall","polycor"),
                     use = "pairwise.complete.obs") {

  method <- match.arg(method)

  # Only numeric columns for Pearson/Spearman/Kendall
  if (method != "polycor") {
    num_data <- data[, sapply(data, is.numeric), drop = FALSE]
    n <- ncol(num_data)
    cor_mat <- matrix(NA, n, n, dimnames = list(names(num_data), names(num_data)))
    p_mat <- matrix(NA, n, n, dimnames = list(names(num_data), names(num_data)))

    for (i in 1:n) {
      for (j in i:n) {
        test <- cor.test(num_data[[i]], num_data[[j]], method = method, use = use)
        cor_mat[i,j] <- cor_mat[j,i] <- test$estimate
        p_mat[i,j] <- p_mat[j,i] <- test$p.value
      }
    }

    return(list(cor = cor_mat, p = p_mat))
  }

  # Polychoric/polyserial correlation requires polycor package
  if (method == "polycor") {
    if (!requireNamespace("polycor", quietly = TRUE))
      stop("Package 'polycor' needed for polycor correlations. Install it first.")

    # polycor::hetcor works for numeric + ordinal/factor
    cor_obj <- polycor::hetcor(data)
    cor_mat <- cor_obj$correlations
    # Note: p-values not provided by polycor
    return(list(cor = cor_mat, p = NULL))
  }
}
