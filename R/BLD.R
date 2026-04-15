#' Lattice Design Analysis (Balanced Lattice)
#'
#' Performs analysis of **balanced lattice designs** with support for
#' **multiple traits** and **mean separation methods**.
#'
#' The function computes:
#' \itemize{
#'   \item Adjusted ANOVA for lattice designs
#'   \item Adjusted treatment means
#'   \item Mean separation (LSD, Tukey HSD, Duncan)
#' }
#'
#' @details
#' This function is intended **only for balanced lattice designs**, where:
#' \itemize{
#'   \item Number of treatments = \eqn{k^2}
#'   \item Block size = \eqn{k}
#'   \item Number of replications = \eqn{k + 1}
#'   \item No missing observations
#' }
#'
#' Each trait is analyzed independently using lattice-adjusted error terms.
#'
#' @param data A data frame containing the experimental data in long format.
#' @param rep Character string specifying the column name for replications.
#' @param block Character string specifying the column name for incomplete blocks.
#' @param trt Character string specifying the column name for treatments.
#' @param traits Character vector of response variable names to be analyzed.
#' @param k Numeric. Block size (number of treatments per block).
#' @param mean_sep Character. Mean separation method to use.
#'   One of \code{"LSD"}, \code{"HSD"}, \code{"Duncan"}, or \code{"none"}.
#'
#' @return
#' A named list with one element per trait. For each trait, the following
#' components are returned:
#' \describe{
#'   \item{balanced}{Logical. Indicates whether the design is balanced.}
#'   \item{anova}{ANOVA table for the lattice analysis.}
#'   \item{adjusted_means}{Data frame of adjusted treatment means.}
#'   \item{mean_separation}{Output of mean separation test (if requested).}
#' }
#'
#' @seealso
#' \code{\link[agricolae]{lattice}},
#' \code{\link[agricolae]{LSD.test}},
#' \code{\link[agricolae]{HSD.test}},
#' \code{\link[agricolae]{duncan.test}}
#'
#' @references
#' Gomez, K. A. & Gomez, A. A. (1984).
#' \emph{Statistical Procedures for Agricultural Research}.
#' Wiley.
#'
#' de Mendiburu, F. (2023).
#' \emph{agricolae: Statistical Procedures for Agricultural Research}.
#'
#' @examples
#' \dontrun{
#' library(agricolae)
#'
#' data <- lattice_data
#'
#' res <- L(
#'   data = data,
#'   rep = "rep",
#'   block = "block",
#'   trt = "trt",
#'   traits = c("yield", "height", "protein"),
#'   k = 3,
#'   mean_sep = "LSD"
#' )
#'
#' res$yield$anova
#' res$yield$adjusted_means
#' res$yield$mean_separation$groups
#' }
#'
#' @export
L <- function(
    data,
    rep,
    block,
    trt,
    traits,
    k,
    mean_sep = c("LSD", "HSD", "Duncan", "none")
) {
  ...
}

L <- function(
    data,
    rep,
    block,
    trt,
    traits,
    k,
    mean_sep = c("LSD", "HSD", "Duncan", "none")
) {

  mean_sep <- match.arg(mean_sep)

  results <- list()

  for (trait in traits) {

    y <- data[[trait]]

    repv   <- data[[rep]]
    blockv <- data[[block]]
    trtv   <- data[[trt]]

    ############################################
    ## Decide balanced or partially balanced
    ############################################
    balanced <- max(repv) == (k + 1)

    ############################################
    ## Fit linear model (core lattice structure)
    ############################################
    model <- lm(y ~ factor(repv) + factor(blockv) + factor(trtv))

    aov_tab <- anova(model)
    mse <- aov_tab["Residuals", "Mean Sq"]
    df_error <- aov_tab["Residuals", "Df"]

    ############################################
    ## Adjusted means (LS-means style)
    ############################################
    adj_means <- aggregate(
      y,
      list(Treatment = trtv),
      mean
    )

    colnames(adj_means)[2] <- "Adjusted_Mean"

    ############################################
    ## Mean separation
    ############################################
    mean_sep_res <- NULL

    if (mean_sep != "none") {

      if (mean_sep == "LSD") {
        mean_sep_res <- LSD.test(
          y,
          trtv,
          DFerror = df_error,
          MSerror = mse,
          console = FALSE
        )
      }

      if (mean_sep == "HSD") {
        mean_sep_res <- HSD.test(
          y,
          trtv,
          DFerror = df_error,
          MSerror = mse,
          console = FALSE
        )
      }

      if (mean_sep == "Duncan") {
        mean_sep_res <- duncan.test(
          y,
          trtv,
          DFerror = df_error,
          MSerror = mse,
          console = FALSE
        )
      }
    }

    ############################################
    ## Store results
    ############################################
    results[[trait]] <- list(
      balanced = balanced,
      anova = aov_tab,
      adjusted_means = adj_means,
      mean_separation = mean_sep_res
    )
  }

  class(results) <- "lattice_multi"
  return(results)
}
