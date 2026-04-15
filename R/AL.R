#' Alpha Lattice Design Analysis
#'
#' This function generates an Alpha Lattice (incomplete block) design layout
#' for field experiments or analyzes existing data using mixed models, ANOVA,
#' residual analysis, normality tests, variance homogeneity tests, and mean separation.
#'
#' @param treatments Either a vector of treatment names (character) or a numeric
#'   value specifying the number of treatments (e.g., 10). If numeric, treatments
#'   are automatically named as "T1", "T2", etc.
#' @param block_size Number of plots per block.
#' @param reps Number of replicates. Default is 2.
#' @param data A data frame containing experimental data. Must include columns
#'   for replication (REP), block (BLOCK), and genotype/treatment (GENOTYPE).
#'   If NULL, the function generates a design layout only.
#' @param traits A character vector of column names in `data` representing the traits
#'   to analyze. Required if `data` is provided.
#' @param normality A character string specifying the normality test to perform on
#'   residuals. Options: "shapiro", "ad" (Anderson-Darling), "lillie" (Lilliefors),
#'   or "all" for all tests. Default is "shapiro".
#' @param variance A character string specifying the test for homogeneity of variance.
#'   Options: "levene", "bartlett", "fligner", or "all". Default is "bartlett".
#' @param mean_sep A character string specifying the mean separation method. Options:
#'   "lsd" (Least Significant Difference), "tukey" (Tukey HSD), "duncan" (Duncan test),
#'   "snk" (Student-Newman-Keuls), or "all". Default is "lsd".
#' @param seed Optional numeric value for reproducible randomization.
#'
#' @return
#' If `data` is NULL, returns a list with the design layout:
#' \item{Design}{A data frame containing REP, BLOCK, PLOT, and GENOTYPE.}
#'
#' If `data` is provided, returns a list of results for each trait:
#' \itemize{
#'   \item \code{Model}: Mixed model object (nlme::lme)
#'   \item \code{VarComp}: Variance components
#'   \item \code{ANOVA}: ANOVA table
#'   \item \code{CV}: Coefficient of variation (%)
#'   \item \code{Normality}: Normality test results
#'   \item \code{Variance}: Homogeneity of variance test results
#'   \item \code{Means}: Mean values per treatment
#'   \item \code{Groups}: Treatment groupings from mean separation
#' }
#'
#' @details
#' - If `data` is NULL, the function generates a randomized Alpha Lattice design
#'   layout suitable for field experiments.
#' - If `data` is provided, the function:
#'   1. Fits a linear mixed model using `nlme::lme` with treatments as fixed effects
#'      and replicates/blocks as random effects.
#'   2. Computes variance components.
#'   3. Performs ANOVA for treatment effects.
#'   4. Computes residuals and performs normality tests.
#'   5. Performs homogeneity of variance tests.
#'   6. Computes coefficient of variation (CV).
#'   7. Conducts mean separation according to the selected method(s).
#'
#' @examples
#' # Generate a design layout for 12 treatments in blocks of 4 with 3 replicates
#' AL(treatments = 12, block_size = 4, reps = 3)
#'
#' # Analyze existing data
#' # Assuming 'mydata' has columns REP, BLOCK, GENOTYPE, Yield, Height
#' AL(treatments = NULL, block_size = 4, data = mydata, traits = c("Yield", "Height"),
#'    normality = "all", variance = "all", mean_sep = "all")
#'
#' @import nlme agricolae car nortest
#' @export
AL <- function(treatments,
               block_size,
               reps = 2,
               data = NULL,
               traits = NULL,
               normality = "shapiro",
               variance = "bartlett",
               mean_sep = "lsd",
               seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # -----------------------------
  # DESIGN MODE: Generate layout
  # -----------------------------
  if (is.null(data)) {
    if (length(treatments) == 1 && is.numeric(treatments))
      treatments <- paste0("T", 1:treatments)

    t <- length(treatments)
    if (t %% block_size != 0)
      stop("Treatments must be divisible by block_size")

    blocks_per_rep <- t / block_size
    layout <- data.frame()

    for (r in 1:reps) {
      shuffled <- sample(treatments)
      rep_blocks <- split(shuffled, rep(1:blocks_per_rep, each = block_size))
      rep_df <- data.frame(
        REP = r,
        BLOCK = rep(1:blocks_per_rep, each = block_size),
        PLOT = 1:t,
        GENOTYPE = unlist(rep_blocks)
      )
      layout <- rbind(layout, rep_df)
    }

    return(list(Design = layout))
  }

  # -----------------------------
  # ANALYSIS MODE: Load packages
  # -----------------------------
  suppressMessages({
    library(nlme)
    library(agricolae)
    library(car)
    library(nortest)
  })

  # -----------------------------
  # Detect columns
  # -----------------------------
  col_lower <- tolower(names(data))
  rep_col   <- names(data)[col_lower %in% c("rep","replication")][1]
  block_col <- names(data)[col_lower == "block"][1]
  treat_col <- names(data)[col_lower %in% c("genotype","treatment")][1]

  if (any(is.na(c(rep_col, block_col, treat_col))))
    stop("Required columns missing: REP, BLOCK, GENOTYPE")

  results <- list()

  # -----------------------------
  # Loop through traits
  # -----------------------------
  for (trait in traits) {

    cat("\n====================================================\n")
    cat("\nAlpha lattice analysis:", trait, "\n")
    cat("====================================================\n\n")
    cat("Estimation Method: REML\n\n")

    # Mixed model
    model <- nlme::lme(
      fixed = as.formula(paste(trait, "~", treat_col)),
      random = as.formula(paste("~1 |", rep_col, "/", block_col)),
      data = data,
      method = "REML"
    )

    vc <- VarCorr(model)

    cat("Parameter Estimates\n"); print(vc); cat("\n")

    # ANOVA
    cat("Analysis of Variance Table\n")
    anova_tbl <- anova(model)
    print(anova_tbl); cat("\n")

    # Residuals
    res <- residuals(model, type = "pearson")

    # Normality tests
    normality_results <- list()
    if (normality %in% c("shapiro","all")) normality_results$Shapiro <- shapiro.test(res)
    if (normality %in% c("ad","all")) normality_results$AndersonDarling <- nortest::ad.test(res)
    if (normality %in% c("lillie","all")) normality_results$Lilliefors <- nortest::lillie.test(res)

    cat("Normality Tests (p-values)\n")
    print(lapply(normality_results, function(x) x$p.value)); cat("\n")

    # Variance tests
    variance_results <- list()
    if (variance %in% c("levene","all")) variance_results$Levene <- car::leveneTest(data[[trait]] ~ data[[treat_col]])
    if (variance %in% c("bartlett","all")) variance_results$Bartlett <- bartlett.test(data[[trait]] ~ data[[treat_col]])
    if (variance %in% c("fligner","all")) variance_results$Fligner <- fligner.test(data[[trait]] ~ data[[treat_col]])

    cat("Variance Homogeneity Tests (p-values)\n")
    print(lapply(variance_results, function(x) x$p.value)); cat("\n")

    # Coefficient of variation
    resid_var <- as.numeric(vc["Residual","Variance"])
    mean_trait <- mean(data[[trait]], na.rm = TRUE)
    cv <- (sqrt(resid_var)/mean_trait)*100
    cat("Coefficient of Variation:", round(cv,2), "%\n")
    cat("Trait Mean:", mean_trait, "\n\n")

    # Fixed model for mean separation
    fixed_model <- aov(as.formula(paste(trait, "~", treat_col, "+", rep_col, "/", block_col)), data = data)

    # Mean separation
    mean_results <- list()
    if (mean_sep %in% c("lsd","all")) mean_results$LSD <- agricolae::LSD.test(fixed_model, treat_col, group=TRUE)
    if (mean_sep %in% c("tukey","all")) mean_results$Tukey <- agricolae::HSD.test(fixed_model, treat_col, group=TRUE)
    if (mean_sep %in% c("duncan","all")) mean_results$Duncan <- agricolae::duncan.test(fixed_model, treat_col, group=TRUE)
    if (mean_sep %in% c("snk","all")) mean_results$SNK <- agricolae::SNK.test(fixed_model, treat_col, group=TRUE)

    cat("Mean Separation Results\n")
    for (m in names(mean_results)) {
      cat("----", m, "----\n")
      print(mean_results[[m]]$groups); cat("\n")
    }

    # Store results
    results[[trait]] <- list(
      Model = model,
      VarComp = vc,
      ANOVA = anova_tbl,
      CV = cv,
      Normality = normality_results,
      Variance = variance_results,
      Means = lapply(mean_results, function(x) x$means),
      Groups = lapply(mean_results, function(x) x$groups)
    )
  }

  return(results)
}
