#' Latin Square Analysis (LS Design)
#'
#' Perform ANOVA, assumption checks, and mean separation for Latin Square experiments.
#'
#' @param data A data frame containing the experimental data.
#' @param trait Character vector of trait column names to analyze.
#' @param treatment Column name for the treatment factor.
#' @param row Column name for row factor.
#' @param column Column name for column factor.
#' @param var_test Variance homogeneity test: "bartlett" (default) or "levene".
#' @param norm_test Normality test: "shapiro" (default), "ks", or "ad".
#' @param mean_method Method for mean separation: "lsd" (default), "tukey", "duncan", "snk", "hsd".
#' @param check_normality Logical, whether to perform normality test (default TRUE).
#' @param check_variance Logical, whether to perform variance homogeneity test (default TRUE).
#' @param alpha Significance level (default 0.05).
#'
#' @return A list containing ANOVA model, ANOVA table, assumption tests, and mean separation for each trait.
#' @export
#'
#' @examples
#' \dontrun{
#' library(WakjiraAgriStats)
#' data("sweetpotato")
#' res_ls <- LS(
#'   data = sweetpotato,
#'   trait = c("yield", "PlantHeight"),
#'   treatment = "virus",
#'   row = "Row",
#'   column = "Column",
#'   check_normality = TRUE,
#'   check_variance = TRUE,
#'   mean_method = "duncan"
#' )
#' }
LS <- function(data, trait, treatment, row, column,
               var_test = "bartlett",
               norm_test = "shapiro",
               mean_method = "lsd",
               check_normality = TRUE,
               check_variance = TRUE,
               alpha = 0.05) {

  results <- list()

  for(tr in trait){
    cat("\n=====================================\n")
    cat("LATIN SQUARE ANALYSIS — Trait:", tr, "\n")

    formula <- as.formula(paste(tr, "~", treatment, "+", row, "+", column))
    model <- aov(formula, data = data)
    cat("\n--- ANOVA ---\n")
    print(summary(model))

    res <- residuals(model)

    # Normality
    if(check_normality){
      normality <- switch(norm_test,
                          shapiro = if(length(res) >=3 & length(res) <= 5000) shapiro.test(res) else "Too few residuals",
                          ks     = ks.test(scale(res), "pnorm"),
                          ad     = nortest::ad.test(res),
                          "Invalid normality test")
      cat("\n--- Normality Test ---\n")
      print(normality)
    } else normality <- "Skipped"

    # Variance homogeneity
    if(check_variance){
      var_test_result <- switch(var_test,
                                bartlett = bartlett.test(res ~ data[[treatment]]),
                                levene   = car::leveneTest(res ~ data[[treatment]]),
                                "Invalid variance test")
      cat("\n--- Homogeneity of Variance ---\n")
      print(var_test_result)
    } else var_test_result <- "Skipped"

    # Mean separation
    cat("\n--- Mean Separation ---\n")
    ms <- mean_separation(model, method = mean_method, alpha = alpha)
    print(ms)

    results[[tr]] <- list(
      model = model,
      anova = summary(model),
      normality = normality,
      variance_test = var_test_result,
      mean_sep = ms
    )
  }

  invisible(results)
}
