#' Row–Column Design Analysis (lmer-based)
#'
#' Performs replicated Row–Column analysis using a mixed model with
#' rows and columns nested within replication. Supports multiple traits,
#' diagnostics, and mean separation.
#'
#' @param data Data frame containing the experiment.
#' @param traits Character vector of trait names.
#' @param normality "shapiro", "ad", "lillie", or "all".
#' @param variance "levene", "bartlett", "fligner", or "all".
#' @param mean_sep "lsd", "tukey", or "all".
#'
#' @return A list of results per trait.
#' @export
RCD<- function(data,
               traits,
               normality = "shapiro",
               variance = "bartlett",
               mean_sep = "tukey") {

  suppressMessages({
    library(lme4)
    library(lmerTest)
    library(emmeans)
    library(multcomp)
    library(multcompView)
    library(car)
    library(nortest)
    library(dplyr)
    library(stringr)
  })

  # ---------------------------
  # COLUMN DETECTION
  # ---------------------------
  col_lower <- tolower(names(data))
  rep_col   <- names(data)[col_lower %in% c("rep","replication")][1]
  row_col   <- names(data)[col_lower == "row"][1]
  col_col   <- names(data)[col_lower %in% c("col","column")][1]
  treat_col <- names(data)[col_lower %in% c("gen","genotype","treat","treatment")][1]

  if (any(is.na(c(rep_col, row_col, col_col, treat_col))))
    stop("Required columns missing: REP, ROW, COL, TREAT/GEN.")

  # ---------------------------
  # FACTOR CONVERSION
  # ---------------------------
  data <- data %>%
    mutate(
      rep   = as.factor(.data[[rep_col]]),
      row   = as.factor(.data[[row_col]]),
      col   = as.factor(.data[[col_col]]),
      treat = as.factor(.data[[treat_col]])
    )

  results <- list()

  # ---------------------------
  # TRAIT LOOP
  # ---------------------------
  for (trait in traits) {

    cat("\n========================================\n")
    cat("Row–Column Analysis:", trait, "\n")
    cat("========================================\n\n")

    # ---------------------------
    # MIXED MODEL
    # ---------------------------
    model_formula <- as.formula(
      paste(trait,
            "~ treat + rep +",
            "(1 | rep:row) +",
            "(1 | rep:col)")
    )

    mod <- lmer(model_formula, data = data)

    # ---------------------------
    # VARIANCE COMPONENTS (FIXED SELECT ISSUE)
    # ---------------------------
    vc <- VarCorr(mod) %>%
      as.data.frame() %>%
      dplyr::select(grp, vcov)

    vc$grp <- vc$grp %>%
      str_replace_all("rep:row", "rep:row") %>%
      str_replace_all("rep:col", "rep:col") %>%
      str_replace_all("Residual", "Residual")

    cat("Variance Components\n")
    print(vc)
    cat("\n")

    # ---------------------------
    # ANOVA (Kenward–Roger)
    # ---------------------------
    anova_tbl <- anova(mod, ddf = "Kenward-Roger")
    cat("ANOVA Table\n")
    print(anova_tbl)
    cat("\n")

    # ---------------------------
    # RESIDUALS & NORMALITY TESTS
    # ---------------------------
    res <- residuals(mod)
    normality_results <- list()

    if (normality %in% c("shapiro","all"))
      normality_results$Shapiro <- shapiro.test(res)
    if (normality %in% c("ad","all"))
      normality_results$AndersonDarling <- nortest::ad.test(res)
    if (normality %in% c("lillie","all"))
      normality_results$Lilliefors <- nortest::lillie.test(res)

    cat("Normality Tests (p-values)\n")
    print(sapply(normality_results, `[[`, "p.value"))
    cat("\n")

    # ---------------------------
    # HOMOGENEITY OF VARIANCE
    # ---------------------------
    variance_results <- list()

    if (variance %in% c("levene","all"))
      variance_results$Levene <- tryCatch(
        car::leveneTest(data[[trait]] ~ data$treat),
        error = function(e) NA
      )

    if (variance %in% c("bartlett","all"))
      variance_results$Bartlett <- tryCatch(
        bartlett.test(data[[trait]] ~ data$treat),
        error = function(e) NA
      )

    if (variance %in% c("fligner","all"))
      variance_results$Fligner <- tryCatch(
        fligner.test(data[[trait]] ~ data$treat),
        error = function(e) NA
      )

    cat("Homogeneity of Variance (p-values)\n")
    for (v in names(variance_results)) {
      test_result <- variance_results[[v]]
      if (is.list(test_result)) {
        if ("Pr(>F)" %in% names(test_result))
          cat(v, "p-value:", test_result$`Pr(>F)`[1], "\n")
        else if ("p.value" %in% names(test_result))
          cat(v, "p-value:", test_result$p.value, "\n")
        else
          cat(v, "p-value: NA\n")
      } else {
        cat(v, "p-value: NA\n")
      }
    }
    cat("\n")

    # ---------------------------
    # COEFFICIENT OF VARIATION
    # ---------------------------
    resid_var <- attr(VarCorr(mod), "sc")^2
    mean_trait <- mean(data[[trait]], na.rm = TRUE)
    cv <- sqrt(resid_var) / mean_trait * 100

    cat("CV (%):", round(cv, 2), "\n\n")

    # ---------------------------
    # MEAN SEPARATION (Fixed CLD Issue)
    # ---------------------------
    emm <- emmeans(mod, specs = "treat", lmer.df = "kenward-roger")
    mean_results <- list()

    if (mean_sep %in% c("lsd","all"))
      mean_results$LSD <- multcomp::cld(emm, adjust = "none", Letters = letters)

    if (mean_sep %in% c("tukey","all"))
      mean_results$Tukey <- multcomp::cld(emm, adjust = "tukey", Letters = letters)

    cat("Mean Separation Results\n")
    for (m in names(mean_results)) {
      cat("----", m, "----\n")
      print(mean_results[[m]]$groups)
      cat("\n")
    }

    # ---------------------------
    # STORE RESULTS
    # ---------------------------
    results[[trait]] <- list(
      Model      = mod,
      VarComp    = vc,
      ANOVA      = anova_tbl,
      CV         = cv,
      Normality  = normality_results,
      Variance   = variance_results,
      Emmeans    = emm,
      Groups     = mean_results
    )
  }

  return(results)
}
