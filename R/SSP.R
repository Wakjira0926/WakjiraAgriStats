SSP <- function(data, responses, Replication, A, B, C,
                alpha = 0.05,
                norm_test = "shapiro",
                homogeneity = "bartlett",
                mean_sep_test = "LSD") {

  library(agricolae)
  library(car)
  library(nortest)

  # Convert to factors
  data[[Replication]] <- as.factor(data[[Replication]])
  data[[A]] <- as.factor(data[[A]])
  data[[B]] <- as.factor(data[[B]])
  data[[C]] <- as.factor(data[[C]])

  # Create interaction terms
  data$AB <- interaction(data[[A]], data[[B]])
  data$AC <- interaction(data[[A]], data[[C]])
  data$BC <- interaction(data[[B]], data[[C]])
  data$ABC <- interaction(data[[A]], data[[B]], data[[C]])

  results <- list()

  for (response in responses) {
    # --- Correct split-split-plot model
    fmla <- as.formula(
      paste(response, "~", Replication, "+", A, "*", B, "*", C,
            "+ Error(", Replication, "/", A, "/", B, ")")
    )
    model <- aov(fmla, data = data)
    anova_list <- summary(model)

    # --- Extract CVs
    MSEa <- anova_list[[2]][[1]]["Residuals", "Mean Sq"]
    MSEb <- anova_list[[3]][[1]]["Residuals", "Mean Sq"]
    MSEc <- anova_list[[4]][[1]]["Residuals", "Mean Sq"]

    grandmean <- mean(data[[response]], na.rm = TRUE)
    CVa <- sqrt(MSEa) / grandmean * 100
    CVb <- sqrt(MSEb) / grandmean * 100
    CVc <- sqrt(MSEc) / grandmean * 100

    # --- Combined ANOVA table (like your desired format)
    lm_model <- lm(as.formula(paste(response, "~", Replication, "+", A, "*", B, "*", C)), data = data)
    anova_combined <- anova(lm_model) # Standard compact ANOVA table

    # --- Normality test
    resid_vals <- residuals(lm_model)
    norm_result <- switch(norm_test,
                          shapiro = shapiro.test(resid_vals),
                          ad = ad.test(resid_vals),
                          ks = ks.test(resid_vals, "pnorm",
                                       mean(resid_vals), sd(resid_vals)))

    # --- Homogeneity test
    inter_fac <- interaction(data[[A]], data[[B]], data[[C]])
    hom_result <- switch(homogeneity,
                         bartlett = bartlett.test(data[[response]], inter_fac),
                         levene = leveneTest(data[[response]] ~ inter_fac),
                         fligner = fligner.test(data[[response]], inter_fac))

    # --- Mean separation tests on compact model factors
    mean_sep <- list()
    run_sep <- function(fac) {
      sep_formula <- as.formula(paste(response, "~", fac))
      tmp <- aov(sep_formula, data = data)
      switch(mean_sep_test,
             LSD = LSD.test(tmp, fac, alpha = alpha, group = TRUE),
             HSD = HSD.test(tmp, fac, alpha = alpha, group = TRUE),
             SNK = SNK.test(tmp, fac, alpha = alpha, group = TRUE),
             Duncan = duncan.test(tmp, fac, alpha = alpha, group = TRUE),
             Tukey = HSD.test(tmp, fac, alpha = alpha, group = TRUE))
    }

    # Run separations
    mean_sep[[A]] <- run_sep(A)
    mean_sep[[B]] <- run_sep(B)
    mean_sep[[C]] <- run_sep(C)
    mean_sep$AB <- run_sep("AB")
    mean_sep$AC <- run_sep("AC")
    mean_sep$BC <- run_sep("BC")
    mean_sep$ABC <- run_sep("ABC")

    # --- Store results
    results[[response]] <- list(
      trait = response,
      design = "Split-Split Plot",
      anova_full = anova_list,      # The full (multi-error) version
      anova_combined = anova_combined,  # Your desired compact version
      CV = c(CVa = CVa, CVb = CVb, CVc = CVc),
      mean_value = grandmean,
      normality_test = norm_result,
      homogeneity_test = hom_result,
      mean_separation = mean_sep,
      model = model
    )
  }

  return(results)
}
#' Split-Split Plot Analysis for Multiple Traits
#'
#' Performs split-split plot ANOVA for one or multiple response variables,
#' computes coefficients of variation (CV) for main plots, subplots, and sub-subplots,
#' tests normality and homogeneity of variances, and conducts mean separation for
#' main effects and interactions.
#'
#' @param data A data frame containing the experimental data.
#' @param responses A character vector of response variable names.
#' @param Replication Name of the replication/block factor (string).
#' @param A Name of the main plot factor (string).
#' @param B Name of the subplot factor (string).
#' @param C Name of the sub-subplot factor (string).
#' @param alpha Significance level for ANOVA and mean separation tests (default 0.05).
#' @param norm_test Normality test to use. Options: "shapiro" (default), "ad", "ks".
#' @param homogeneity Test for homogeneity of variance. Options: "bartlett" (default), "levene", "fligner".
#' @param mean_sep_test Method for mean separation. Options: "LSD" (default), "HSD", "SNK", "Duncan", "Tukey".
#'
#' @return A list containing results for each response variable, including:
#' \describe{
#'   \item{trait}{Name of the response variable.}
#'   \item{design}{Design type: "Split-Split Plot".}
#'   \item{anova_full}{Full ANOVA table using multi-error terms.}
#'   \item{anova_combined}{Compact ANOVA table with standard F-tests.}
#'   \item{CV}{Named vector of coefficients of variation: CVa, CVb, CVc.}
#'   \item{mean_value}{Grand mean of the response variable.}
#'   \item{normality_test}{Result of the selected normality test on residuals.}
#'   \item{homogeneity_test}{Result of the selected homogeneity of variance test.}
#'   \item{mean_separation}{List of mean separation results for main effects (A, B, C) and interactions (AB, AC, BC, ABC).}
#'   \item{model}{Fitted aov object.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage
#' results <- SSP(
#'   data = mydata,
#'   responses = c("yield", "height"),
#'   Replication = "block",
#'   A = "nitrogen",
#'   B = "management",
#'   C = "variety",
#'   alpha = 0.05,
#'   norm_test = "shapiro",
#'   homogeneity = "bartlett",
#'   mean_sep_test = "HSD"
#' )
#'
#' # Access ANOVA for yield
#' results$yield$anova_combined
#'
#' # Access interaction mean separation
#' results$yield$mean_separation$AB
#' }
#'
#' @export
SSP <- function(data, responses, Replication, A, B, C,
                alpha = 0.05,
                norm_test = "shapiro",
                homogeneity = "bartlett",
                mean_sep_test = "LSD") {
  # ... [Your existing SSP function code here]
}
