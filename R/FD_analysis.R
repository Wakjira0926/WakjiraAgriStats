#' Factorial Design ANOVA Analysis
#'
#' Perform ANOVA for factorial designs (CRD or RCBD) including normality tests,
#' homogeneity of variance, coefficient of variation, and mean separation for
#' main effects and interactions.
#'
#' @param data A data.frame containing the dataset.
#' @param response Character. Name of the response variable.
#' @param factors Character vector. Names of factor variables.
#' @param design Character. Experimental design: "CRD" (default) or "RCBD".
#' @param alpha Numeric. Significance level for mean separation (default 0.05).
#' @param norm_test Character. Normality test: "shapiro" (default), "ad", or "ks".
#' @param homogeneity Character. Homogeneity test: "bartlett" (default), "levene", or "fligner".
#' @param mean_sep_test Character. Mean separation test: "LSD" (default), "HSD", or "SNK".
#'
#' @return A list containing:
#' \item{design}{Experimental design used.}
#' \item{anova}{ANOVA table.}
#' \item{CV}{Coefficient of variation (%).}
#' \item{mean_value}{Overall mean of the response.}
#' \item{normality_test}{Result of normality test.}
#' \item{homogeneity_test}{Result of homogeneity of variance test.}
#' \item{mean_separation}{List of mean separation results for main effects and interactions.}
#' \item{model}{Fitted aov model object.}
#'
#' @examples
#' # CRD with LSD mean separation
#' res <- FD(data = mydata, response = "yield", factors = c("A", "B"))
#' res$anova
#' res$mean_separation$A$groups
#'
#' # RCBD with HSD and Anderson-Darling normality test
#' res <- FD(data = mydata, response = "yield", factors = c("Block","A","B"),
#'           design = "RCBD", mean_sep_test = "HSD", norm_test = "ad")
#' @export
FD <- function(data, response, factors, design="CRD",
               alpha=0.05,
               norm_test="shapiro",
               homogeneity="bartlett",
               mean_sep_test="LSD") {

  library(agricolae)
  library(car)
  library(nortest)

  # Prepare formula
  fmla <- as.formula(
    paste(response, "~", paste(factors, collapse="*"))
  )

  # Fit ANOVA
  if(design=="CRD"){
    model <- aov(fmla, data=data)
  } else if(design=="RCBD"){
    block <- factors[1]
    fmla <- as.formula(
      paste(response, "~", block, "+", paste(factors[-1], collapse="*"))
    )
    model <- aov(fmla, data=data)
  } else {
    stop("Only CRD or RCBD supported")
  }

  # Compute residuals
  resid_vals <- residuals(model)

  # Normality test
  norm_result <- switch(norm_test,
                        shapiro = shapiro.test(resid_vals),
                        ad = ad.test(resid_vals),
                        ks = ks.test(resid_vals, "pnorm", mean(resid_vals), sd(resid_vals))
  )

  # Homogeneity of variance
  hom_result <- switch(homogeneity,
                       bartlett = bartlett.test(data[[response]], data[[factors[1]]]),
                       levene = leveneTest(fmla, data=data),
                       fligner = fligner.test(data[[response]], data[[factors[1]]])
  )

  # Coefficient of variation
  MSres <- summary(model)[[length(summary(model))]]["Residuals","Mean Sq"]
  CV <- sqrt(MSres)/mean(data[[response]])*100

  # Mean separation
  mean_sep <- list()
  all_terms <- attr(model$terms,"term.labels")

  for (trm in all_terms) {
    # If interaction, create temporary factor and refit model
    if (grepl(":", trm)) {
      facs <- unlist(strsplit(trm, ":"))
      int_name <- paste(facs, collapse=":")
      data[[int_name]] <- interaction(data[, facs], sep=":")
      tmp_model <- aov(data[[response]] ~ data[[int_name]])

      mean_sep[[trm]] <- switch(mean_sep_test,
                                LSD = LSD.test(tmp_model, "data[[int_name]]", alpha=alpha, group=TRUE),
                                HSD = HSD.test(tmp_model, "data[[int_name]]", alpha=alpha, group=TRUE),
                                SNK = SNK.test(tmp_model, "data[[int_name]]", alpha=alpha, group=TRUE)
      )

    } else {
      mean_sep[[trm]] <- switch(mean_sep_test,
                                LSD = LSD.test(model, trm, alpha=alpha, group=TRUE),
                                HSD = HSD.test(model, trm, alpha=alpha, group=TRUE),
                                SNK = SNK.test(model, trm, alpha=alpha, group=TRUE)
      )
    }
  }

  # Return results
  list(
    design = design,
    anova = summary(model),
    CV = CV,
    mean_value = mean(data[[response]]),
    normality_test = norm_result,
    homogeneity_test = hom_result,
    mean_separation = mean_sep,
    model = model
  )
}
