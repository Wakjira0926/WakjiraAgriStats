FD <- function(data, responses, factors, design="CRD",
               alpha=0.05,
               norm_test="shapiro",
               homogeneity="bartlett",
               mean_sep_test="LSD") {

  # Load required packages
  require(agricolae)
  require(car)
  require(nortest)

  # Allow single response as character string
  if (length(responses) == 1) responses <- as.character(responses)

  # Function for a single response
  FD_single <- function(response) {

    # Build formula
    fmla <- if(design=="CRD") {
      as.formula(paste(response, "~", paste(factors, collapse="*")))
    } else if(design=="RCBD") {
      block <- factors[1]
      as.formula(paste(response, "~", block, "+", paste(factors[-1], collapse="*")))
    } else stop("Only CRD or RCBD supported")

    # Fit ANOVA
    model <- aov(fmla, data=data)
    resid_vals <- residuals(model)

    # Normality test
    norm_result <- switch(norm_test,
                          shapiro = shapiro.test(resid_vals),
                          ad = ad.test(resid_vals),
                          ks = ks.test(resid_vals, "pnorm", mean(resid_vals), sd(resid_vals))
    )

    # Homogeneity test
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

    for(trm in all_terms){
      # Interaction term handling
      if(grepl(":", trm)){
        facs <- unlist(strsplit(trm, ":"))
        int_name <- paste(facs, collapse=":")
        data[[int_name]] <- interaction(data[,facs], sep=":")

        # Temporary model for the interaction
        tmp_model <- aov(data[[response]] ~ data[[int_name]])

        mean_sep[[trm]] <- switch(mean_sep_test,
                                  LSD = LSD.test(tmp_model, "data[[int_name]]", alpha=alpha, group=TRUE),
                                  HSD = HSD.test(tmp_model, "data[[int_name]]", alpha=alpha, group=TRUE),
                                  SNK = SNK.test(tmp_model, "data[[int_name]]", alpha=alpha, group=TRUE)
        )
      } else {
        # Main effect
        mean_sep[[trm]] <- switch(mean_sep_test,
                                  LSD = LSD.test(model, trm, alpha=alpha, group=TRUE),
                                  HSD = HSD.test(model, trm, alpha=alpha, group=TRUE),
                                  SNK = SNK.test(model, trm, alpha=alpha, group=TRUE)
        )
      }
    }

    # Return single response result
    list(
      response = response,
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

  # Apply to all responses
  res_list <- lapply(responses, FD_single)
  names(res_list) <- responses
  return(res_list)
}
