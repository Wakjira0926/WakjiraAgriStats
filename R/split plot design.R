SP <- function(data,
               trait,
               main,
               sub,
               block,
               var_test = "bartlett",
               norm_test = "shapiro",
               mean_method = "lsd",
               check_normality = TRUE,
               check_variance = TRUE,
               alpha = 0.05) {

  results <- list()

  for (tr in trait) {

    cat("\n================ Trait:", tr, "================\n")

    # Build model
    form <- as.formula(paste(tr, "~", block, "+", main, "*", sub, "+ Error(", block, "/", main, ")"))
    model <- aov(form, data = data)

    sum_model <- summary(model)

    block_tab  <- sum_model[[1]][[1]]
    main_tab   <- sum_model[[2]][[1]]
    within_tab <- sum_model[[3]][[1]]

    # Classical ANOVA table
    classic_table <- rbind(
      cbind(Source = "Block", block_tab[1,]),
      cbind(Source = "Main Plot", main_tab[main,]),
      cbind(Source = "Error A", main_tab["Residuals",]),
      cbind(Source = "Sub Plot", within_tab[sub,]),
      cbind(Source = "Main × Sub", within_tab[paste(main, sub, sep=":"),]),
      cbind(Source = "Error B", within_tab["Residuals",])
    )

    classic_table <- as.data.frame(classic_table)

    cat("\n--- Classical Split-Plot ANOVA ---\n")
    print(classic_table)

    # Extract residuals
    res_main <- residuals(model$`Error: bloco:parc`)
    res_sub  <- residuals(model$`Error: Within`)

    # Assumption tests
    assumptions <- list()

    if (check_normality) {
      assumptions$normality_main <- tryCatch({
        switch(norm_test,
               shapiro = shapiro.test(res_main),
               ks = ks.test(res_main, "pnorm"),
               ad = nortest::ad.test(res_main))
      }, error = function(e) "Too few residuals")

      assumptions$normality_sub <- tryCatch({
        switch(norm_test,
               shapiro = shapiro.test(res_sub),
               ks = ks.test(res_sub, "pnorm"),
               ad = nortest::ad.test(res_sub))
      }, error = function(e) "Too few residuals")
    }

    if (check_variance) {
      form_var <- as.formula(paste(tr, "~", main))

      assumptions$variance <- tryCatch({
        switch(var_test,
               bartlett = bartlett.test(form_var, data = data),
               levene = car::leveneTest(form_var, data = data))
      }, error = function(e) "Variance test failed")
    }

    cat("\n--- Assumption Checks ---\n")
    print(assumptions)

    # CV calculations
    cv_main <- sqrt(main_tab["Residuals", "Mean Sq"]) / mean(data[[tr]], na.rm = TRUE) * 100
    cv_sub  <- sqrt(within_tab["Residuals", "Mean Sq"]) / mean(data[[tr]], na.rm = TRUE) * 100

    cat("\n--- Coefficient of Variation ---\n")
    cat("Main Plot CV:", round(cv_main, 2), "%\n")
    cat("Sub Plot CV:", round(cv_sub, 2), "%\n")

    # Mean separation
    mean_sep_main <- tryCatch({
      mean_separation(model, effect = main, method = mean_method, alpha = alpha)
    }, error = function(e) NULL)

    mean_sep_sub <- tryCatch({
      mean_separation(model, effect = sub, method = mean_method, alpha = alpha)
    }, error = function(e) NULL)

    mean_sep_inter <- tryCatch({
      mean_separation(model, effect = paste(main, sub, sep=":"), method = mean_method, alpha = alpha)
    }, error = function(e) NULL)

    cat("\n--- Mean Separation ---\n")
    cat("\nMain Plot Means:\n")
    print(mean_sep_main)
    cat("\nSub Plot Means:\n")
    print(mean_sep_sub)
    cat("\nInteraction Means:\n")
    print(mean_sep_inter)

    # Save results
    results[[tr]] <- list(
      model = model,
      anova = classic_table,
      assumptions = assumptions,
      cv_main = cv_main,
      cv_sub = cv_sub,
      mean_sep_main = mean_sep_main,
      mean_sep_sub = mean_sep_sub,
      mean_sep_inter = mean_sep_inter
    )
  }

  invisible(results)
}
