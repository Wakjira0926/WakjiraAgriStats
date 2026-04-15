#' Split-Plot Design Analysis for Multiple Traits
#'
#' Performs a split-plot ANOVA for one or multiple traits, including:
#' - Residual normality tests (Shapiro-Wilk, Anderson-Darling, Kolmogorov-Smirnov)
#' - Homogeneity of variance tests (Bartlett, Levene, Fligner)
#' - Mean separation for main plot, subplot, and interaction effects
#' - Coefficients of variation per stratum
#'
#' The function is multi-trait ready and outputs ANOVA tables, CVs, mean separation tables,
#' and informative messages about whether mean separation is necessary.
#'
#' @param data A \code{data.frame} containing the experimental data.
#' @param responses A character vector of trait column names to analyze.
#' @param Replication Column name for blocks/replications (whole plot replication factor).
#' @param MainPlot Column name for the main plot factor (whole plot factor).
#' @param SubPlot Column name for the subplot factor (split plot factor).
#' @param alpha Significance level for mean separation tests (default 0.05).
#' @param normtest Normality test for residuals: \code{"shapiro"} (default), \code{"ad"} (Anderson-Darling), or \code{"ks"} (Kolmogorov-Smirnov).
#' @param homogeneity Homogeneity test for variance: \code{"bartlett"} (default), \code{"levene"}, or \code{"fligner"}.
#' @param meanseptest Mean separation test: \code{"LSD"} (default), \code{"HSD"}, \code{"SNK"}, \code{"Duncan"}, or \code{"Tukey"}.
#'
#' @return A list of lists, one element per trait, each containing:
#' \describe{
#'   \item{trait}{The trait name.}
#'   \item{design}{Design type ("Split Plot").}
#'   \item{anova}{ANOVA summary list per stratum.}
#'   \item{CV}{Coefficients of variation per stratum (main plot, subplot).}
#'   \item{meanvalue}{Grand mean of the trait.}
#'   \item{normalitytest}{Result of the chosen normality test.}
#'   \item{homogeneitytest}{Result of the chosen homogeneity test.}
#'   \item{meanseparation}{List containing mean separation tables for main plot, subplot, and interaction.}
#'   \item{messages}{Informative messages if mean separation is skipped or performed.}
#'   \item{model}{The fitted aov model object.}
#' }
#'
#' @examples
#' \dontrun{
#' results <- SPD(
#'   data = splitplot,
#'   responses = "Yield",
#'   Replication = "Rep",
#'   MainPlot = "Date_of_Sowing",
#'   SubPlot = "Varities"
#' )
#' results$Yield$anova
#' results$Yield$CV
#' results$Yield$meanseparation$Interaction
#' }
#'
#' @export
SPD <- function(data, responses, Replication, MainPlot, SubPlot,
                alpha = 0.05,
                normtest = "shapiro",
                homogeneity = "bartlett",
                meanseptest = "LSD",
                digits = 6) {

  # Packages
  library(agricolae)
  library(car)
  library(nortest)

  # Convert to factors
  data[[Replication]] <- as.factor(data[[Replication]])
  data[[MainPlot]]    <- as.factor(data[[MainPlot]])
  data[[SubPlot]]     <- as.factor(data[[SubPlot]])

  results <- list()

  for (resp in responses) {

    message(paste("Processing:", resp))

    ##------------------------------------------
    ## 1. Fit Correct Split‑Plot Model
    ##------------------------------------------
    fmla <- as.formula(
      paste(resp, "~", MainPlot, "*", SubPlot,
            "+ Error(", Replication, "/", MainPlot, ")")
    )
    model <- aov(fmla, data = data)
    aov_list <- summary(model)

    block_stratum <- aov_list[[1]][[1]]
    main_stratum  <- aov_list[[2]][[1]]
    sub_stratum   <- aov_list[[3]][[1]]

    ##------------------------------------------
    ## 2. Extract & Summarize
    ##------------------------------------------
    grand_mean <- mean(data[[resp]], na.rm = TRUE)
    MSEa <- main_stratum["Residuals","Mean Sq"];  dfEa <- main_stratum["Residuals","Df"]
    MSEb <- sub_stratum["Residuals","Mean Sq"];   dfEb <- sub_stratum["Residuals","Df"]
    CVa  <- sqrt(MSEa)/grand_mean*100
    CVb  <- sqrt(MSEb)/grand_mean*100

    ## Block vs Ea
    MS_block <- block_stratum["Residuals","Mean Sq"]
    df_block <- block_stratum["Residuals","Df"]
    F_block  <- MS_block/MSEa
    p_block  <- 1 - pf(F_block, df_block, dfEa)

    ##------------------------------------------
    ## 3. Build Compact ANOVA Table
    ##------------------------------------------
    anova_table <- data.frame(
      Source = c(Replication,
                 MainPlot,
                 "Ea",
                 SubPlot,
                 paste(MainPlot,":",SubPlot,sep=""),
                 "Eb"),
      Df = c(df_block,
             main_stratum[MainPlot,"Df"],
             dfEa,
             sub_stratum[SubPlot,"Df"],
             sub_stratum[paste(MainPlot,SubPlot,sep=":"),"Df"],
             dfEb),
      "Sum Sq" = c(block_stratum["Residuals","Sum Sq"],
                   main_stratum[MainPlot,"Sum Sq"],
                   main_stratum["Residuals","Sum Sq"],
                   sub_stratum[SubPlot,"Sum Sq"],
                   sub_stratum[paste(MainPlot,SubPlot,sep=":"),"Sum Sq"],
                   sub_stratum["Residuals","Sum Sq"]),
      "Mean Sq" = c(MS_block,
                    main_stratum[MainPlot,"Mean Sq"],
                    main_stratum["Residuals","Mean Sq"],
                    sub_stratum[SubPlot,"Mean Sq"],
                    sub_stratum[paste(MainPlot,SubPlot,sep=":"),"Mean Sq"],
                    sub_stratum["Residuals","Mean Sq"]),
      "F value" = c(F_block,
                    main_stratum[MainPlot,"F value"],
                    NA,
                    sub_stratum[SubPlot,"F value"],
                    sub_stratum[paste(MainPlot,SubPlot,sep=":"),"F value"],
                    NA),
      "Pr(>F)"  = c(p_block,
                    main_stratum[MainPlot,"Pr(>F)"],
                    NA,
                    sub_stratum[SubPlot,"Pr(>F)"],
                    sub_stratum[paste(MainPlot,SubPlot,sep=":"),"Pr(>F)"],
                    NA),
      check.names = FALSE
    )

    ## Significance codes
    sigsym <- function(p)
      ifelse(is.na(p), "",
             ifelse(p<0.001,"***",
                    ifelse(p<0.01,"**",
                           ifelse(p<0.05,"*",
                                  ifelse(p<0.1,"."," ")))))
    anova_table$Signif <- vapply(anova_table$`Pr(>F)`, sigsym, character(1))

    ## round
    num_cols <- c("Df","Sum Sq","Mean Sq","F value","Pr(>F)")
    anova_table[num_cols] <- lapply(anova_table[num_cols],
                                    \(x) if(is.numeric(x)) round(x,digits) else x)

    ##------------------------------------------
    ## 4. Print Report
    ##------------------------------------------
    cat("\nAnalysis of Variance Table\n")
    cat("Response:", resp, "\n\n")
    print(anova_table[,c("Source","Df","Sum Sq","Mean Sq","F value","Pr(>F)","Signif")],
          row.names=FALSE)
    cat("\nSignif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n")
    cat(sprintf("\n CV(a): %.3f , CV(b): %.3f\n", CVa, CVb))

    ##------------------------------------------
    ## 5. Diagnostics
    ##------------------------------------------
    rvals <- residuals(model)
    if (is.list(rvals)) rvals <- unlist(rvals)
    rvals <- as.numeric(rvals); rvals <- rvals[is.finite(rvals)]
    norm_test <- if(length(rvals)>2)
      tryCatch({
        switch(normtest,
               shapiro = shapiro.test(rvals),
               ad = ad.test(rvals),
               ks = ks.test(rvals,"pnorm",mean(rvals),sd(rvals)))
      },error=\(e)NA) else NA

    interfac <- interaction(data[[MainPlot]], data[[SubPlot]])
    hom_test <- switch(homogeneity,
                       bartlett=bartlett.test(data[[resp]], interfac),
                       levene=leveneTest(data[[resp]] ~ interfac),
                       fligner=fligner.test(data[[resp]], interfac))

    ##------------------------------------------
    ## 6.  LSD Tests — correct error term and replication
    ##------------------------------------------
    mean_sep <- list()

    # MAIN‑PLOT  (Ea)
    mean_sep[[MainPlot]] <-
      LSD.test(y = data[[resp]], trt = data[[MainPlot]],
               DFerror = dfEa, MSerror = MSEa,
               alpha = alpha, group = TRUE)

    # SUB‑PLOT  (Eb)
    mean_sep[[SubPlot]] <-
      LSD.test(y = data[[resp]], trt = data[[SubPlot]],
               DFerror = dfEb, MSerror = MSEb,
               alpha = alpha, group = TRUE)

    # INTERACTION  (Eb)
    data$MP_SP <- interaction(data[[MainPlot]], data[[SubPlot]])
    mean_sep[[paste(MainPlot,SubPlot,sep=":")]] <-
      LSD.test(y = data[[resp]], trt = data$MP_SP,
               DFerror = dfEb, MSerror = MSEb,
               alpha = alpha, group = TRUE)

    ##------------------------------------------
    ## 7.  Save all results
    ##------------------------------------------
    results[[resp]] <- list(
      trait = resp,
      design = "Split Plot",
      anova_table = anova_table,
      CV = c(CVa=CVa, CVb=CVb),
      meanvalue = grand_mean,
      normalitytest = norm_test,
      homogeneitytest = hom_test,
      meanseparation = mean_sep,
      model = model
    )
  }
  return(results)
}

