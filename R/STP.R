#' Strip-Plot ANOVA with Mean Separation, Normality, and Homogeneity Tests
#'
#' Performs a full strip-plot analysis of variance (ANOVA) on a dataset,
#' computes mean separation for factors, tests normality of residuals, and
#' evaluates homogeneity of variance. The function supports multiple options
#' for normality tests, homogeneity tests, and mean separation methods, but
#' only computes and displays the ones selected.
#'
#' @param data A \code{data.frame} containing the experimental data.
#' @param responses Character vector of column names for the response variable(s).
#' @param Replication Column name for the replication/block factor.
#' @param Column Column name for the column factor (main strip factor).
#' @param Row Column name for the row factor (main strip factor).
#' @param alpha Significance level for mean separation tests and confidence intervals (default 0.05).
#' @param normtest Character vector specifying which normality test(s) to run.
#'   Options: \code{"shapiro"} (Shapiro-Wilk), \code{"lillie"} (Lilliefors),
#'   \code{"adf"} (Anderson-Darling). Default: \code{"shapiro"}.
#' @param homogeneity Character vector specifying which homogeneity of variance test(s) to run.
#'   Options: \code{"bartlett"}, \code{"levene"}, \code{"fligner"}. Default: \code{"bartlett"}.
#' @param meanseptest Character vector specifying which mean separation test(s) to perform.
#'   Options: \code{"LSD"}, \code{"HSD"}, \code{"Duncan"}, \code{"SNK"}. Default: \code{"LSD"}.
#'
#' @return A list of results, one element per response variable. Each element is a list containing:
#' \describe{
#'   \item{ANOVA}{A data.frame of combined strip-plot ANOVA with terms: block, Column, Ea, Row, Eb, Interaction, Ec.}
#'   \item{CV}{Coefficients of variation (CV) for each stratum: CV(a), CV(b), CV(c).}
#'   \item{R2}{R-squared value for the model.}
#'   \item{Normality}{A list of selected normality test results (e.g., Shapiro-Wilk).}
#'   \item{Homogeneity}{A list of selected homogeneity test results (e.g., Bartlett test).}
#'   \item{MeanSeparation}{A list of mean separation results for Column, Row, and Interaction factors using selected methods.}
#' }
#'
#' @details
#' The function fits a strip-plot model:
#' \code{response ~ Replication + Column + Row + Column:Row + Error(Replication/Column + Replication/Row)}
#' and correctly calculates the block F-value using the within-stratum residual (Ec).
#' You can select which normality test, homogeneity test, and mean separation method(s) to run.
#' Only the selected tests are computed and displayed.
#' Mean separation uses the \code{agricolae} package and supports LSD, HSD, Duncan, and SNK methods.
#'
#' @examples
#' \dontrun{
#' # Single response with Shapiro normality, Bartlett homogeneity, LSD mean separation
#' res <- STP(SPData,
#'            responses="Yield",
#'            Replication="Rep",
#'            Column="Date_of_Sowing",
#'            Row="Varities",
#'            alpha=0.05,
#'            normtest="shapiro",
#'            homogeneity="bartlett",
#'            meanseptest="LSD")
#'
#' # Access combined ANOVA table
#' res$Yield$ANOVA
#'
#' # Access Shapiro-Wilk test
#' res$Yield$Normality$shapiro
#'
#' # Access LSD groups for Column factor
#' res$Yield$MeanSeparation$Date_of_Sowing$LSD$groups
#' }
#'
#' @importFrom agricolae LSD.test HSD.test duncan.test SNK.test
#' @importFrom car leveneTest
#' @importFrom nortest lillie.test adf.test
#' @export
STP <- function(data, responses, Replication, Column, Row,
                alpha = 0.05,
                normtest = c("shapiro","lillie","adf"),
                homogeneity = c("bartlett","levene","fligner"),
                meanseptest = c("LSD","HSD","Duncan","SNK")) {
  # Function body here (same as the final version we built)
}
###################
STP <- function(data, responses, Replication, Column, Row,
                alpha = 0.05,
                normtest = c("shapiro","lillie","adf"),
                homogeneity = c("bartlett","levene","fligner"),
                meanseptest = c("LSD","HSD","Duncan","SNK")) {

  # -----------------------------
  # Required packages
  # -----------------------------
  if(!requireNamespace("agricolae", quietly=TRUE)) stop("Install 'agricolae'")
  if(!requireNamespace("car", quietly=TRUE)) stop("Install 'car'")
  if(!requireNamespace("nortest", quietly=TRUE)) stop("Install 'nortest'")

  # -----------------------------
  # Convert columns to factors
  # -----------------------------
  data[[Replication]] <- as.factor(data[[Replication]])
  data[[Column]]      <- as.factor(data[[Column]])
  data[[Row]]         <- as.factor(data[[Row]])

  # -----------------------------
  # Helper to safely extract term
  # -----------------------------
  get_term <- function(tbl, pat, col){
    if(is.null(tbl)) return(NA_real_)
    rn <- rownames(tbl)
    id <- grep(pat, rn, ignore.case=TRUE)
    if(length(id)>0) return(as.numeric(tbl[id[1], col]))
    NA_real_
  }

  # -----------------------------
  # Mean separation helper
  # -----------------------------
  run_mean_sep <- function(fac_name, df_err, ms_err, factor_vec, alpha, method, response_name, data){
    if(length(unique(factor_vec))==1) return(NULL)
    res <- list()
    if("LSD" %in% method){
      grp <- agricolae::LSD.test(y=data[[response_name]], trt=factor_vec,
                                 DFerror=df_err, MSerror=ms_err,
                                 alpha=alpha, group=TRUE)
      res$LSD <- list(statistics=grp$statistics, groups=grp$groups)
    }
    if("HSD" %in% method){
      grp <- agricolae::HSD.test(y=data[[response_name]], trt=factor_vec,
                                 DFerror=df_err, MSerror=ms_err,
                                 alpha=alpha, group=TRUE)
      res$HSD <- list(statistics=grp$statistics, groups=grp$groups)
    }
    if("Duncan" %in% method){
      grp <- agricolae::duncan.test(y=data[[response_name]], trt=factor_vec,
                                    DFerror=df_err, MSerror=ms_err,
                                    alpha=alpha, group=TRUE)
      res$Duncan <- list(statistics=grp$statistics, groups=grp$groups)
    }
    if("SNK" %in% method){
      grp <- agricolae::SNK.test(y=data[[response_name]], trt=factor_vec,
                                 DFerror=df_err, MSerror=ms_err,
                                 alpha=alpha, group=TRUE)
      res$SNK <- list(statistics=grp$statistics, groups=grp$groups)
    }
    return(res)
  }

  # -----------------------------
  # Results list
  # -----------------------------
  results <- list()

  for(response in responses){
    cat("\nProcessing:", response, "\n")

    # --- Fit strip-plot model ---
    fmla <- as.formula(paste0(
      response," ~ ",Replication," + ",Column," + ",Row," + ",
      Column,":",Row," + Error(",Replication,"/",Column," + ",Replication,"/",Row,")"))
    model <- aov(fmla, data=data)
    sm <- summary(model)

    # --- Stratum tables ---
    tab_block  <- sm[[1]][[1]]
    tab_col    <- sm[[2]][[1]]
    tab_row    <- sm[[3]][[1]]
    tab_within <- sm[[4]][[1]]

    # --- Residuals ---
    MSEa <- get_term(tab_col,"Residual","Mean Sq")
    DFEa <- get_term(tab_col,"Residual","Df")
    MSEb <- get_term(tab_row,"Residual","Mean Sq")
    DFEb <- get_term(tab_row,"Residual","Df")
    MSEc <- get_term(tab_within,"Residual","Mean Sq")
    DFEc <- get_term(tab_within,"Residual","Df")

    # --- Block F using Ec ---
    MS_block <- get_term(tab_block,Replication,"Mean Sq")
    DF_block <- get_term(tab_block,Replication,"Df")
    F_block  <- MS_block / MSEc
    P_block  <- pf(F_block, DF_block, DFEc, lower.tail=FALSE)

    # --- Combined ANOVA table ---
    anova_tbl <- data.frame(
      Df=c(DF_block,
           get_term(tab_col,Column,"Df"),
           DFEa,
           get_term(tab_row,Row,"Df"),
           DFEb,
           get_term(tab_within,paste(Column,Row,sep=":"),"Df"),
           DFEc),
      `Sum Sq`=c(get_term(tab_block,Replication,"Sum Sq"),
                 get_term(tab_col,Column,"Sum Sq"),
                 get_term(tab_col,"Residual","Sum Sq"),
                 get_term(tab_row,Row,"Sum Sq"),
                 get_term(tab_row,"Residual","Sum Sq"),
                 get_term(tab_within,paste(Column,Row,sep=":"),"Sum Sq"),
                 get_term(tab_within,"Residual","Sum Sq")),
      `Mean Sq`=c(MS_block,
                  get_term(tab_col,Column,"Mean Sq"),
                  MSEa,
                  get_term(tab_row,Row,"Mean Sq"),
                  MSEb,
                  get_term(tab_within,paste(Column,Row,sep=":"),"Mean Sq"),
                  MSEc),
      `F value`=c(F_block,
                  get_term(tab_col,Column,"F value"),
                  NA,
                  get_term(tab_row,Row,"F value"),
                  NA,
                  get_term(tab_within,paste(Column,Row,sep=":"),"F value"),
                  NA),
      `Pr(>F)`=c(P_block,
                 get_term(tab_col,Column,"Pr(>F)"),
                 NA,
                 get_term(tab_row,Row,"Pr(>F)"),
                 NA,
                 get_term(tab_within,paste(Column,Row,sep=":"),"Pr(>F)"),
                 NA),
      row.names=c("block","Column","Ea","Row","Eb","Interaction","Ec"),
      check.names=FALSE
    )

    # Round for display
    anova_tbl$Df <- round(anova_tbl$Df,0)
    anova_tbl$`Sum Sq` <- round(anova_tbl$`Sum Sq`,4)
    anova_tbl$`Mean Sq` <- round(anova_tbl$`Mean Sq`,4)
    anova_tbl$`F value` <- ifelse(is.na(anova_tbl$`F value`), NA, round(anova_tbl$`F value`,7))
    anova_tbl$`Pr(>F)` <- ifelse(is.na(anova_tbl$`Pr(>F)`), NA, round(anova_tbl$`Pr(>F)`,7))

    print(anova_tbl, right=TRUE)

    # --- CV ---
    meanY <- mean(data[[response]], na.rm=TRUE)
    CVa <- sqrt(MSEa)/meanY*100
    CVb <- sqrt(MSEb)/meanY*100
    CVc <- sqrt(MSEc)/meanY*100
    cat("\nCV(a):",round(CVa,3),", CV(b):",round(CVb,3),", CV(c):",round(CVc,3),"\n")

    # --- R² ---
    totalSS <- sum(anova_tbl$`Sum Sq`, na.rm=TRUE)
    residSS <- anova_tbl["Ec","Sum Sq"]
    R2 <- 1 - residSS/totalSS
    cat("R Square:",round(R2,3),"\n")

    # --- Normality test ---
    norm_results <- NULL
    resid_within <- tryCatch(residuals(model$Within), error=function(e) NULL)
    if(!is.null(resid_within)){
      if("shapiro" %in% normtest){
        norm_results$shapiro <- shapiro.test(resid_within)
        print(norm_results$shapiro)
      }
      if("lillie" %in% normtest){
        norm_results$lillie <- nortest::lillie.test(resid_within)
        print(norm_results$lillie)
      }
      if("adf" %in% normtest){
        norm_results$adf <- nortest::adf.test(resid_within)
        print(norm_results$adf)
      }
    }

    # --- Homogeneity test ---
    hom_results <- NULL
    interfac <- interaction(data[[Column]], data[[Row]])
    if("bartlett" %in% homogeneity) hom_results$bartlett <- bartlett.test(data[[response]] ~ interfac)
    if("levene"  %in% homogeneity) hom_results$levene  <- car::leveneTest(data[[response]] ~ interfac)
    if("fligner" %in% homogeneity) hom_results$fligner <- fligner.test(data[[response]] ~ interfac)
    print(hom_results)

    # --- Mean separation ---
    lsd_results <- list()
    lsd_results[[Column]] <- run_mean_sep(Column, DFEa, MSEa, data[[Column]], alpha, meanseptest, response, data)
    lsd_results[[Row]]    <- run_mean_sep(Row, DFEb, MSEb, data[[Row]], alpha, meanseptest, response, data)
    intert <- interaction(data[[Column]], data[[Row]])
    lsd_results[["Interaction"]] <- run_mean_sep("Interaction", DFEc, MSEc, intert, alpha, meanseptest, response, data)

    # --- Store results ---
    results[[response]] <- list(
      ANOVA=anova_tbl,
      CV=c(CVa=CVa, CVb=CVb, CVc=CVc),
      R2=R2,
      Normality=norm_results,
      Homogeneity=hom_results,
      MeanSeparation=lsd_results
    )
  }

  invisible(results)
}

