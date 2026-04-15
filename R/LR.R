#' Logistic Regression with Publication-Ready Plot
#'
#' @description
#' Fits a binary logistic regression model and optionally produces a
#' publication-quality ggplot with observed data and fitted probability curve.
#'
#' @param data A data.frame containing variables.
#' @param x Character vector. Predictor variable(s).
#' @param y Character. Binary response variable (0/1 or factor).
#' @param plot Logical. If TRUE and length(x)==1, produce a plot.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param main Character. Plot title.
#' @param point_color Color of observed points.
#' @param line_color Color of fitted logistic curve.
#' @param eq_x Numeric. X position for equation annotation.
#' @param eq_y Numeric. Y position for equation annotation.
#' @param eq_size Numeric. Size of equation text.
#' @param axis_title_size Numeric. Axis title size.
#' @param jitter_width Numeric. Jitter width for points.
#' @param theme_type Character. "bw", "classic", or "minimal".
#' @param return_model Logical. If TRUE, returns model object.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item Model: glm object
#'   \item Summary: coefficient table
#'   \item OddsRatio: odds ratios with CI
#'   \item Plot: ggplot object (if plot = TRUE)
#' }
#'
#' @examples
#' LR(ALdata, x = "N", y = "Response", plot = TRUE)
#'
#' @export

LR <- function(data, x, y,
               plot = TRUE,
               xlab = NULL, ylab = NULL, main = NULL,
               point_color = "black", line_color = "blue",
               eq_x = NULL, eq_y = NULL, eq_size = 5,
               axis_title_size = 12, jitter_width = 0.1,
               theme_type = "bw", return_model = FALSE) {

  library(ggplot2)

  # Formula
  form <- as.formula(paste(y, "~", paste(x, collapse = " + ")))

  # Fit model
  model <- glm(form, data = data, family = binomial)

  summ <- summary(model)

  # Coefficient table
  coef_tab <- data.frame(
    Estimate = summ$coefficients[,1],
    Std.Error = summ$coefficients[,2],
    z.value = summ$coefficients[,3],
    Pr = summ$coefficients[,4]
  )

  # Odds ratios
  OR <- exp(cbind(
    OR = coef(model),
    confint.default(model)
  ))

  # Pseudo R2 (McFadden)
  r2 <- 1 - model$deviance / model$null.deviance

  # ---- Plot only if single predictor ----
  p <- NULL
  if (plot && length(x) == 1) {

    data$.prob <- predict(model, type = "response")

    eq_text <- paste0(
      "logit(p) = ",
      round(coef(model)[1], 3),
      " + ",
      round(coef(model)[2], 3), "x\n",
      "McFadden R² = ", round(r2, 3)
    )

    p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
      geom_jitter(height = jitter_width, color = point_color, size = 2) +
      geom_smooth(method = "glm",
                  method.args = list(family = "binomial"),
                  se = FALSE, color = line_color, linetype = "solid") +
      labs(
        x = ifelse(is.null(xlab), x, xlab),
        y = ifelse(is.null(ylab), y, ylab),
        title = main
      ) +
      theme(axis.title = element_text(size = axis_title_size))

    if (!is.null(eq_x) && !is.null(eq_y)) {
      p <- p + annotate("text", x = eq_x, y = eq_y,
                        label = eq_text, size = eq_size, hjust = 0)
    }

    if (theme_type == "bw") p <- p + theme_bw() + theme(panel.grid = element_blank())
    if (theme_type == "classic") p <- p + theme_classic()
    if (theme_type == "minimal") p <- p + theme_minimal() + theme(panel.grid = element_blank())
  }

  out <- list(
    Model = model,
    Summary = coef_tab,
    OddsRatio = OR,
    R2 = r2,
    Plot = p
  )

  if (return_model) return(out)
  return(out)
}
