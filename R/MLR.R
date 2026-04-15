#' Multiple Linear Regression with Publication-Ready Plot
#'
#' @description
#' \code{ML()} fits a multiple linear regression model and produces a
#' publication-quality ggplot showing observed values, a solid regression
#' line, regression equation, and R-squared value. The function allows
#' extensive customization of colors, themes, text size, and annotation
#' placement, and can optionally return the fitted model.
#'
#' @usage
#' ML(
#'   data,
#'   x,
#'   y,
#'   xlab = NULL,
#'   ylab = NULL,
#'   main = NULL,
#'   point_color = "black",
#'   line_color = "blue",
#'   eq_x = NULL,
#'   eq_y = NULL,
#'   eq_size = 5,
#'   axis_title_size = 12,
#'   jitter_width = 0.2,
#'   theme_type = "bw",
#'   return_model = FALSE
#' )
#'
#' @arguments
#' \describe{
#'   \item{data}{A data frame containing response and predictor variables.}
#'   \item{x}{Character vector of predictor variable names.}
#'   \item{y}{Character name of the response variable.}
#'   \item{xlab}{X-axis label. If NULL, default is used.}
#'   \item{ylab}{Y-axis label. If NULL, default is used.}
#'   \item{main}{Main title of the plot.}
#'   \item{point_color}{Color of observed data points.}
#'   \item{line_color}{Color of the regression line.}
#'   \item{eq_x}{X-coordinate for regression equation placement.}
#'   \item{eq_y}{Y-coordinate for regression equation placement.}
#'   \item{eq_size}{Font size of regression equation text.}
#'   \item{axis_title_size}{Font size of axis titles.}
#'   \item{jitter_width}{Amount of horizontal jitter for points.}
#'   \item{theme_type}{Plot theme: one of "bw", "classic", "minimal", "void".}
#'   \item{return_model}{Logical; if TRUE, returns the fitted lm model.}
#' }
#'
#' @details
#' The regression is fitted using \code{lm()}.
#' For visualization, the first predictor in \code{x} is used on the x-axis.
#' The plot contains only x and y axes (no rectangular frame or grid lines),
#' making it suitable for publication.
#'
#' @value
#' Invisibly returns a ggplot object.
#' If \code{return_model = TRUE}, an object of class \code{lm} is returned.
#'
#' @examples
#' # Simple multiple regression
#' ML(
#'   data = ALdata,
#'   x = c("N", "P", "K"),
#'   y = "GY"
#' )
#'
#' # Customized plot
#' ML(
#'   data = ALdata,
#'   x = c("N", "P"),
#'   y = "GY",
#'   xlab = "Nitrogen (kg ha-1)",
#'   ylab = "Grain Yield (kg ha-1)",
#'   main = "Multiple Linear Regression of Yield",
#'   point_color = "darkgreen",
#'   line_color = "red",
#'   eq_x = 2,
#'   eq_y = 8600,
#'   eq_size = 5,
#'   axis_title_size = 14,
#'   theme_type = "classic"
#' )
#'
#' # Extract model results
#' fit <- ML(
#'   data = ALdata,
#'   x = c("N", "P", "K"),
#'   y = "GY",
#'   return_model = TRUE
#' )
#' summary(fit)
#'
#' @seealso
#' \code{\link{lm}}, \code{\link{ggplot}}
#'
#' @export

ML <- function(data, x, y,
               xlab = "Predicted",
               ylab = "Observed",
               main = NULL,
               point_color = "black",
               line_color = "red",
               eq_x = NULL,
               eq_y = NULL,
               eq_size = 5,
               axis_title_size = 12,
               axis_text_size = 10,
               axis_title_face = "bold",
               jitter_width = 0.2,
               theme_type = "bw",
               add_layers = NULL,
               return_model = FALSE) {

  library(ggplot2)

  # -----------------------------
  # Model
  # -----------------------------
  formula_txt <- paste(y, "~", paste(x, collapse = " + "))
  model <- lm(as.formula(formula_txt), data = data)

  data$.pred <- predict(model, newdata = data)

  # -----------------------------
  # Equation
  # -----------------------------
  coefs <- coef(model)
  coef_names <- names(coefs)

  eq_terms <- c(
    sprintf("%.3f", coefs[1]),
    sapply(2:length(coefs), function(i) {
      sign <- ifelse(coefs[i] >= 0, " + ", " − ")
      paste0(sign, sprintf("%.3f", abs(coefs[i])), "·", coef_names[i])
    })
  )

  eq_text <- paste0(
    "y = ", paste(eq_terms, collapse = ""),
    "\nR² = ", round(summary(model)$r.squared, 3)
  )

  # -----------------------------
  # Base plot
  # -----------------------------
  p <- ggplot(data, aes(x = .pred, y = .data[[y]])) +
    geom_jitter(
      width = jitter_width,
      color = point_color,
      size = 2
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = line_color,
      linewidth = 1
    ) +
    labs(
      x = xlab,
      y = ylab,
      title = main
    )

  # -----------------------------
  # Theme control
  # -----------------------------
  if (theme_type == "bw") {
    p <- p + theme_bw()
  } else if (theme_type == "classic") {
    p <- p + theme_classic()
  } else if (theme_type == "minimal") {
    p <- p + theme_minimal()
  }

  p <- p +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title = element_text(
        size = axis_title_size,
        face = axis_title_face
      ),
      axis.text = element_text(size = axis_text_size)
    )

  # -----------------------------
  # Equation annotation
  # -----------------------------
  if (!is.null(eq_x) && !is.null(eq_y)) {
    p <- p + annotate(
      "text",
      x = eq_x,
      y = eq_y,
      label = eq_text,
      size = eq_size,
      hjust = 0
    )
  }

  # -----------------------------
  # User-added ggplot2 layers
  # -----------------------------
  if (!is.null(add_layers)) {
    p <- p + add_layers
  }

  print(p)

  if (return_model) {
    return(summary(model))
  }
}
