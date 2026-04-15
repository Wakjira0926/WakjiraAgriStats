#' Simple Linear Regression with Publication-Ready Plot
#'
#' @description
#' Fits a simple linear regression model and produces a publication-quality
#' ggplot showing observed data points, a solid regression line, regression
#' equation, and R-squared value. Users can customize colors, themes, text size,
#' and annotation placement. Optionally returns the model summary.
#'
#' @param data A data.frame containing the variables.
#' @param x Character. Name of the predictor variable.
#' @param y Character. Name of the response variable.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param main Character. Plot title.
#' @param point_color Color of observation points.
#' @param line_color Color of regression line.
#' @param eq_x Numeric. X position for regression equation.
#' @param eq_y Numeric. Y position for regression equation.
#' @param eq_size Numeric. Text size of regression equation.
#' @param axis_title_size Numeric. Axis title font size.
#' @param jitter_width Numeric. Jitter width for points.
#' @param theme_type Character. One of "bw", "classic", or "minimal".
#' @param return_model Logical. If TRUE, returns model summary.
#'
#' @return
#' A ggplot object. If return_model = TRUE, returns a list with:
#' \itemize{
#'   \item plot: ggplot object
#'   \item model: lm summary
#' }
#'
#' @examples
#' SL(ALdata, x = "N", y = "GY", eq_x = 2, eq_y = 8600)
#'
#' @export
SL <- function(data, x, y, xlab = NULL, ylab = NULL, main = NULL,
               point_color = "black", line_color = "blue",
               eq_x = NULL, eq_y = NULL, eq_size = 5, axis_title_size = 12,
               jitter_width = 0.2, theme_type = "bw", return_model = FALSE) {

  # Fit model
  formula <- as.formula(paste(y, "~", x))
  model <- lm(formula, data = data)
  summ <- summary(model)

  # Extract equation and R2
  eq <- paste0("y = ", round(coef(model)[1],2),
               ifelse(coef(model)[2]>=0," + "," - "),
               round(abs(coef(model)[2]),2),
               "*x\nR² = ", round(summ$r.squared,3))

  # Base ggplot
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(color = point_color, size = 2, position = position_jitter(width = jitter_width)) +
    geom_smooth(method = "lm", color = line_color, se = FALSE, linetype = "solid") +
    labs(x = ifelse(is.null(xlab), x, xlab),
         y = ifelse(is.null(ylab), y, ylab),
         title = main) +
    theme(axis.title = element_text(size = axis_title_size))

  # Add equation
  if(!is.null(eq_x) & !is.null(eq_y)) {
    p <- p + annotate("text", x = eq_x, y = eq_y, label = eq, hjust = 0, size = eq_size)
  } else {
    p <- p + annotate("text", x = max(data[[x]])*0.7, y = max(data[[y]])*0.9, label = eq, hjust = 0, size = eq_size)
  }

  # Theme selection
  if(theme_type == "bw") p <- p + theme_bw() + theme(panel.grid = element_blank())
  if(theme_type == "minimal") p <- p + theme_minimal() + theme(panel.grid = element_blank())
  if(theme_type == "classic") p <- p + theme_classic()

  # Return
  if(return_model) return(list(plot = p, model = summ))
  return(p)
}
