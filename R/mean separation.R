# mean_separation.R
mean_separation <- function(model, method = "lsd", alpha = 0.05) {
  if(!requireNamespace("agricolae", quietly = TRUE)) {
    stop("Please install the 'agricolae' package first.")
  }

  tr <- names(model$model)[2]  # treatment column
  y <- model$model[[1]]        # response

  ms <- switch(tolower(method),
               lsd    = agricolae::LSD.test(model, tr, alpha = alpha),
               tukey  = agricolae::HSD.test(model, tr, alpha = alpha),
               duncan = agricolae::duncan.test(model, tr, alpha = alpha),
               snk    = agricolae::SNK.test(model, tr, alpha = alpha),
               hsd    = agricolae::HSD.test(model, tr, alpha = alpha),
               stop("Invalid mean separation method"))

  return(ms$groups)
}
