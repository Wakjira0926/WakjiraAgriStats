mean_separation <- function(model, method = "lsd", alpha = 0.05) {

  trt <- names(model$model)[2]

  switch(method,
         lsd     = agricolae::LSD.test(model, trt, alpha = alpha),
         tukey   = agricolae::HSD.test(model, trt, alpha = alpha),
         hsd     = agricolae::HSD.test(model, trt, alpha = alpha),
         duncan  = agricolae::duncan.test(model, trt, alpha = alpha),
         snk     = agricolae::SNK.test(model, trt, alpha = alpha),
         "Invalid mean separation method"
  )
}
