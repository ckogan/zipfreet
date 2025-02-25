#' @exportS3Method
summary.diseasefree <- function(object, ...)
{
  nrow <- length(object$n)
  cat("time phi_prior     n phi_post p_eff_freedom\n")
  cat("---- --------- ----- -------- -------------\n")
  for (j in 1:nrow)
  {
    cat(sprintf("%4d %9.3f %5d %8.3f %13.3f\n", 
      j, object$phi_prior[j], object$n[j], object$phi_post[j], object$p_eff_freedom[j]))
  }
}
