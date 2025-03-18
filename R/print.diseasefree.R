#' @exportS3Method
print.diseasefree <- function(object, ...)
{
  nrow <- length(object$n)
  cat("time     n p_freedom_prior p_freedom_post p(R+ | D+)\n")
  cat("---- ----- --------------- -------------- ----------\n")
  for (j in 1:nrow)
  {
    cat(sprintf("%4d %5d %15.3f %14.3f %10.3f\n", 
      j, object$n[j], object$p_freedom_prior[j], 
      object$p_freedom_post[j], object$sensitivity[j]))
  }
}
