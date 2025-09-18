#' @exportS3Method
plot.zipfreet <- function(object, type="posterior", time=NULL, ...)
{
  if (length(time) == 0)
  {
    time = seq(1, length(object$n))
  }
  
  x <- seq(0, 1, 0.01)
  
  df <- data.frame(x=x)
  for (t in time)
  {
    colname <- paste0("time.", t)
    df[colname] <- sapply(x, object$f_posterior[[t]]$f)
  }
  df <- tidyr::gather(df, "time", "value", -x)
  
  plt <- ggplot2::ggplot(
                    data=df, 
                    mapping=ggplot2::aes(x,value,color=time)) +
            ggplot2::geom_line()
  
  print( plt )
}
