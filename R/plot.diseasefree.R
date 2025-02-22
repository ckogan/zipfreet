library(ggplot2)
library(tidyr)

plot.diseasefree <- function(object, type="posterior", time=NULL, ...)
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
    df[colname] <- sapply(x, object$f_posterior[[t]])
  }
  df <- gather(df, "time", "value", -x)
  
  plt <- ggplot(data=df, mapping=aes(x,value,color=time)) + geom_line()
  
  print( plt )
}
