is_valid_length <- function(x, min_size, max_size)
{
  return( length(x) == min_size || length(x) == max_size )
}
