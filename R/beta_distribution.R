#' Compute alternative parameterizations of the beta distribution 
#' 
#' @param alpha alpha parameter
#' @param beta beta parameter
#' @param mu mean
#' @param sigma standard deviation of beta distribution
#' @export
beta_distribution_statistics <- function(alpha = NULL, beta = NULL, mu = NULL, sigma = NULL, alpha_range = c(0.0001, 10000)) {
  if(!(
    (!is.null(alpha) & is.null(beta) & !is.null(mu) & is.null(sigma)) |
    (!is.null(alpha) & !is.null(beta) & is.null(mu) & is.null(sigma)) |
    (is.null(alpha) & is.null(beta) & !is.null(mu) & !is.null(sigma)) 
    )) stop("Exactly two of alpha, beta, mu, sigma should be not NULL (we've considered 3 cases)")

  beta_f <- function(alpha, mu) alpha * (1-mu) / mu
  sigma2_f <- function(alpha, beta) (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
  if(!is.null(mu) & !is.null(sigma)) {
      alpha <- uniroot(function(alpha) sigma2_f(alpha, beta_f(alpha, mu)) - sigma^2, alpha_range)$root
  }
  
  if(is.null(beta) & !is.null(alpha) & !is.null(mu)) {
    beta <- beta_f(alpha, mu)
  } 
  
  sigma <- sqrt(sigma2_f(alpha, beta))
  c(alpha = alpha, beta = beta, mu = mu, sigma = sigma)
}