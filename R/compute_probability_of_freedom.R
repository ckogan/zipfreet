source("R/is_valid_length.R")
source("R/prevpdf.R")

#' Compute probability of disease freedom given samples taken per timestep
#' 
#' @param n number of samples taken per timestep
#' @param phi_prior prior probability of freedom
#' @param alpha beta distribution alpha shape parameter
#' @param beta beta distribution beta shape parameter
#' @param p_intro probability of introduction
#' @param growth_rate exponential growth rate
#' @param rho unit test sensitivity
#' @param pi design threshold
#' @param delta_t test time duration
#' @param pi_seq discretization granularity
#' @returns "diseasefree" structure including sample sizes and prior/posterior distributions
#' @examples
#' u <- compute_probability_of_freedom(c(21,4,4,4,4), 0.5, 1, 1, 0.04, 0.01, 0.9, 0.0)
#' summary(u)
#' plot(u)
#' @export
compute_probability_of_freedom <- function(n, phi_prior, alpha, beta, p_intro, growth_rate, rho, pi=0, delta_t=1, pi_seq=1000)
{
  # phi_prior must be length 1
  if (length(phi_prior) > 1)
  {
    stop( simpleError("phi_prior must be length 1"))
  }  
  
  # the size of n drives the number of steps to take
  n_steps = length(n)
  
  # if any input vector has length > 1 but < n_steps, error out
  if (!is_valid_length(alpha, 1, n_steps)) stop(simpleError("Vector length mismatch: 'alpha'"))
  if (!is_valid_length(beta, 1, n_steps)) stop(simpleError("Vector length mismatch: 'beta'"))
  if (!is_valid_length(p_intro, 1, n_steps)) stop(simpleError("Vector length mismatch: 'p_intro'"))
  if (!is_valid_length(rho, 1, n_steps)) stop(simpleError("Vector length mismatch: 'rho'"))
  if (!is_valid_length(growth_rate, 1, n_steps)) stop(simpleError("Vector length mismatch: 'growth_rate'"))
  if (!is_valid_length(pi, 1, n_steps)) stop(simpleError("Vector length mismatch: 'pi'"))
  if (!is_valid_length(delta_t, 1, n_steps)) stop(simpleError("Vector length mismatch: 'delta_t'"))
  
  # make all vectors the same length (matching number of steps)
  if (length(alpha) < n_steps) alpha <- rep(alpha, n_steps)
  if (length(beta) < n_steps) beta <- rep(beta, n_steps)
  if (length(p_intro) < n_steps) p_intro <- rep(p_intro, n_steps)
  if (length(rho) < n_steps) rho <- rep(rho, n_steps)
  if (length(pi) < n_steps) pi <- rep(pi, n_steps)
  if (length(growth_rate) < n_steps) growth_rate <- rep(growth_rate, n_steps)
  if (length(delta_t) < n_steps) delta_t <- rep(delta_t, n_steps)
  
  # Create an instance of PrevPdf
  prevpdf <- PrevPdf$new(
    alpha = alpha[1],
    beta = beta[1],
    phi_prior = phi_prior,
    pi_seq = seq(0, 1, length.out=pi_seq)
  )
  
  p_eff_freedom_prior <- rep(NA, n_steps)
  p_eff_freedom_post  <- rep(NA, n_steps)
  for (j in 1:n_steps)
  {
    prevpdf$update(n[j], alpha[j], beta[j], rho[j], growth_rate[j], delta_t[j], 1 - p_intro[j])
    p_eff_freedom_prior[j] <- prevpdf$compute_posterior_cdf_given_n(pi[j], 0, rho[j])
    p_eff_freedom_post[j] <- prevpdf$compute_cdf(pi[j], step=j)
  }
  
  result <- list(n                   = n,
                 phi_prior           = prevpdf$phi_prior,
                 phi_post            = prevpdf$phi_post,
                 p_eff_freedom_prior = p_eff_freedom_prior,
                 p_eff_freedom_post  = p_eff_freedom_post,
                 f_posterior         = prevpdf$f_posterior_list)
  class(result) <- "diseasefree"
  return( result )  
}
