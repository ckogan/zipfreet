source("R/is_valid_length.R")
source("R/prevpdf.R")

#' Compute the sample size requirements to meet a specified design
#' threshold at a specified confidence.
#' 
#' @param phi_prior prior probability of freedom
#' @param alpha_prior beta distribution alpha parameter for prior prevalence distribution
#' @param beta_prior beta distribution beta parameter for prior prevalence distribution
#' @param alpha_intro beta distribution alpha shape parameter for introduction distribution
#' @param beta_intro beta distribution beta shape parameter for introduction distribution
#' @param p_intro probability of introduction; if provided as a vector, the
#'                initial entry corresponds to the probability of introduction
#'                at time step 2.
#' @param rho unit test sensitivity
#' @param pi design threshold
#' @param dconf desired confidence
#' @param growth_rate exponential growth rate
#' @param delta_t test time duration
#' @param n_steps number of time units to step
#' @param method  "restore" or "maintain"
#' @param pi_seq discretization granularity
#' @param n_max maximum sample size
#' @returns "zipfreet" structure including sample sizes and prior/posterior distributions
#' @examples
#' u <- compute_sample_size(0.5, 1, 1, 1, 1, 0.04, 0.01, 0.9, 0, 0.95, n_steps=5)
#' print(u)
#' plot(u)
#' @export
compute_sample_size <- function(phi_prior, alpha_prior, beta_prior, alpha_intro, beta_intro, p_intro, growth_rate, rho, pi, dconf, delta_t=1, n_steps=1, method="restore", pi_seq=1000, n_max = 1000)
{
  # phi_prior must be length 1
  if (length(phi_prior) > 1)
  {
    stop( simpleError("phi_prior must be length 1"))
  }
  
  # alpha_prior, beta_prior should be length 1
  if (length(alpha_prior) > 1)
  {
    stop( simpleError("alpha_prior must be length 1"))
  }
  
  if (length(beta_prior) > 1)
  {
    stop( simpleError("beta_prior must be length 1"))
  }  

  # make sure the method is either "restore" or "maintain"
  if (method != "restore" && method != "maintain")
  {
    stop( simpleError("Invalid method; must be either 'restore' or 'maintain'"))
  }
  
  # determine max input vector length
  max_len = max( c(length(alpha_intro),
                   length(beta_intro),
                   length(p_intro),
                   length(rho),
                   length(dconf),
                   length(pi),
                   length(growth_rate),
                   length(delta_t)) )
  
  # the max vector size drives the number of steps to take
  if (max_len > 1)
  {
    n_steps = max_len
  }
  
  # if any input vector has length > 1 but < n_steps, error out
  if (!is_valid_length(alpha_intro, 1, n_steps)) stop(simpleError("Vector length mismatch: 'alpha_intro'"))
  if (!is_valid_length(beta_intro, 1, n_steps)) stop(simpleError("Vector length mismatch: 'beta_intro'"))
  if (!is_valid_length(p_intro, 1, n_steps)) stop(simpleError("Vector length mismatch: 'p_intro'"))
  if (!is_valid_length(rho, 1, n_steps)) stop(simpleError("Vector length mismatch: 'rho'"))
  if (!is_valid_length(dconf, 1, n_steps)) stop(simpleError("Vector length mismatch: 'dconf'"))
  if (!is_valid_length(pi, 1, n_steps)) stop(simpleError("Vector length mismatch: 'pi'"))
  if (!is_valid_length(growth_rate, 1, n_steps)) stop(simpleError("Vector length mismatch: 'growth_rate'"))
  if (!is_valid_length(delta_t, 1, n_steps)) stop(simpleError("Vector length mismatch: 'delta_t'"))
  
  # make all vectors the same length (matching number of steps)
  if (length(alpha_intro) < n_steps) alpha_intro <- rep(alpha_intro, n_steps)
  if (length(beta_intro) < n_steps) beta_intro <- rep(beta_intro, n_steps)
  if (length(p_intro) < n_steps) p_intro <- rep(p_intro, n_steps)
  if (length(rho) < n_steps) rho <- rep(rho, n_steps)
  if (length(dconf) < n_steps) dconf <- rep(dconf, n_steps)
  if (length(pi) < n_steps) pi <- rep(pi, n_steps)
  if (length(growth_rate) < n_steps) growth_rate <- rep(growth_rate, n_steps)
  if (length(delta_t) < n_steps) delta_t <- rep(delta_t, n_steps)
  
  # Create an instance of PrevPdf
  prevpdf <- PrevPdf$new(
    alpha = alpha_prior,
    beta = beta_prior,
    phi_prior = phi_prior,
    pi_seq = seq(1e-6, 1-1e-6,, length.out=pi_seq)
  )
  
  result <- prevpdf$compute_sample_size(alpha_intro, beta_intro, p_intro, growth_rate, rho, pi, dconf, delta_t, n_steps, method, n_max)
  
  return( result )
}
