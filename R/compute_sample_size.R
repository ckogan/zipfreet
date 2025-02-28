source("R/is_valid_length.R")
source("R/prevpdf.R")

#' Compute the sample size requirements to meet a specified design
#' threshold at a specified confidence.
#' 
#' @param phi_prior prior probability of freedom
#' @param alpha beta distribution alpha shape parameter
#' @param beta  beta distribution beta shape parameter
#' @param p_intro probability of introduction
#' @param rho unit test sensitivity
#' @param pi design threshold
#' @param dconf desired confidence
#' @param growth_rate exponential growth rate
#' @param delta_t test time duration
#' @param n_steps number of time units to step
#' @param method  "restore" or "maintain" - Note that when using "maintain"
#'                with a time varying p_intro, you must provide p_intro
#'                for one step beyond the calculation horizon.
#' @param pi_seq discretization granularity
#' @returns "diseasefree" structure including sample sizes and prior/posterior distributions
#' @examples
#' u <- compute_sample_size(0.5, 1, 1, 0.04, 0.01, 0.9, 0, 0.95, n_steps=5)
#' summary(u)
#' plot(u)
#' @export
compute_sample_size <- function(phi_prior, alpha, beta, p_intro, growth_rate, rho, pi, dconf, delta_t=1, n_steps=1, method="restore", pi_seq=1000)
{
  # phi_prior must be length 1
  if (length(phi_prior) > 1)
  {
    stop( simpleError("phi_prior must be length 1"))
  }
  
  # determine max input vector length
  max_len = max( c(length(alpha),
                   length(beta),
                   length(p_intro),
                   length(rho),
                   length(dconf),
                   length(pi),
                   length(growth_rate),
                   length(delta_t)) )
  
  # make sure the method is either "restore" or "maintain"
  if (method != "restore" && method != "maintain")
  {
    stop( simpleError("Invalid method; must be either 'restore' or 'maintain'"))
  }
  
  # the max vector size drives the number of steps to take
  if (max_len > 1)
  {
    n_steps = max_len
  }
  
  # if the method is "maintain", then p_intro needs to be one entry longer
  # than the number of steps, since we need p_intro for the subsequent
  # step to make the threshold adjustment
  if (method == "maintain")
  {
    p_intro_len_adj <- 1
  }
  else
  {
    p_intro_len_adj <- 0
  }
  
  # if any input vector has length > 1 but < n_steps, error out
  if (!is_valid_length(alpha, 1, n_steps)) stop(simpleError("Vector length mismatch: 'alpha'"))
  if (!is_valid_length(beta, 1, n_steps)) stop(simpleError("Vector length mismatch: 'beta'"))
  if (!is_valid_length(p_intro, 1, n_steps + p_intro_len_adj)) stop(simpleError("Vector length mismatch: 'p_intro'"))
  if (!is_valid_length(rho, 1, n_steps)) stop(simpleError("Vector length mismatch: 'rho'"))
  if (!is_valid_length(dconf, 1, n_steps)) stop(simpleError("Vector length mismatch: 'dconf'"))
  if (!is_valid_length(pi, 1, n_steps)) stop(simpleError("Vector length mismatch: 'pi'"))
  if (!is_valid_length(growth_rate, 1, n_steps)) stop(simpleError("Vector length mismatch: 'growth_rate'"))
  if (!is_valid_length(delta_t, 1, n_steps)) stop(simpleError("Vector length mismatch: 'delta_t'"))
  
  # make all vectors the same length (matching number of steps)
  if (length(alpha) < n_steps) alpha <- rep(alpha, n_steps)
  if (length(beta) < n_steps) beta <- rep(beta, n_steps)
  if (length(p_intro) < n_steps + p_intro_len_adj) p_intro <- rep(p_intro, n_steps + p_intro_len_adj)
  if (length(rho) < n_steps) rho <- rep(rho, n_steps)
  if (length(dconf) < n_steps) dconf <- rep(dconf, n_steps)
  if (length(pi) < n_steps) pi <- rep(pi, n_steps)
  if (length(growth_rate) < n_steps) growth_rate <- rep(growth_rate, n_steps)
  if (length(delta_t) < n_steps) delta_t <- rep(delta_t, n_steps)
  
  # Create an instance of PrevPdf
  prevpdf <- PrevPdf$new(
    alpha = alpha[1],
    beta = beta[1],
    phi = phi_prior,
    rho = rho[1],
    pi_seq = seq(0, 1, length.out=pi_seq),
    r = growth_rate[1],
    deltaT = delta_t[1],
    p_no_intro = 1 - p_intro[1]
  )
  
  n_required <- rep(NA, n_steps)
  p_eff_freedom <- rep(NA, n_steps)
  for (j in 1:n_steps)
  {
    threshold_quantile <- dconf[j]
    if (method == "maintain") {
      threshold_quantile <- min(c(0.999, threshold_quantile / (1 - p_intro[j+1])))
    }
    n_required[j] <- prevpdf$n_from_cdf(threshold_quantile, pi[j])
    # update needs dynamic rho, r, delta_t, p_no_intro
    prevpdf$update(n_required[j], alpha[j], beta[j])
    p_eff_freedom[j] <- prevpdf$compute_cdf(pi[j], step=j)
  }
  
  result <- list(n             = n_required,
                 phi_prior     = prevpdf$phi_prior,
                 phi_post      = prevpdf$phi_post,
                 p_eff_freedom = p_eff_freedom,
                 f_posterior   = prevpdf$f_posterior_list)
  class(result) <- "diseasefree"
  return( result )
}
