source("prevpdf.R")


is_valid_length <- function(x, min_size, max_size)
{
  return( length(x) == min_size || length(x) == max_size )
}

# phi_prior    - prior probability of freedom
# alpha        - beta distribution alpha shape parameter
# beta         - beta distribution beta shape parameter
# p_intro      - probability of introduction
# rho          - unit test sensitivity
# pi           - design threshold
# dconf        - desired confidence
# growth_rate  - exponential growth rate
# delta_t      - test time duration
# n_steps      - number of time units to step
# method       - "restore" or "maintain"
# pi_seq       - discretization granularity
compute_sample_size <- function(phi_prior, alpha, beta, p_intro, growth_rate, rho, pi, dconf, delta_t=1, n_steps=1, method="restore", pi_seq=1000)
{
  # determine max input vector length
  max_len = max( c(length(phi_prior), 
                   length(alpha),
                   length(beta),
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
  if (!is_valid_length(phi_prior, 1, n_steps)) stop(simpleError("Vector length mismatch: 'phi_prior'"))
  if (!is_valid_length(alpha, 1, n_steps)) stop(simpleError("Vector length mismatch: 'alpha'"))
  if (!is_valid_length(beta, 1, n_steps)) stop(simpleError("Vector length mismatch: 'beta'"))
  if (!is_valid_length(p_intro, 1, n_steps)) stop(simpleError("Vector length mismatch: 'p_intro'"))
  if (!is_valid_length(rho, 1, n_steps)) stop(simpleError("Vector length mismatch: 'rho'"))
  if (!is_valid_length(dconf, 1, n_steps)) stop(simpleError("Vector length mismatch: 'dconf'"))
  if (!is_valid_length(pi, 1, n_steps)) stop(simpleError("Vector length mismatch: 'pi'"))
  if (!is_valid_length(growth_rate, 1, n_steps)) stop(simpleError("Vector length mismatch: 'growth_rate'"))
  if (!is_valid_length(delta_t, 1, n_steps)) stop(simpleError("Vector length mismatch: 'delta_t'"))
  
  # make all vectors the same length (matching number of steps)
  if (length(phi_prior) < n_steps) phi_prior <- rep(phi_prior, n_steps)
  if (length(alpha) < n_steps) alpha <- rep(alpha, n_steps)
  if (length(beta) < n_steps) beta <- rep(beta, n_steps)
  if (length(p_intro) < n_steps) p_intro <- rep(p_intro, n_steps)
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
    pi_seq = pi_seq,
    r = growth_rate[1],
    deltaT = delta_t[1],
    p_no_intro = 1 - p_intro[1]
  )
  
  n_required <- rep(NA, n_steps)
  for (j in 1:n_steps)
  {
    n_required[j] <- prevpdf$n_from_cdf(dconf[j], pi[j])
    if (j < n_steps)
    {
      # update needs dynamic rho, r, delta_t, p_no_intro
      prevpdf$update(n_required[j], alpha[j+1], beta[j+1])
    }
  }
  
  result <- list(n_required=n_required)
  class(result) <- "diseasefree"
  return( result )
}

compute_probability_of_freedom <- function()
{
  
}