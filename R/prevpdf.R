library(R6)

PrevPdf <- R6Class("PrevPdf",
                   public = list(
                     ts = NULL,
                     N = NULL,
                     alpha = NULL,
                     beta = NULL,
                     theta = NULL,
                     phi_prior = NULL,
                     phi_post = NULL,
                     rho = NULL,
                     f_prior = NULL,
                     f_prior_list = NULL,
                     f_posterior = NULL,
                     f_posterior_list = NULL,
                     pi_seq = NULL,
                     p_no_intro = NULL,
                     f_intro = NULL,
                     r = NULL,
                     deltaT = NULL,
                     initial_params = NULL,

                     # Helper function to expand parameters
                     .expand_param = function(param, num_steps, default_value) {
                       if (is.null(param)) {
                         # Use default value for all steps
                         rep(default_value, num_steps)
                       } else if (length(param) == 1) {
                         # Use the same value for all steps
                         rep(param, num_steps)
                       } else if (length(param) == num_steps) {
                         # Use provided values
                         param
                       } else {
                         stop("Parameter length must be 1 or equal to the number of time steps.")
                       }
                     },

                     initialize_state = function() {
                       # Reset fields to initial values
                       self$N <- integer(0)
                       self$alpha <- c(self$initial_params$alpha)
                       self$beta <- c(self$initial_params$beta)
                       self$phi_prior <- c(self$initial_params$phi_prior)
                       self$phi_post <- numeric(0)
                       self$rho <- numeric(0)
                       self$theta <- numeric(0)
                       self$pi_seq <- self$initial_params$pi_seq
                       self$r <- numeric(0)
                       self$deltaT <- numeric(0)
                       self$p_no_intro <- c(self$phi_prior)

                       # Initialize functions that depend on parameters
  
                       self$f_prior_list <- list(self$store(self$f_pi_pos(self$alpha[1], self$beta[1], self$phi_prior[1])))

                       self$f_posterior <- NULL

                       # Initialize the list to store posterior distributions for each step
                       self$f_posterior_list <- list()
                       
                       # Set time step
                       self$ts <- 1
                     },

                     initialize = function(alpha, beta, phi_prior, pi_seq) {
                       # Store initial parameters for reset
                       self$initial_params <- list(
                         alpha = alpha,
                         beta = beta,
                         phi_prior = phi_prior,
                         pi_seq = pi_seq
                       )

                       # Initialize the object's state
                       self$initialize_state()
                     },

                     reset = function() {
                       # Re-initialize the object's state
                       self$initialize_state()
                     },

                     update_params = function(params) {
                       valid_params <- c("alpha", "beta", "phi", "pi_seq")

                       if (!is.list(params) || length(params) == 0) {
                         stop("Please provide a list of parameters to update.")
                       }

                       for (param_name in names(params)) {
                         if (param_name %in% valid_params) {
                           param_value <- params[[param_name]]
                           # Update initial_params with the new value
                           self$initial_params[[param_name]] <- param_value
                         } else {
                           stop(sprintf("Unknown parameter name: %s", param_name))
                         }
                       }

                       # Re-initialize the object's state using the updated initial_params
                       self$initialize_state()
                     },

                     deep_copy = function() {
                       return(self$clone(deep = TRUE))
                     },


                     store = function(f) {
                       splinefun(self$pi_seq, f(self$pi_seq))
                     },

                     f_pi_pos = function(alpha, beta, phi) {
                       function(pi) (1 - phi) * dbeta(pi, alpha, beta)
                     },

                     likelihood = function(N, rho) {
                       function(pi) (1 - rho * pi) ^ N
                     },

                     get_theta = function(f, lik) {
                       integrate(function(pi) lik(pi) * f(pi), 0, 1)$value
                     },

                     phi_bayes_update = function(phi, theta_val) {
                       phi / (phi + theta_val)
                     },

                     f_pi_bayes_update = function(theta_val, phi_t, lik, f) {
                       function(pi) lik(pi) * f(pi) / (phi_t + theta_val)
                     },

                     Ginv = function(pi, r, deltaT) {
                       1 / (1 + (1 / pi - 1) * exp(r * deltaT))
                     },

                     dGinv = function(pi, r, deltaT) {
                       exp(r * deltaT) / (pi + (1 - pi) * exp(r * deltaT)) ^ 2
                     },

                     f_pi_time_update = function(f, r, deltaT) {
                       function(pi) self$dGinv(pi, r, deltaT) * f(self$Ginv(pi, r, deltaT))
                     },

                     f_pi_introduction_update = function(f_pi, phi_t, f_intro, p_no_intro) {
                       f <- function(pi_t_prime) {
                         # Attempt integration and catch any errors
                         integral_value <- tryCatch({
                           integrate(
                             function(x) f_pi(x) * f_intro((pi_t_prime - x) / (1 - x)) / (1 - x),
                             0,
                             pi_t_prime-5e-03
                           )$value
                         }, error = function(e) {
                           # Provide a more helpful error message
                           stop(
                             "Integration failed in f_pi_introduction_update. ",
                             "This may occur if pi_seq is too coarse or parameters lead to a non-integrable scenario.\n",
                             "Suggestion: Increase length.out in pi_seq or adjust model parameters.\n",
                             "Original error: ", conditionMessage(e)
                           )
                         })
                         # If integration succeeded, compute the final value
                         p_no_intro * f_pi(pi_t_prime) +
                           phi_t * f_intro(pi_t_prime) +
                           integral_value
                       }
                       function(x) vapply(x, f, 0)
                     },
                     update = function(N, alpha = NULL, beta = NULL, rho = NULL, r = NULL, deltaT = NULL, p_no_intro = NULL) {
                       if(length(N) > 1) stop("N must be length 1")
                       if(length(alpha) > 1) stop("alpha must be length 1")
                       if(length(beta) > 1) stop("beta must be length 1")
                       if (length(rho) > 1) stop("rho_i must be length 1")
                       if (length(r) > 1)   stop("r_i must be length 1")
                       if (length(deltaT) > 1) stop("deltaT_i must be length 1")
                       if (length(p_no_intro) > 1) stop("p_no_intro_i must be length 1")
                       
                       ts <- self$ts
                       
                       # Use values from previous timestep if not passed
                       alpha <- private$last_if_null(alpha, "alpha")
                       beta <- private$last_if_null(beta, "beta")
                       rho <- private$last_if_null(rho, "rho")
                       r <- private$last_if_null(r, "r")
                       deltaT <- private$last_if_null(deltaT, "deltaT")
                       p_no_intro <- private$last_if_null(p_no_intro, "p_no_intro")
                       
                       lik <- self$likelihood(N, rho)

                       self$theta[ts] <- self$get_theta(self$f_prior_list[[ts]], lik)

                       self$phi_post[ts] <- self$phi_bayes_update(self$phi_prior[ts], self$theta[ts])

                       # Bayes update
                       self$f_posterior <- self$store(
                         self$f_pi_bayes_update(self$theta[ts], self$phi_prior[ts], lik, self$f_prior_list[[ts]])
                       )
                       
                       # Record the Bayesian update params
                       self$N[ts] <- N
                       self$rho[ts] <- rho
                       
                       # Store f_posterior for this step
                       self$f_posterior_list[[ts]] <- self$f_posterior
                       
                       # Time (logistic) update
                       f_preintro <- self$store(self$f_pi_time_update(self$f_posterior, r, deltaT))
                       
                       # Record growth params
                       self$r[ts] <- r
                       self$deltaT[ts] <- deltaT

                       # Increment time step at the end
                       self$ts <- self$ts + 1
                       ts <- self$ts

                       # Create introduction distribution
                       f_intro <- self$f_pi_pos(alpha, beta, p_no_intro)

                       

                       # Apply introduction update
                       self$f_prior_list[[ts]] <- self$store(
                         self$f_pi_introduction_update(
                           f_preintro,
                           self$phi_post[ts - 1],
                           f_intro,
                           p_no_intro
                         )
                       )
                       self$phi_prior[ts] <- self$phi_post[ts - 1] * p_no_intro
                       
                       # Store parameters for introduction update
                       self$alpha[ts] <- alpha
                       self$beta[ts] <- beta
                       self$p_no_intro[ts] <- p_no_intro

                     },

                     compute_cdf = function(pi, step = NULL, prior=F) {
                       if (is.null(step)) step <- self$ts
                       # get zero inflation part
                       if (prior)
                       {
                         phi <- self$phi_prior[step]
                         fn <- self$f_prior_list[[step]]
                       }
                       else
                       {
                         phi <- self$phi_post[step]
                         fn <- self$f_posterior_list[[step]]
                       }
                       phi + integrate(fn, 0, pi)$value
                     },


                     compute_posterior_cdf_given_n = function(pi, N, rho = NULL) {
                       
                       rho <- private$last_if_null(rho, "rho")
                       
                       ts <- self$ts

                       if (ts == 1) {
                         phi_prior <- self$phi_prior[1]
                       } else {
                         phi_prior <- self$phi_prior[ts - 1]
                       }

                       lik <- self$likelihood(N, rho)
                       theta <- self$get_theta(self$f_prior_list[[ts]], lik)
                       phi_prior <- self$phi_prior[ts]
                       phi_post <- self$phi_bayes_update(phi_prior, theta)

                       # Bayes update
                       f_posterior <- self$store(
                         self$f_pi_bayes_update(theta, phi_prior, lik, self$f_prior_list[[ts]])
                       )

                       # Compute cumulative density
                       cdf <- phi_post + integrate(f_posterior, 0, pi)$value

                       return(cdf)
                     },
                     
                     # NOTE that all parameters up to "n_steps" must be vectors
                     # of length n_steps. 
                     # !!!! No input validation is done at this stage !!!!
                     compute_sample_size = function(alpha_intro, beta_intro, p_intro, growth_rate, rho, pi, dconf, delta_t, n_steps, method) {
                       n_required <- rep(NA, n_steps)
                       p_eff_freedom_prior <- rep(NA, n_steps)
                       p_eff_freedom_post  <- rep(NA, n_steps)
                       sensitivity <- rep(NA, n_steps)
                       for (j in 1:n_steps)
                       {
                         threshold_quantile <- dconf[j]
                         if (method == "maintain") {
                           threshold_quantile <- min(c(0.999, threshold_quantile / (1 - p_intro[j])))
                         }
                         n_required[j] <- self$n_from_cdf(threshold_quantile, pi[j], rho[j])
                         self$update(n_required[j], alpha_intro[j], beta_intro[j], rho[j], growth_rate[j], delta_t[j], 1 - p_intro[j])
                         p_eff_freedom_prior[j] <- self$compute_cdf(pi[j], step=j, prior=T)
                         p_eff_freedom_post[j] <- self$compute_cdf(pi[j], step=j)
                         sensitivity[j] <- self$sensitivity(j)
                       }
                       
                       result <- list(n                   = n_required,
                                      pi                  = pi,
                                      phi_prior           = self$phi_prior,
                                      phi_post            = self$phi_post,
                                      p_freedom_prior     = p_eff_freedom_prior,
                                      p_freedom_post      = p_eff_freedom_post,
                                      sensitivity         = sensitivity,
                                      f_prior             = self$f_prior_list,
                                      f_posterior         = self$f_posterior_list)
                       class(result) <- "diseasefree"
                       return( result )
                     },
                     
                     # NOTE all parameters must be vectors of equal length
                     # !!!! No input validation is done at this stage !!!!
                     compute_probability_of_freedom = function(n, alpha_intro, beta_intro, p_intro, growth_rate, rho, pi, delta_t) {
                       n_steps <- length(n)
                       p_eff_freedom_prior <- rep(NA, n_steps)
                       p_eff_freedom_post  <- rep(NA, n_steps)
                       sensitivity <- rep(NA, n_steps)
                       for (j in 1:n_steps)
                       {
                         self$update(n[j], alpha_intro[j], beta_intro[j], rho[j], growth_rate[j], delta_t[j], 1 - p_intro[j])
                         p_eff_freedom_prior[j] <- self$compute_cdf(pi[j], step=j, prior=T)
                         p_eff_freedom_post[j] <- self$compute_cdf(pi[j], step=j)
                         sensitivity[j] <- self$sensitivity(j)
                       }
                       
                       result <- list(n                   = n,
                                      pi                  = pi,
                                      phi_prior           = self$phi_prior,
                                      phi_post            = self$phi_post,
                                      p_freedom_prior     = p_eff_freedom_prior,
                                      p_freedom_post      = p_eff_freedom_post,
                                      sensitivity         = sensitivity,
                                      f_prior             = self$f_prior_list,
                                      f_posterior         = self$f_posterior_list)
                       class(result) <- "diseasefree"   
                       return( result )
                     },

                     n_from_cdf = function(q, pi_value, rho = NULL, n_max = 1000) {
                       
                       rho <- private$last_if_null(rho, "rho")
                       
                       low <- 1
                       high <- n_max

                       while (low <= high) {
                         mid <- floor((low + high) / 2)

                         # Compute the CDF at pi_value for mid samples
                         cdf_value <- self$compute_posterior_cdf_given_n(pi_value, mid, rho)

                         if (cdf_value < q) {
                           low <- mid + 1
                         } else {
                           high <- mid - 1
                         }
                       }

                       return(low)
                     },
                     
                     sensitivity = function(ts = NULL) {
                       if (is.null(ts)) {
                         sens <- 1-self$theta / (1-self$phi_prior)
                       } else {
                         sens <- 1-self$theta[ts] / (1-self$phi_prior[ts])
                       }
                       sens
                     }                     
                   ),
                   
                   private = list(
                     last_if_null = function(x, param) {
                       if(is.null(x)) {
                         if(length(self[[param]]) > 0) {
                           # Use last param value
                           self[[param]][self$ts]
                         } else {
                           # Initial step and no param specified
                           stop(paste0(param, " must be specified"))
                         }
                       } else {
                         # Param specified so return it
                         x
                       }
                       
                     }                       
                   )
)



if(F) {

  # Example usage ----
  # Define the parameters
  phi_1 <- 0.5
  alpha <- 1
  beta <- 8.06
  rho <- 1
  deltaT <- 1
  r <- 0
  p_intro <- 0#0.001
  pi_seq <- seq(0, 1, length.out = 1000)
  N_samp <- c(146,0,0)  # Example sample sizes at each time step

  # Create an instance of PrevPdf
  prevpdf <- PrevPdf$new(
    alpha = alpha,
    beta = beta,
    phi_prior = phi_1,
    pi_seq = pi_seq
  )

  prevpdf$compute_posterior_cdf_given_n(.1, 10, rho)
  # Update the instance over time steps
  for (N in N_samp) {
    prevpdf$update(N, alpha = alpha, beta = beta, rho = rho, r = r, deltaT = deltaT, p_no_intro = 1-p_intro)
  }
  print(prevpdf)

  # Reset the prev pdf timeseries
  prevpdf$reset()

  # Desired CDF value (confidence level) and prevalence value pi
  desired_cdf_level <- 0.95
  pi_value <- 0.0  # For example, testing for disease freedom


  for (N in 1:3) {
    n_required <- prevpdf$n_from_cdf(desired_cdf_level, pi_value, rho)
    prevpdf$update(n_required, alpha = alpha, beta = beta, rho = rho, r = r, deltaT = deltaT, p_no_intro = 1-p_intro)
  }
  print(prevpdf)
  # Find the number of samples required
  n_required <- prevpdf$n_from_cdf(desired_cdf_level, pi_value, rho)

  cat(sprintf("Number of samples required to achieve CDF >= %.2f at pi = %.2f: %d\n", desired_cdf_level, pi_value, n_required))

}
