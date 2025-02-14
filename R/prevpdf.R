library(R6)
library(dplyr)
library(ggplot2)
library(gganimate)
library(gifski)

PrevPdf <- R6Class("PrevPdf",
                   public = list(
                     ts = 0,
                     N = NULL,
                     alpha = NULL,
                     beta = NULL,
                     alpha_i = NULL,
                     beta_i = NULL,
                     theta = NULL,
                     phi_prior = NULL,
                     phi_post = NULL,
                     rho = NULL,
                     f_prior = NULL,
                     f_posterior = NULL,
                     f_posterior_list = NULL,
                     f_preintro = NULL,
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
                       self$N <- c()
                       self$alpha <- self$initial_params$alpha
                       self$beta <- self$initial_params$beta
                       # Initialize alpha_i and beta_i with initial alpha and beta
                       self$alpha_i <- c(self$alpha)
                       self$beta_i <- c(self$beta)
                       self$phi_prior <- c(self$initial_params$phi)
                       self$phi_post <- c()
                       self$rho <- self$initial_params$rho
                       self$theta <- c()
                       self$pi_seq <- self$initial_params$pi_seq
                       self$r <- self$initial_params$r
                       self$deltaT <- self$initial_params$deltaT
                       self$p_no_intro <- self$initial_params$p_no_intro

                       # Initialize functions that depend on parameters
                       self$f_intro <- self$f_pi_pos(self$alpha, self$beta, self$p_no_intro)
                       self$f_prior <- self$store(self$f_pi_pos(self$alpha, self$beta, self$phi_prior[1]))
                       self$f_posterior <- NULL
                       self$f_preintro <- NULL

                       # Initialize the list to store posterior distributions for each step
                       self$f_posterior_list <- list()
                       
                       # Set time step
                       self$ts <- 1
                     },

                     initialize = function(alpha, beta, phi, rho, pi_seq, r, deltaT, p_no_intro) {
                       # Store initial parameters for reset
                       self$initial_params <- list(
                         alpha = alpha,
                         beta = beta,
                         phi = phi,
                         rho = rho,
                         pi_seq = pi_seq,
                         r = r,
                         deltaT = deltaT,
                         p_no_intro = p_no_intro
                       )

                       # Initialize the object's state
                       self$initialize_state()
                     },

                     reset = function() {
                       # Re-initialize the object's state
                       self$initialize_state()
                     },

                     update_params = function(params) {
                       valid_params <- c("alpha", "beta", "phi", "rho", "pi_seq", "r", "deltaT", "p_no_intro")

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

                     likelihood = function(N) {
                       function(pi) (1 - self$rho * pi) ^ N
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

                     Ginv = function(pi) {
                       1 / (1 + (1 / pi - 1) * exp(self$r * self$deltaT))
                     },

                     dGinv = function(pi) {
                       exp(self$r * self$deltaT) / (pi + (1 - pi) * exp(self$r * self$deltaT)) ^ 2
                     },

                     f_pi_time_update = function(f) {
                       function(pi) self$dGinv(pi) * f(self$Ginv(pi))
                     },

                     f_pi_introduction_update = function(f_pi, phi_t, f_intro) {
                       f <- function(pi_t_prime) {
                         # Attempt integration and catch any errors
                         integral_value <- tryCatch({
                           integrate(
                             function(x) f_pi(x) * f_intro((pi_t_prime - x) / (1 - x)) / (1 - x),
                             0,
                             pi_t_prime
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
                         self$p_no_intro * f_pi(pi_t_prime) +
                           phi_t * f_intro(pi_t_prime) +
                           integral_value
                       }
                       function(x) vapply(x, f, 0)
                     },
                     

                     update = function(N, alpha_i = NULL, beta_i = NULL) {

                       ts <- self$ts
                       # Record the sample size
                       self$N <- c(self$N, N)

                       # Use provided alpha_i and beta_i or default to prior values
                       if (is.null(alpha_i)) {
                         alpha_i <- self$alpha
                       }
                       if (is.null(beta_i)) {
                         beta_i <- self$beta
                       }

                       lik <- self$likelihood(N)

                       self$theta[ts] <- self$get_theta(self$f_prior, lik)
                       self$phi_post[ts] <- self$phi_bayes_update(self$phi_prior[ts], self$theta[ts])

                       # Bayes update
                       self$f_posterior <- self$store(
                         self$f_pi_bayes_update(self$theta[ts], self$phi_prior[ts], lik, self$f_prior)
                       )
                       
                       # Store f_posterior for this step
                       self$f_posterior_list[[ts]] <- self$f_posterior
                       
                       # Time (logistic) update
                       self$f_preintro <- self$store(self$f_pi_time_update(self$f_posterior))

                       # Increment time step at the end
                       self$ts <- self$ts + 1
                       ts <- self$ts

                       # Create introduction distribution
                       f_intro_i <- self$f_pi_pos(alpha_i, beta_i, self$p_no_intro)

                       # Store alpha_i and beta_i for this time step
                       self$alpha_i[ts] <- alpha_i
                       self$beta_i[ts] <- beta_i

                       # Apply introduction update
                       self$f_prior <- self$store(
                         self$f_pi_introduction_update(
                           self$f_preintro,
                           self$phi_post[ts - 1],
                           f_intro_i
                         )
                       )
                       self$phi_prior[ts] <- self$phi_post[ts - 1] * self$p_no_intro

                     },

                     compute_cdf = function(x, step = NULL) {
                       if (is.null(step)) step <- self$ts
                       # get zero inflation part
                       phi <- self$phi_post[step]
                       f_posterior_current <- self$f_posterior_list[[step]]

                       numerator <- phi + integrate(f_posterior_current, 0, x)$value
                       denominator <- phi + integrate(f_posterior_current, 0, 1)$value
                       
                       numerator / denominator
                     },


                     compute_posterior_cdf_given_n = function(pi, N) {
                       ts <- self$ts

                       if (ts == 1) {
                         phi_prior <- self$phi_prior[1]
                       } else {
                         phi_prior <- self$phi_prior[ts - 1]
                       }

                       lik <- self$likelihood(N)
                       theta <- self$get_theta(self$f_prior, lik)
                       phi_prior <- self$phi_prior[ts]
                       phi_post <- self$phi_bayes_update(phi_prior, theta)

                       # Bayes update
                       f_posterior <- self$store(
                         self$f_pi_bayes_update(theta, phi_prior, lik, self$f_prior)
                       )

                       # Compute cumulative density
                       cdf <- phi_post + integrate(f_posterior, 0, pi)$value

                       return(cdf)
                     },

                     n_from_cdf = function(q, pi_value, n_max = 1000) {
                       low <- 1
                       high <- n_max

                       while (low <= high) {
                         mid <- floor((low + high) / 2)

                         # Compute the CDF at pi_value for mid samples
                         cdf_value <- self$compute_posterior_cdf_given_n(pi_value, mid)

                         if (cdf_value < q) {
                           low <- mid + 1
                         } else {
                           high <- mid - 1
                         }
                       }

                       return(low)
                     },

                     apply_sampling_schedule = function(N_samp, alpha_i = NULL, beta_i = NULL) {
                       num_steps <- length(N_samp)

                       # Handle alpha_i and beta_i inputs
                       alpha_list <- self$.expand_param(alpha_i, num_steps, self$alpha)
                       beta_list <- self$.expand_param(beta_i, num_steps, self$beta)

                       for (i in seq_len(num_steps)) {
                         N <- N_samp[i]
                         self$update(N, alpha_i = alpha_list[i], beta_i = beta_list[i])
                       }
                     },

                     determine_schedule_counts = function(desired_cdf_level, pi_value, sampling_timing = NULL, num_steps = NULL, alpha_i = NULL, beta_i = NULL) {
                       if (is.null(sampling_timing)) {
                         if (is.null(num_steps)) {
                           stop("Please provide 'num_steps' when 'samp_sched' is NULL.")
                         }
                         # Sampling occurs at every time point
                         sampling_timing <- rep(1, num_steps)  # Sample at every time step
                       }

                       num_steps <- length(sampling_timing)

                       # Handle alpha_i and beta_i inputs
                       alpha_list <- self$.expand_param(alpha_i, num_steps, self$alpha)
                       beta_list <- self$.expand_param(beta_i, num_steps, self$beta)

                       for (i in seq_len(num_steps)) {
                         if (sampling_timing[i] == 0) {
                           self$update(0, alpha_i = alpha_list[i], beta_i = beta_list[i])
                         } else {
                           n_required <- self$n_from_cdf(desired_cdf_level, pi_value)
                           self$update(n_required, alpha_i = alpha_list[i], beta_i = beta_list[i])
                         }
                       }

                     },

                     get_freedom = function(only_positive = FALSE, design_prevalence = NULL) {
                       if(length(self$phi_post) == 0) {
                         df <- data.frame(
                           `time step` = integer(0),
                           N           = integer(0),
                           `p(free,rr)` = numeric(0),
                           `p(free)` = numeric(0)
                         )
                       } else {
                         timeseq <- seq_along(self$phi_post)
                         df <- data.frame(
                           `time step` = timeseq,
                           N           = self$N[timeseq],
                           `p(free,rr)` = self$phi_prior[timeseq],
                           `p(free)` = self$phi_post[timeseq]
                         )
                         if (only_positive) {
                           df <- df[df$N > 0, ]
                         }
                         
                         if (!is.null(design_prevalence)) {
                           # Compute probability that disease ≤ design_prevalence at each step
                           cdf_values <- sapply(df$time.step, function(step) {
                             self$compute_cdf(design_prevalence, step = step)
                           })
                           df$`p(disease ≤ design_prevalence)` <- cdf_values
                         }
                       }
                       
                       return(df)
                     },
                     
                     prevalence_pdf_seq_df = function() {
                       
                       make_df <- function(i, pi_seq) {
                         data.frame(
                           value = self$pi_seq,
                           time = i,
                           density = self$f_posterior_list[[i]](self$pi_seq)
                         )
                       }
                       
                       prev_df <- lapply(seq_along(self$f_posterior_list), function(i) make_df(i, pi_seq)) %>% bind_rows()
                       
                       prev_df
                     },
                     
                     prevalence_pdf_seq_animation = function() {
                       prev_df <- self$prevalence_pdf_seq_df()
                       
                       p <- ggplot(prev_df, aes(x = value, y = density, group = time)) +
                         geom_line() +
                         labs(title = 'Time: {frame_time}') +
                         transition_time(time) + theme_minimal() + 
                         labs(y = "Probability density", x = "Prevalence") + 
                         scale_y_continuous(limits = c(0, 1))
                       
                       # Animate the plot
                       anim <- animate(p, nframes = length(self$f_posterior_list), renderer = gifski_renderer())
                       anim
                     },
                     
                     prevalence_pdf_mean = function() {
                       mean_prev <- sapply(self$f_posterior_list, function(f) {
                         integrate(function(x) x * f(x), 0, 1)$value
                       })
                       
                       data.frame(
                         time = seq_along(mean_prev),
                         mean_prev = mean_prev
                       )
                     }

                   #   print = function(...) {
                   #     # Get the freedom data
                   #     freedom_data <- self$get_freedom()
                   #
                   #     if (nrow(freedom_data) == 0) {
                   #       cat("PrevPdf Object Summary:\n")
                   #       cat("No time steps have been performed yet.\n")
                   #     } else {
                   #       # Print the data frame
                   #       print(freedom_data)
                   #     }
                   #
                   #     # Return the object invisibly
                   #     invisible(self)
                   #   }
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
    phi = phi_1,
    rho = rho,
    pi_seq = pi_seq,
    r = r,
    deltaT = deltaT,
    p_no_intro = 1 - p_intro
  )

  prevpdf$compute_posterior_cdf_given_n(.1, 10)
  # Update the instance over time steps
  for (N in N_samp) {
    prevpdf$update(N)
  }
  print(prevpdf)

  # Reset the prev pdf timeseries
  prevpdf$reset()

  # Desired CDF value (confidence level) and prevalence value pi
  desired_cdf_level <- 0.95
  pi_value <- 0.0  # For example, testing for disease freedom


  for (N in 1:3) {
    n_required <- prevpdf$n_from_cdf(desired_cdf_level, pi_value)
    prevpdf$update(n_required)

  }
  print(prevpdf)
  # Find the number of samples required
  n_required <- prevpdf$n_from_cdf(desired_cdf_level, pi_value)

  cat(sprintf("Number of samples required to achieve CDF >= %.2f at pi = %.2f: %d\n", desired_cdf_level, pi_value, n_required))

}
