library(R6)
library(dplyr)
library(ggplot2)
library(gganimate)
library(gifski)

source("R/prevpdf.R")

PrevPdfExt <- R6Class("PrevPdfExt",
                   inherit=PrevPdf,
                   public=list(
                     apply_sampling_schedule = function(N_samp, alpha = NULL, beta = NULL, rho = NULL, r = NULL, deltaT = NULL, p_no_intro = NULL) {
                       num_steps <- length(N_samp)
                       
                       alpha <- private$last_if_null(alpha, "alpha")
                       beta <- private$last_if_null(beta, "beta")
                       rho <- private$last_if_null(rho, "rho")
                       r <- private$last_if_null(r, "r")
                       deltaT <- private$last_if_null(deltaT, "deltaT")
                       p_no_intro <- private$last_if_null(p_no_intro, "p_no_intro")
                       
                       # Handle alpha and beta inputs
                       alpha_list <- self$.expand_param(alpha, num_steps, NA)
                       beta_list <- self$.expand_param(beta, num_steps, NA)
                       rho_list <- self$.expand_param(rho, num_steps, NA)
                       r_list <- self$.expand_param(r, num_steps, NA)
                       deltaT_list <- self$.expand_param(deltaT, num_steps, NA)
                       p_no_intro_list <- self$.expand_param(p_no_intro, num_steps, NA)
                       
                       for (i in seq_len(num_steps)) {
                         N <- N_samp[i]
                         self$update(N, alpha = alpha_list[i], beta = beta_list[i], rho = rho_list[i], r = r_list[i], deltaT = deltaT_list[i], p_no_intro = p_no_intro_list[i])
                       }
                     },
                     
                     determine_schedule_counts = function(desired_cdf_level, pi_value, sampling_timing = NULL, num_steps = NULL, alpha = NULL, beta = NULL, rho = NULL, r = NULL, deltaT = NULL, p_no_intro = NULL) {
                       if (is.null(sampling_timing)) {
                         if (is.null(num_steps)) {
                           stop("Please provide 'num_steps' when 'samp_sched' is NULL.")
                         }
                         # Sampling occurs at every time point
                         sampling_timing <- rep(1, num_steps)  # Sample at every time step
                       }
                       num_steps <- length(sampling_timing)
                       
                       alpha <- private$last_if_null(alpha, "alpha")
                       beta <- private$last_if_null(beta, "beta")
                       rho <- private$last_if_null(rho, "rho")
                       r <- private$last_if_null(r, "r")
                       deltaT <- private$last_if_null(deltaT, "deltaT")
                       p_no_intro <- private$last_if_null(p_no_intro, "p_no_intro")
                       
                       # Handle alpha and beta inputs
                       alpha_list <- self$.expand_param(alpha, num_steps, NA)
                       beta_list <- self$.expand_param(beta, num_steps, NA)
                       rho_list <- self$.expand_param(rho, num_steps, NA)
                       r_list <- self$.expand_param(r, num_steps, NA)
                       deltaT_list <- self$.expand_param(deltaT, num_steps, NA)
                       p_no_intro_list <- self$.expand_param(p_no_intro, num_steps, NA)
                       
                       for (i in seq_len(num_steps)) {
                         if (sampling_timing[i] == 0) {
                           n_required <- 0
                         } else {
                           n_required <- self$n_from_cdf(desired_cdf_level, pi_value, rho_list[i])
                         }
                         self$update(n_required, alpha = alpha_list[i], beta = beta_list[i], rho = rho_list[i], r = r_list[i], deltaT = deltaT_list[i], p_no_intro = p_no_intro_list[i])
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
  N_samp <- c(100,50,20)  # Example sample sizes at each time step
  
  # Create an instance of PrevPdf
  prevpdf <- PrevPdfExt$new(
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

  prevpdf$prevalence_pdf_seq_animation()
}
