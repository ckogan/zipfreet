library(testthat)
library(here)
context("Testing PrevPdf Class")

test_that("Correct phi_post values over 3 time points", {
  prevpdf_test <- PrevPdfExt$new(
    alpha = 1,
    beta = 8,
    phi = 0.5,
    pi_seq = seq(0, 1, length.out = 1000)
  )

  # Check initial time step
  expect_equal(prevpdf_test$ts, 1)
  expect_equal(prevpdf_test$phi_prior[1], 0.5)

  # Perform updates with known sample sizes
  N_samp <- c(146,  25, 15)
  prevpdf_test$apply_sampling_schedule(N_samp, alpha = 1, beta = 8, rho = 0.92, r = 1.2, deltaT = 1, p_no_intro = 1-0.03)

  # Expected phi_pre values
  # dput(prevpdf_test$phi_prior)
  expected_phi_prior <- c(0.5, 0.918350671522885, 0.928201564259836, 0.935491512489105)

  # Expected phi_post values
  # dput(prevpdf_test$phi_post)
  expected_phi_post <- c(0.946753269611221, 0.956908829133852, 0.96442423967949)

  # Compare phi_pre to expected values
  expect_equal(prevpdf_test$phi_prior, expected_phi_prior, tolerance = 1e-6)

  # Compare phi_post to expected values
  expect_equal(prevpdf_test$phi_post, expected_phi_post, tolerance = 1e-6)

})

test_that("determine_schedule_counts produces correct sample counts and phi_post", {
  # Initialize
  prevpdf_test <- PrevPdfExt$new(
    alpha = 1,
    beta = 8,
    phi = 0.5,
    pi_seq = seq(0, 1, length.out = 1000)
  )

  # Define sampling timing
  sampling_timing <- c(1, 0, 1)
  desired_cdf_level <- 0.95
  pi_value <- 0.0

  # Apply determine_schedule_counts
  prevpdf_test$determine_schedule_counts(
    desired_cdf_level = desired_cdf_level,
    pi_value = pi_value,
    sampling_timing = sampling_timing,
    alpha = 1, beta = 8, rho = 0.92, r = 1.2, deltaT = 1, p_no_intro = 1-0.03
  )

  # Expected sample counts (n_required values)
  expected_N <- c(157, 0, 11)  # Assuming n_from_cdf returns 146

  # Compare N to expected sample counts
  expect_equal(prevpdf_test$N, expected_N)

  # Expected phi_post values
  # dput(prevpdf_test$phi_post)
  expected_phi_post <- c(0.950113398975309, 0.921609936062774, 0.951045294135913)

  # Compare phi_post to expected values
  expect_equal(prevpdf_test$phi_post, expected_phi_post, tolerance = 1e-6)
})

test_that("Special Case: Unit sensitivity, no growth, alpha=1", {
  # Initialize
  prevpdf_test <- PrevPdfExt$new(
    alpha = 1,
    beta = 8,
    phi = 0.5,
    pi_seq = seq(0, 1, length.out = 1000)
  )
  
  prevpdf_test$apply_sampling_schedule(150, rho = 1, r = 0, deltaT = 1, p_no_intro = 1)
  
  
  # Expected phi_post values
  phi1 <- 0.5
  b <- 8
  Nsamp <- 150
  
  # shell.exec(here('obsidian/P119 Derivation of Posterior k Assume Beta 1 b 2023-11-30.pdf')) # see pg2
  expected_phi_post <- phi1 / (phi1 + (1-phi1) * beta(1,b+Nsamp)/beta(1,b))

  # Compare phi_post to expected value
  expect_equal(prevpdf_test$phi_post, expected_phi_post)

})

