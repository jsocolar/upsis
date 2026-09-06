test_that("brmsfit updates run end to end when explicitly enabled", {
  skip_if_not(
    identical(tolower(Sys.getenv("RUN_BRMS_TESTS")), "true"),
    "Set RUN_BRMS_TESTS=true to run brms integration tests."
  )
  skip_if_not_installed("brms")

  set.seed(1)
  initial_data <- data.frame(x = rnorm(8))
  initial_data$y <- rnorm(nrow(initial_data), initial_data$x)
  added_data <- data.frame(x = rnorm(3))
  added_data$y <- rnorm(nrow(added_data), added_data$x)

  fit <- brms::brm(
    y ~ x,
    data = initial_data,
    chains = 1,
    iter = 200,
    warmup = 100,
    refresh = 0,
    seed = 1
  )

  result <- upsis(fit, data_add = added_data)

  expect_s3_class(result, "upsis")
  expect_s3_class(result$updated_model, "brmsfit")
  expect_true(is.numeric(result$pareto_k))
  expect_length(result$weights, brms::ndraws(fit))
})
