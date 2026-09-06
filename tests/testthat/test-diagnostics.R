test_that("log-ratio validation requires finite numeric input", {
  expect_invisible(upsis:::validate_log_ratios(c(0, 1, -1)))

  expect_error(
    upsis:::validate_log_ratios(character()),
    "non-empty numeric vector"
  )
  expect_error(
    upsis:::validate_log_ratios(c(0, NA_real_)),
    "finite and non-missing"
  )
  expect_error(
    upsis:::validate_log_ratios(c(0, Inf)),
    "finite and non-missing"
  )
})

test_that("upsis result prints compact diagnostics", {
  result <- upsis:::new_upsis_result(updated_model = list(), pareto_k = NA_real_)

  expect_s3_class(result, "upsis")
  expect_output(print(result), "upsis update")
  expect_output(print(result), "Pareto k: not computed")
})
