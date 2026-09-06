test_that("weights are split by chain in draw order", {
  split_weights <- upsis:::split_weights_by_chain(1:6, n_chains = 2)

  expect_equal(split_weights, list(1:3, 4:6))
})

test_that("weight splitting validates chain compatibility", {
  expect_error(
    upsis:::split_weights_by_chain(1:5, n_chains = 2),
    "divisible by the number of chains"
  )
  expect_error(
    upsis:::split_weights_by_chain(1:5, n_chains = 0),
    "positive scalar"
  )
})

test_that("warmup draws per chain are inferred from stored and post-warmup draws", {
  draws <- draws_list_fixture(n_iterations = 5, n_chains = 2)

  expect_equal(
    upsis:::n_warmup_draws_per_chain(draws, n_post_warmup_draws = 6),
    2
  )
})

test_that("warmup accounting validates draw counts", {
  draws <- draws_list_fixture(n_iterations = 5, n_chains = 2)

  expect_error(
    upsis:::n_warmup_draws_per_chain(draws, n_post_warmup_draws = 5),
    "evenly divisible across chains"
  )
  expect_error(
    upsis:::n_warmup_draws_per_chain(draws, n_post_warmup_draws = 12),
    "exceeds the stored draw count"
  )
})

test_that("one-chain resampling preserves warmup draws and draw shape", {
  set.seed(1)
  draws <- draws_list_fixture(n_iterations = 5, n_chains = 1)
  resampled <- upsis:::resample_one_chain(
    draws,
    n_warmup_draws = 2,
    weights = rep(1 / 3, 3)
  )

  expect_type(resampled, "list")
  expect_named(resampled, c("alpha", "sigma"))
  expect_equal(lengths(resampled), c(alpha = 5, sigma = 5))
  expect_equal(resampled$alpha[1:2], c(1, 2))
})

test_that("resampled post-warmup draws are sorted by original draw id", {
  resampled_draws <- posterior::as_draws_df(data.frame(
    alpha = c(4, 1, 4, 2),
    .chain = rep(1, 4),
    .iteration = 1:4,
    .draw = 1:4,
    .upsis_draw_id = c(4, 1, 4, 2)
  ))

  sorted <- upsis:::sort_resampled_draws(resampled_draws)

  expect_equal(sorted$alpha, c(1, 2, 4, 4))
  expect_false(".upsis_draw_id" %in% names(sorted))
  expect_equal(sorted$.draw, 1:4)
  expect_equal(sorted$.iteration, 1:4)
})

test_that("one-chain resampling validates post-warmup weight length", {
  draws <- draws_list_fixture(n_iterations = 5, n_chains = 1)

  expect_error(
    upsis:::resample_one_chain(draws, n_warmup_draws = 2, weights = c(0.5, 0.5)),
    "number of weights"
  )
})

test_that("multi-chain resampling returns one resampled element per chain", {
  set.seed(1)
  draws <- draws_list_fixture(n_iterations = 4, n_chains = 2)
  resampled <- upsis:::resample_by_chain(
    draws,
    n_warmup = 1,
    weights = list(rep(1 / 3, 3), rep(1 / 3, 3))
  )

  expect_length(resampled, 2)
  expect_true(all(vapply(resampled, is.list, logical(1))))
  expect_equal(vapply(resampled, function(x) length(x$alpha), integer(1)), c(4, 4))
})
