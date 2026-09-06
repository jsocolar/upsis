#' Update a model using PSIS
#' @param model a model object (currently only brmsfit supported)
#' @param data_add additional data collected after the model was fit
#' @param data_remove data included in model fitting whose influence to remove
#' @export
upsis <- function(model, data_add = NULL, data_remove = NULL) {
  UseMethod("upsis")
}

#' Update a brms model using PSIS
#' @param model a brmsfit object
#' @param data_add additional data collected after the model was fit
#' @param data_remove data included in model fitting whose influence to remove
#' @export
#' @method upsis brmsfit
upsis.brmsfit <- function(model, data_add = NULL, data_remove = NULL) {
  if (is.null(data_add) && is.null(data_remove)) {
    warning(
      "No data added or removed; returning a result with the original model.",
      call. = FALSE
    )
    return(new_upsis_result(model))
  }

  ll_change <- brms_log_lik_change(model, data_add, data_remove)
  psis_out <- psis_update(ll_change)
  draw_weights <- loo::weights.importance_sampling(psis_out, log = FALSE)
  draw_weights_list <- split_weights_by_chain(
    draw_weights,
    n_chains = brms::nchains(model)
  )

  updated_model <- model
  draws_list <- posterior::as_draws_list(model$fit@sim$samples)
  n_warmup <- n_warmup_draws_per_chain(
    draws_list,
    n_post_warmup_draws = brms::ndraws(model)
  )

  updated_model$fit@sim$samples <-
    resample_by_chain(
      draws_list, n_warmup, draw_weights_list
    )

  for (i in seq_along(updated_model$fit@sim$samples)) {
    attributes(updated_model$fit@sim$samples[[i]]) <-
      attributes(model$fit@sim$samples[[i]])
  }

  new_upsis_result(
    updated_model = updated_model,
    psis = psis_out,
    pareto_k = psis_out$diagnostics$pareto_k,
    weights = draw_weights
  )
}

new_upsis_result <- function(updated_model, psis = NULL, pareto_k = NA_real_,
                             weights = NULL) {
  out <- list(
    updated_model = updated_model,
    psis = psis,
    pareto_k = pareto_k,
    weights = weights
  )
  class(out) <- "upsis"
  out
}

brms_log_lik_change <- function(model, data_add = NULL, data_remove = NULL) {
  n_draws <- brms::ndraws(model)
  ll_add <- ll_remove <- rep(0, n_draws)

  if (!is.null(data_add)) {
    ll_add <- log_lik_rowsums(model, data_add, "data_add")
  }
  if (!is.null(data_remove)) {
    ll_remove <- log_lik_rowsums(model, data_remove, "data_remove")
  }

  ll_add - ll_remove
}

log_lik_rowsums <- function(model, newdata, data_arg) {
  log_lik <- brms::log_lik(model, newdata = newdata)
  n_draws <- brms::ndraws(model)

  if (length(dim(log_lik)) != 2 || nrow(log_lik) != n_draws) {
    stop(
      sprintf(
        "`brms::log_lik(model, newdata = %s)` must return a matrix with %d rows.",
        data_arg,
        n_draws
      ),
      call. = FALSE
    )
  }

  rowSums(log_lik)
}

psis_update <- function(log_ratios, pareto_k_threshold = 0.7) {
  message("Performing PSIS.")
  psis_out <- loo::psis(log_ratios)
  pareto_k <- psis_out$diagnostics$pareto_k
  message(sprintf("Pareto k is %s.", format(pareto_k, digits = 3)))

  if (pareto_k > pareto_k_threshold) {
    warning(
      sprintf(
        paste0(
          "Pareto k is greater than %s, indicating that the PSIS ",
          "approximation to the updated posterior may be unreliable."
        ),
        pareto_k_threshold
      ),
      call. = FALSE
    )
  }

  psis_out
}

split_weights_by_chain <- function(weights, n_chains) {
  assertthat::assert_that(length(n_chains) == 1)
  assertthat::assert_that(n_chains > 0)
  assertthat::assert_that(length(weights) %% n_chains == 0)

  draws_per_chain <- length(weights) / n_chains
  chain_id <- rep(seq_len(n_chains), each = draws_per_chain)
  unname(split(weights, chain_id))
}

n_warmup_draws_per_chain <- function(draws_list, n_post_warmup_draws) {
  assertthat::assert_that("draws_list" %in% class(draws_list))
  assertthat::assert_that(length(n_post_warmup_draws) == 1)
  assertthat::assert_that(n_post_warmup_draws >= 0)

  n_chains <- posterior::nchains(draws_list)
  assertthat::assert_that(n_post_warmup_draws %% n_chains == 0)

  n_warmup_total <- posterior::ndraws(draws_list) - n_post_warmup_draws
  assertthat::assert_that(n_warmup_total >= 0)
  assertthat::assert_that(n_warmup_total %% n_chains == 0)

  n_warmup_total / n_chains
}

resample_by_chain <- function(draws_list, n_warmup, weights) {
  if (length(n_warmup) == 1) {
    n_warmup <- rep(n_warmup, posterior::nchains(draws_list))
  }
  assertthat::assert_that("draws_list" %in% class(draws_list))
  assertthat::assert_that(length(draws_list) == length(weights))
  assertthat::assert_that(length(draws_list) == posterior::nchains(draws_list))

  out <- list()
  for (i in seq_along(draws_list)) {
    out[[i]] <- resample_one_chain(draws_list[i], n_warmup[i], weights[[i]])
  }

  out
}

resample_one_chain <- function(draws_list, n_warmup_draws, weights) {
  draws_list <- posterior::as_draws_list(draws_list)
  assertthat::assert_that(posterior::nchains(draws_list) == 1)
  assertthat::assert_that(length(n_warmup_draws) == 1)
  n_draws <- posterior::ndraws(draws_list)
  assertthat::assert_that(n_draws > n_warmup_draws)
  assertthat::assert_that((n_draws - n_warmup_draws) == length(weights))
  draws_df <- posterior::as_draws_df(draws_list)
  post_warmup_draws <- draws_df[(n_warmup_draws + 1):n_draws, ]
  resampled_draws <- posterior::resample_draws(
    post_warmup_draws,
    weights = weights,
    method = "stratified"
  )
  if (n_warmup_draws > 0) {
    warmup_draws <- draws_df[1:n_warmup_draws, ]
    resampled_draws <- rbind(warmup_draws, resampled_draws)
    resampled_draws$.chain <- 1
    resampled_draws$.draw <- 1:n_draws
    resampled_draws$.iteration <- 1:n_draws
  }
  posterior::as_draws_list(resampled_draws)[[1]]
}
