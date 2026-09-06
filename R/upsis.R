#' Update a model using PSIS
#'
#' Update a fitted model after adding or removing observations by using
#' Pareto-smoothed importance sampling (PSIS).
#'
#' @param model A model object. Currently only [brms::brmsfit] objects are
#'   supported.
#' @param data_add Additional data collected after the model was fit.
#' @param data_remove Data included in model fitting whose influence should be
#'   removed.
#'
#' @return An object of class `"upsis"` containing the updated model, the raw
#'   PSIS object, the Pareto `k` diagnostic, and the PSIS weights.
#' @export
upsis <- function(model, data_add = NULL, data_remove = NULL) {
  UseMethod("upsis")
}

#' Update a brms model using PSIS
#'
#' `upsis.brmsfit()` computes the change in log likelihood associated with
#' `data_add` and `data_remove`, uses PSIS to smooth the resulting importance
#' weights, and stratified-resamples the post-warmup draws inside the
#' `brmsfit` object's underlying Stan fit. The returned `updated_model` can then
#' be passed to standard `brms` post-processing functions.
#'
#' @inheritParams upsis
#'
#' @section Assumptions:
#' This update is an approximation to refitting the model on the changed data.
#' It is intended for data updates that do not radically change the
#' posterior and do not change the fitted model structure. In particular, it
#' does not update data-dependent priors, spline knots, basis expansions, or
#' other model components chosen when the original model was fit. Updates that
#' introduce unsupported new grouping levels or other new parameters may require
#' a full refit instead.
#'
#' @section Diagnostics:
#' The Pareto `k` diagnostic is returned as `pareto_k` and printed with the
#' result. Values above 0.7 trigger a warning and indicate that the PSIS
#' approximation may be unreliable.
#'
#' @return An object of class `"upsis"` with elements:
#' \describe{
#'   \item{`updated_model`}{A `brmsfit` object with resampled posterior draws.}
#'   \item{`psis`}{The object returned by [loo::psis()].}
#'   \item{`pareto_k`}{The scalar Pareto `k` diagnostic.}
#'   \item{`weights`}{The smoothed importance weights used for resampling.}
#' }
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(upsis)
#'
#' fit <- brm(y ~ x, data = old_data)
#' result <- upsis(fit, data_add = new_data)
#'
#' summary(result$updated_model)
#' result$pareto_k
#' }
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
  draw_weights <- psis_weights(psis_out)
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
    pareto_k = pareto_k(psis_out),
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
  class(out) <- c("upsis", "list")
  out
}

#' @export
print.upsis <- function(x, ...) {
  cat("upsis update\n")
  if (!is.null(x$psis)) {
    cat(sprintf("Pareto k: %s\n", format(x$pareto_k, digits = 3)))
  } else {
    cat("Pareto k: not computed\n")
  }
  invisible(x)
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
  validate_log_ratios(log_ratios)

  message("Performing PSIS.")
  psis_out <- loo::psis(log_ratios)
  k <- pareto_k(psis_out)
  message(sprintf("Pareto k is %s.", format(k, digits = 3)))

  if (k > pareto_k_threshold) {
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

validate_log_ratios <- function(log_ratios) {
  if (!is.numeric(log_ratios) || length(log_ratios) == 0) {
    stop("PSIS log ratios must be a non-empty numeric vector.", call. = FALSE)
  }
  if (anyNA(log_ratios) || any(!is.finite(log_ratios))) {
    stop("PSIS log ratios must all be finite and non-missing.", call. = FALSE)
  }
  invisible(log_ratios)
}

pareto_k <- function(psis_out) {
  k <- psis_out$diagnostics$pareto_k
  if (!is.numeric(k) || length(k) != 1 || is.na(k)) {
    stop("PSIS output did not include a scalar Pareto k diagnostic.", call. = FALSE)
  }
  k
}

psis_weights <- function(psis_out) {
  weights <- loo::weights.importance_sampling(psis_out, log = FALSE)
  if (!is.numeric(weights) || length(weights) == 0) {
    stop("PSIS output did not produce a non-empty numeric weight vector.", call. = FALSE)
  }
  if (anyNA(weights) || any(weights < 0)) {
    stop("PSIS weights must be non-missing and non-negative.", call. = FALSE)
  }
  weights
}

split_weights_by_chain <- function(weights, n_chains) {
  if (length(n_chains) != 1 || n_chains < 1) {
    stop("`n_chains` must be a positive scalar.", call. = FALSE)
  }
  if (length(weights) %% n_chains != 0) {
    stop("The number of PSIS weights must be divisible by the number of chains.", call. = FALSE)
  }

  draws_per_chain <- length(weights) / n_chains
  chain_id <- rep(seq_len(n_chains), each = draws_per_chain)
  unname(split(weights, chain_id))
}

n_warmup_draws_per_chain <- function(draws_list, n_post_warmup_draws) {
  assertthat::assert_that("draws_list" %in% class(draws_list))
  if (length(n_post_warmup_draws) != 1 || n_post_warmup_draws < 0) {
    stop("`n_post_warmup_draws` must be a non-negative scalar.", call. = FALSE)
  }

  n_chains <- posterior::nchains(draws_list)
  if (n_post_warmup_draws %% n_chains != 0) {
    stop("Post-warmup draws must be evenly divisible across chains.", call. = FALSE)
  }

  n_warmup_total <- posterior::ndraws(draws_list) - n_post_warmup_draws
  if (n_warmup_total < 0) {
    stop("The requested post-warmup draw count exceeds the stored draw count.", call. = FALSE)
  }
  if (n_warmup_total %% n_chains != 0) {
    stop("Warmup draws must be evenly divisible across chains.", call. = FALSE)
  }

  n_warmup_total / n_chains
}

resample_by_chain <- function(draws_list, n_warmup, weights) {
  if (length(n_warmup) == 1) {
    n_warmup <- rep(n_warmup, posterior::nchains(draws_list))
  }
  assertthat::assert_that("draws_list" %in% class(draws_list))
  assertthat::assert_that(length(draws_list) == length(weights))
  assertthat::assert_that(length(draws_list) == posterior::nchains(draws_list))
  assertthat::assert_that(length(n_warmup) == length(weights))

  out <- list()
  for (i in seq_along(draws_list)) {
    out[[i]] <- resample_one_chain(draws_list[i], n_warmup[i], weights[[i]])
  }

  out
}

resample_one_chain <- function(draws_list, n_warmup_draws, weights) {
  draws_list <- posterior::as_draws_list(draws_list)
  assertthat::assert_that(posterior::nchains(draws_list) == 1)
  if (length(n_warmup_draws) != 1 || n_warmup_draws < 0) {
    stop("`n_warmup_draws` must be a non-negative scalar.", call. = FALSE)
  }
  n_draws <- posterior::ndraws(draws_list)
  if (n_draws <= n_warmup_draws) {
    stop("Each chain must contain at least one post-warmup draw.", call. = FALSE)
  }
  if ((n_draws - n_warmup_draws) != length(weights)) {
    stop(
      "The number of weights must equal the number of post-warmup draws.",
      call. = FALSE
    )
  }
  draws_df <- posterior::as_draws_df(draws_list)
  post_warmup_draws <- draws_df[(n_warmup_draws + 1):n_draws, ]
  post_warmup_draws$.upsis_draw_id <- seq_len(nrow(post_warmup_draws))
  resampled_draws <- posterior::resample_draws(
    post_warmup_draws,
    weights = weights,
    method = "stratified"
  )
  resampled_draws <- sort_resampled_draws(resampled_draws)
  if (n_warmup_draws > 0) {
    warmup_draws <- draws_df[1:n_warmup_draws, ]
    resampled_draws <- rbind(warmup_draws, resampled_draws)
    resampled_draws$.chain <- 1
    resampled_draws$.draw <- 1:n_draws
    resampled_draws$.iteration <- 1:n_draws
  }
  posterior::as_draws_list(resampled_draws)[[1]]
}

sort_resampled_draws <- function(resampled_draws) {
  if (!".upsis_draw_id" %in% names(resampled_draws)) {
    stop("Resampled draws are missing the original draw id.", call. = FALSE)
  }

  resampled_draws <- resampled_draws[order(resampled_draws$.upsis_draw_id), ]
  resampled_draws$.upsis_draw_id <- NULL
  resampled_draws$.draw <- seq_len(nrow(resampled_draws))
  resampled_draws$.iteration <- seq_len(nrow(resampled_draws))
  resampled_draws
}
