draws_list_fixture <- function(n_iterations = 5, n_chains = 2) {
  values <- seq_len(n_iterations * n_chains * 2)
  draws <- array(values, dim = c(n_iterations, n_chains, 2))
  dimnames(draws) <- list(
    iteration = NULL,
    chain = NULL,
    variable = c("alpha", "sigma")
  )
  posterior::as_draws_list(posterior::as_draws_array(draws))
}
