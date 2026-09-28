.poc_lp_system <- function(model, max_variables = 10000L, lp_backend = NULL) {
  lp_backend <- .poc_lp_backend(lp_backend, model)
  if (lp_backend == "reference") return(.poc_lp_system_reference(model, max_variables))
  if (!is.numeric(max_variables) || length(max_variables) != 1L ||
      !is.finite(max_variables) || max_variables < 1)
    stop("Invalid max_variables.", call. = FALSE)
  m <- nrow(model$o)
  z <- poc_static_system_zig(m, ncol(model$o), as.numeric(model$o), as.numeric(model$a),
    as.integer(min(floor(max_variables), .Machine$integer.max)))
  nvar <- z[1L]; rows <- z[2L]
  cursor <- 3L
  take <- function(n) {
    values <- z[cursor + seq_len(n)]
    cursor <<- cursor + n
    values
  }
  response <- matrix(as.integer(take(nvar * m)), nvar, m,
    dimnames = list(NULL, paste0("Var", seq_len(m))))
  factual_x <- as.integer(take(nvar))
  factual_y <- as.integer(take(nvar))
  A <- matrix(take(rows * nvar), rows, nvar)
  b <- take(rows)
  list(A = A, b = b, response = response, factual_x = factual_x,
    factual_y = factual_y, nvar = nvar)
}

.poc_lp_system_reference <- function(model, max_variables = 10000L) {
  m <- nrow(model$o)
  n <- ncol(model$o)
  ntypes <- n^m
  nvar <- m * ntypes
  if (!is.finite(nvar) || nvar > max_variables) {
    stop(sprintf("The exact LP needs %.0f variables, above the configured limit.", nvar),
         call. = FALSE)
  }
  types <- as.matrix(do.call(expand.grid, rep(list(seq_len(n)), m)))
  storage.mode(types) <- "integer"
  response <- types[rep(seq_len(ntypes), times = m), , drop = FALSE]
  factual_x <- rep(seq_len(m), each = ntypes)
  factual_y <- response[cbind(seq_len(nvar), factual_x)]
  A <- matrix(0, 1L + 2L * m * n, nvar)
  b <- numeric(nrow(A))
  A[1L, ] <- 1
  b[1L] <- 1
  row <- 1L
  for (x in seq_len(m)) for (y in seq_len(n)) {
    row <- row + 1L
    A[row, ] <- as.numeric(factual_x == x & factual_y == y)
    b[row] <- model$o[x, y]
  }
  for (x in seq_len(m)) for (y in seq_len(n)) {
    row <- row + 1L
    A[row, ] <- as.numeric(response[, x] == y)
    b[row] <- model$a[x, y]
  }
  list(A = A, b = b, response = response, factual_x = factual_x,
       factual_y = factual_y, nvar = nvar)
}

#' Exact finite response-type probability bounds
#'
#' Solves two linear programs over the complete distribution of observed
#' treatment and potential-outcome response types. Bounds are sharp for the
#' supplied finite-model observational and interventional margins.
#'
#' @param model A model returned by [poc_model()].
#' @param query A query returned by [poc_query()].
#' @param max_variables Maximum LP variables before declining exact solution.
#' @param lp_backend NULL inherits the model's LP solver; "zig" or "reference"
#'   explicitly overrides it. Older models use "reference". No fallback.
#' @return A one-row data frame with lower, upper, status and method.
#' @export
poc_exact <- function(model, query, max_variables = 10000L, lp_backend = NULL) {
  lp_backend <- .poc_lp_backend(lp_backend, model)
  idx <- .poc_indices(model, query)
  denominator <- .poc_denominator(model, query, idx)
  if (denominator <= model$tolerance && query$conditional) {
    return(.poc_result(NA_real_, NA_real_, "exact_lp", "sharp",
                       "undefined_condition", denominator, query))
  }
  system <- .poc_lp_system(model, max_variables, lp_backend)
  obj <- .poc_event_indicator(system, idx, lp_backend)
  directions <- rep("=", nrow(system$A))
  fit_lo <- .poc_lp("min", as.numeric(obj), system$A,
                     directions, system$b, lp_backend)
  fit_hi <- .poc_lp("max", as.numeric(obj), system$A,
                     directions, system$b, lp_backend)
  if (fit_lo$status != 0L || fit_hi$status != 0L) {
    return(.poc_result(NA_real_, NA_real_, "exact_lp", "sharp",
                       "incompatible_margins", denominator, query))
  }
  .poc_result(fit_lo$objval / denominator, fit_hi$objval / denominator,
              "exact_lp", "sharp", "ok", denominator, query)
}
