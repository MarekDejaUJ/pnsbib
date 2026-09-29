#' Sensitivity to bounded interventional margins
#'
#' For each nonnegative delta, allows every interventional cell probability to
#' vary within the absolute band `[a - delta, a + delta]`, subject to each row
#' remaining a probability distribution and to the fixed observational table.
#' An exact response-type linear program optimizes the requested event jointly
#' over every compatible interventional table in the band. The interval is an
#' identification/sensitivity bound, not a confidence interval. Delta is a
#' hypothetical departure from the supplied interventional margins; it does
#' not by itself quantify an amount of unmeasured confounding.
#'
#' @param model A model returned by [poc_model()].
#' @param query A query returned by [poc_query()].
#' @param delta Nonnegative absolute probability departures, at most one.
#' @param max_variables Maximum number of response-type LP variables.
#' @inheritParams poc_exact
#' @return A data frame with one exact interval per delta.
#' @export
poc_sensitivity <- function(model, query, delta = c(0, 0.025, 0.05, 0.1),
                            max_variables = 10000L, lp_backend = NULL) {
  lp_backend <- .poc_lp_backend(lp_backend, model)
  if (!inherits(model, "pnsbib_model")) stop("Expected a pnsbib_model.", call. = FALSE)
  if (!is.numeric(delta) || !length(delta) || anyNA(delta) ||
      any(!is.finite(delta)) || any(delta < 0 | delta > 1)) {
    stop("delta must contain finite values in [0, 1].", call. = FALSE)
  }
  idx <- .poc_indices(model, query)
  denominator <- .poc_denominator(model, query, idx)
  if (query$conditional && denominator <= model$tolerance) {
    return(data.frame(delta = delta, lower = NA_real_, upper = NA_real_,
                      status = "undefined_condition", denominator = denominator,
                      assumption = "absolute_intervention_probability_band"))
  }
  system <- .poc_lp_system(model, max_variables = max_variables, lp_backend = lp_backend)
  objective <- .poc_event_indicator(system, idx, lp_backend)
  observed_rows <- seq_len(1L + length(model$o))
  intervention_rows <- seq.int(max(observed_rows) + 1L, nrow(system$A))
  constraints <- rbind(system$A[observed_rows, , drop = FALSE],
                       system$A[intervention_rows, , drop = FALSE],
                       system$A[intervention_rows, , drop = FALSE])
  directions <- c(rep("=", length(observed_rows)),
                  rep(">=", length(intervention_rows)),
                  rep("<=", length(intervention_rows)))
  answers <- lapply(delta, function(d) {
    lower_a <- pmax(model$a - d, 0)
    upper_a <- pmin(model$a + d, 1)
    rhs <- c(system$b[observed_rows], as.vector(t(lower_a)),
             as.vector(t(upper_a)))
    lo <- .poc_lp("min", as.numeric(objective), constraints,
                   directions, rhs, lp_backend)
    hi <- .poc_lp("max", as.numeric(objective), constraints,
                   directions, rhs, lp_backend)
    status <- if (lo$status == 0L && hi$status == 0L) "ok" else
      "incompatible_margin_band"
    data.frame(delta = d,
               lower = if (status == "ok") lo$objval / denominator else NA_real_,
               upper = if (status == "ok") hi$objval / denominator else NA_real_,
               status = status, denominator = denominator,
               assumption = "absolute_intervention_probability_band")
  })
  do.call(rbind, answers)
}

#' Sensitivity to simultaneous observational and interventional cell bands
#'
#' Allows each observational joint cell and each interventional outcome cell
#' to vary by a specified absolute probability amount around the supplied
#' model. The response-type LP preserves normalization and consistency while
#' optimizing jointly over all compatible input tables. Conditional queries
#' use a linear-fractional transformation and range over compatible tables
#' where the factual condition has positive probability. A positive
#' `minimum_condition_probability` can exclude nearly empty factual strata.
#' Cell bands are hypothetical input departures; they are not error rates,
#' confidence intervals or an identified confounding model.
#'
#' @param model A model returned by [poc_model()].
#' @param query A query returned by [poc_query()].
#' @param observed_delta Scalar absolute band around every observational cell.
#' @param interventional_delta Scalar absolute band around every intervention cell.
#' @param minimum_condition_probability Optional positive lower bound on the
#'   factual condition probability for conditional queries.
#' @param max_variables Maximum number of response-type LP variables.
#' @inheritParams poc_exact
#' @return A one-row data frame with sharp finite-model bounds under the bands.
#' @export
poc_input_sensitivity <- function(model, query, observed_delta = 0,
                                  interventional_delta = 0,
                                  minimum_condition_probability = NULL,
                                  max_variables = 10000L, lp_backend = NULL) {
  lp_backend <- .poc_lp_backend(lp_backend, model)
  if (!inherits(model, "pnsbib_model")) {
    stop("Expected a pnsbib_model.", call. = FALSE)
  }
  check_delta <- function(value, label) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 0 || value > 1) {
      stop(sprintf("%s must be one finite value in [0, 1].", label),
           call. = FALSE)
    }
  }
  check_delta(observed_delta, "observed_delta")
  check_delta(interventional_delta, "interventional_delta")
  if (!is.null(minimum_condition_probability) &&
      (!is.numeric(minimum_condition_probability) ||
       length(minimum_condition_probability) != 1L ||
       is.na(minimum_condition_probability) ||
       !is.finite(minimum_condition_probability) ||
       minimum_condition_probability <= 0 ||
       minimum_condition_probability > 1)) {
    stop("minimum_condition_probability must be in (0, 1].", call. = FALSE)
  }
  idx <- .poc_indices(model, query)
  if (!query$conditional && !is.null(minimum_condition_probability)) {
    stop("minimum_condition_probability requires a conditional query.",
         call. = FALSE)
  }
  system <- .poc_lp_system(model, max_variables = max_variables, lp_backend = lp_backend)
  ncells <- length(model$o)
  obs_rows <- seq.int(2L, 1L + ncells)
  int_rows <- seq.int(2L + ncells, 1L + 2L * ncells)
  Aobs <- system$A[obs_rows, , drop = FALSE]
  Aint <- system$A[int_rows, , drop = FALSE]
  lo_obs <- pmax(as.vector(t(model$o)) - observed_delta, 0)
  hi_obs <- pmin(as.vector(t(model$o)) + observed_delta, 1)
  lo_int <- pmax(as.vector(t(model$a)) - interventional_delta, 0)
  hi_int <- pmin(as.vector(t(model$a)) + interventional_delta, 1)
  A <- rbind(system$A[1L, , drop = FALSE], Aobs, Aobs, Aint, Aint)
  directions <- c("=", rep(">=", ncells), rep("<=", ncells),
                  rep(">=", ncells), rep("<=", ncells))
  rhs <- c(1, lo_obs, hi_obs, lo_int, hi_int)
  event <- .poc_event_indicator(system, idx, lp_backend)
  fit_feasible <- .poc_lp("min", rep(0, system$nvar), A,
                           directions, rhs, lp_backend)
  result <- function(lower, upper, status, condition_min = NA_real_,
                     condition_max = NA_real_) {
    data.frame(
      observed_delta = observed_delta,
      interventional_delta = interventional_delta,
      lower = lower, upper = upper, status = status,
      method = "exact_response_type_lp",
      condition_probability_min = condition_min,
      condition_probability_max = condition_max,
      minimum_condition_probability = if (is.null(minimum_condition_probability))
        NA_real_ else minimum_condition_probability,
      assumption = "joint_absolute_input_cell_bands")
  }
  if (fit_feasible$status != 0L) {
    return(result(NA_real_, NA_real_, "incompatible_margin_bands"))
  }
  if (!query$conditional) {
    lo <- .poc_lp("min", as.numeric(event), A, directions, rhs, lp_backend)
    hi <- .poc_lp("max", as.numeric(event), A, directions, rhs, lp_backend)
    if (lo$status != 0L || hi$status != 0L) {
      return(result(NA_real_, NA_real_, "solver_failure"))
    }
    return(result(lo$objval, hi$objval, "ok"))
  }

  condition <- idx
  condition$x <- condition$y <- integer()
  factual <- .poc_event_indicator(system, condition, lp_backend)
  condition_lo <- .poc_lp("min", as.numeric(factual), A, directions, rhs, lp_backend)
  condition_hi <- .poc_lp("max", as.numeric(factual), A, directions, rhs, lp_backend)
  if (condition_lo$status != 0L || condition_hi$status != 0L) {
    return(result(NA_real_, NA_real_, "solver_failure"))
  }
  pmin_condition <- condition_lo$objval
  pmax_condition <- condition_hi$objval
  if (pmax_condition <= model$tolerance ||
      (!is.null(minimum_condition_probability) &&
       pmax_condition + model$tolerance < minimum_condition_probability)) {
    return(result(NA_real_, NA_real_, "undefined_condition",
                  pmin_condition, pmax_condition))
  }

  # Charnes-Cooper: z = q / P(factual), t = 1 / P(factual).
  # The augmented LP has nonnegative (z, t), factual(z) = 1,
  # sum(z) = t, and each input band multiplied by t.
  nvar <- system$nvar
  augmented <- rbind(
    c(rep(1, nvar), -1),
    c(as.numeric(factual), 0),
    cbind(Aobs, -lo_obs), cbind(Aobs, -hi_obs),
    cbind(Aint, -lo_int), cbind(Aint, -hi_int))
  aug_directions <- c("=", "=", rep(">=", ncells), rep("<=", ncells),
                      rep(">=", ncells), rep("<=", ncells))
  aug_rhs <- c(0, 1, rep(0, 4L * ncells))
  if (!is.null(minimum_condition_probability)) {
    augmented <- rbind(augmented, c(rep(0, nvar), 1))
    aug_directions <- c(aug_directions, "<=")
    aug_rhs <- c(aug_rhs, 1 / minimum_condition_probability)
  }
  objective <- c(as.numeric(event), 0)
  lo <- .poc_lp("min", objective, augmented, aug_directions, aug_rhs, lp_backend)
  hi <- .poc_lp("max", objective, augmented, aug_directions, aug_rhs, lp_backend)
  if (lo$status != 0L || hi$status != 0L) {
    return(result(NA_real_, NA_real_, "solver_failure",
                  pmin_condition, pmax_condition))
  }
  result(lo$objval, hi$objval, "ok", pmin_condition, pmax_condition)
}
