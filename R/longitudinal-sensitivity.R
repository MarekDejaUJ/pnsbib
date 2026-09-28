#' Sensitivity to longitudinal observed and intervention margin bands
#'
#' Relax each supplied observed history/path cell, intervention path cell,
#' and intervention horizon/state probability by its absolute probability band.
#' One latent response-type distribution preserves normalization, consistency,
#' allowed paths, absorption and no anticipation jointly. Missing intervention
#' margins remain free. Bands are hypothetical departures, not confidence
#' intervals, source-label error rates or an identified confounding model.
#'
#' PN and PS condition on the complete observed history or its prefix, as
#' specified in the query. Prefix conditions pool all matching continuations.
#' With uncertain observed cells, their denominator also varies. A
#' Charnes-Cooper linear-fractional transformation optimizes the ratio over
#' distributions with positive conditioning probability. The optional floor
#' excludes nearly empty conditions. Reported condition probability ranges
#' are computed before imposing this floor. The model supplied as the band
#' center must itself be feasible; this function does not repair invalid input.
#'
#' @param model A model returned by [poc_longitudinal_model()].
#' @param query A query returned by [poc_longitudinal_query()].
#' @param observed_delta Scalar band width for observed joint history/path cells.
#' @param path_delta Scalar band width for supplied interventional path cells.
#' @param time_delta Scalar band width for supplied intervention horizon/state
#'   probability cells, including the legacy binary prob1 schema.
#' @param minimum_condition_probability Optional positive conditioning floor
#'   for PN or PS, at most one.
#' @return A one-row pnsbib_longitudinal_result. The denominator is NA when
#'   it varies over the bands; condition_probability_min/max report its range.
#' @md
#' @examples
#' regimes <- matrix(0:1, ncol = 1, dimnames = list(c("0", "1"), "t1"))
#' observed <- data.frame(history = c("0", "0", "1", "1"),
#'   path = c("0", "1", "0", "1"), prob = c(.3, .2, .1, .4))
#' margins <- data.frame(regime = observed$history, path = observed$path,
#'   prob = c(.7, .3, .4, .6))
#' model <- poc_longitudinal_model(observed, regimes, margins)
#' query <- poc_longitudinal_query("pn", "1", "0", 1)
#' poc_longitudinal_sensitivity(model, query, observed_delta = .05,
#'   path_delta = .1, minimum_condition_probability = .1)
#' @export
#' @inheritParams poc_exact
poc_longitudinal_sensitivity <- function(model, query, observed_delta = 0,
                                         path_delta = 0, time_delta = 0,
                                         minimum_condition_probability = NULL,
                                         lp_backend = NULL) {
  lp_backend <- .poc_lp_backend(lp_backend, model)
  check_delta <- function(value, label) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 0 || value > 1)
      stop(sprintf("%s must be one finite value in [0, 1].", label), call. = FALSE)
  }
  check_delta(observed_delta, "observed_delta")
  check_delta(path_delta, "path_delta")
  check_delta(time_delta, "time_delta")
  if (!is.null(minimum_condition_probability) &&
      (!is.numeric(minimum_condition_probability) ||
       length(minimum_condition_probability) != 1L ||
       is.na(minimum_condition_probability) ||
       !is.finite(minimum_condition_probability) ||
       minimum_condition_probability <= 0 || minimum_condition_probability > 1))
    stop("minimum_condition_probability must be in (0, 1].", call. = FALSE)
  specification <- .poc_long_event(model, query)
  if (!specification$conditional && !is.null(minimum_condition_probability))
    stop("minimum_condition_probability requires PN or PS.", call. = FALSE)
  system <- specification$system
  width <- c(rep(observed_delta, length(model$observed)),
             rep(path_delta, nrow(model$interventional_paths)),
             rep(time_delta, nrow(model$time_margins)))
  cells <- system$A[-1L, , drop = FALSE]
  lower <- pmax(system$b[-1L] - width, 0)
  upper <- pmin(system$b[-1L] + width, 1)
  A <- rbind(system$A[1L, , drop = FALSE], cells, cells)
  direction <- c("=", rep(">=", nrow(cells)), rep("<=", nrow(cells)))
  rhs <- c(1, lower, upper)
  solve <- function(sense, objective, matrix = A, directions = direction,
                    right = rhs) {
    .poc_lp(sense, objective, matrix, directions, right, lp_backend)
  }
  result <- function(lo, hi, status, condition_min = NA_real_,
                     condition_max = NA_real_) {
    denominator <- if (!specification$conditional) 1 else
      if (is.finite(condition_min) && is.finite(condition_max) &&
          abs(condition_max - condition_min) <= model$tolerance)
        (condition_min + condition_max) / 2 else NA_real_
    structure(data.frame(
      estimand = toupper(query$kind), event = specification$label,
      regime = query$regime, comparison = query$reference, horizon = query$horizon,
      lower = lo, upper = hi, denominator = denominator,
      condition_on = specification$condition_on,
      condition_horizon = specification$condition_horizon,
      regime_outcome = specification$regime_outcome,
      reference_outcome = specification$reference_outcome,
      outcome_order = specification$outcome_order,
      method = "longitudinal_exact_lp", sharpness = "sharp_given_bands",
      status = status, observed_delta = observed_delta, path_delta = path_delta,
      time_delta = time_delta, condition_probability_min = condition_min,
      condition_probability_max = condition_max,
      minimum_condition_probability = if (is.null(minimum_condition_probability))
        NA_real_ else minimum_condition_probability,
      assumptions = paste0("joint_absolute_longitudinal_cell_bands;consistency;absorbing=",
        model$absorbing, ";no_anticipation=", model$no_anticipation),
      stringsAsFactors = FALSE), class = c("pnsbib_longitudinal_result", "data.frame"))
  }
  feasible <- solve("min", rep(0, system$nvar))
  if (feasible$status != 0L)
    return(result(NA_real_, NA_real_, if (feasible$status == 2L)
      "incompatible_margin_bands" else "solver_failure"))
  condition_min <- condition_max <- NA_real_
  objective <- specification$event
  if (specification$conditional) {
    condition_lo <- solve("min", specification$condition)
    condition_hi <- solve("max", specification$condition)
    if (condition_lo$status != 0L || condition_hi$status != 0L)
      return(result(NA_real_, NA_real_, "solver_failure"))
    condition_min <- condition_lo$objval
    condition_max <- condition_hi$objval
    if (condition_max <= model$tolerance ||
        (!is.null(minimum_condition_probability) &&
         condition_max + model$tolerance < minimum_condition_probability))
      return(result(NA_real_, NA_real_, "undefined_condition", condition_min, condition_max))
    # z = q / d(q), s = 1 / d(q). Every original row is scaled by s;
    # the normalization row becomes sum(z) = s. d(z) = 1 ensures s > 0.
    A <- rbind(cbind(A, -rhs), c(specification$condition, 0))
    direction <- c(direction, "=")
    rhs <- c(rep(0, length(rhs)), 1)
    if (!is.null(minimum_condition_probability)) {
      A <- rbind(A, c(rep(0, system$nvar), 1))
      direction <- c(direction, "<=")
      rhs <- c(rhs, 1 / minimum_condition_probability)
    }
    objective <- c(objective, 0)
  }
  lo <- solve("min", objective)
  hi <- solve("max", objective)
  if (lo$status != 0L || hi$status != 0L)
    return(result(NA_real_, NA_real_, "solver_failure", condition_min, condition_max))
  result(max(0, lo$objval), min(1, hi$objval), "ok", condition_min, condition_max)
}
