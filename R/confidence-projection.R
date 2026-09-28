#' Project a simultaneous margin confidence region through the finite LP
#'
#' Optimizes over every latent distribution compatible with the confidence
#' cells and structural template. All numerical template probabilities are
#' ignored. A common input-region coverage event covers the entire identified
#' sets of queries on this same fixed model. This guarantee requires the
#' region's sampling assumptions and correct causal/structural restrictions.
#'
#' Conditional targets range over positive conditioning probability. A zero
#' lower denominator is reported, not removed by a data-dependent floor.
#' Empty regions and solver failures return explicit statuses and NA endpoints;
#' do not exclude those samples when assessing coverage. A population target
#' with zero conditioning probability is undefined.
#'
#' @param region A [poc_confidence_region()].
#' @param query A compatible static or longitudinal query.
#' @param max_variables Maximum LP variables for a static template; longitudinal
#'   templates retain their own declared cap.
#' @param witnesses Attach normalized latent endpoint distributions for audit.
#' @param lp_backend NULL inherits the template's LP solver; "zig" or
#'   "reference" overrides it. Older templates use "reference". No fallback.
#' @return A one-row static or longitudinal result with confidence metadata.
#'   sharpness describes projection of the supplied region, not an optimal
#'   statistical confidence procedure. Optional witnesses are attributes.
#' @export
poc_confidence_bounds <- function(region, query, max_variables = 10000L,
                                  witnesses = FALSE, lp_backend = NULL) {
  if (!inherits(region, "pnsbib_confidence_region") || !identical(region$assume_sampling, TRUE))
    stop("Expected a declared pnsbib_confidence_region.", call. = FALSE)
  if (!is.logical(witnesses) || length(witnesses) != 1L || is.na(witnesses) ||
      !is.numeric(max_variables) || length(max_variables) != 1L ||
      !is.finite(max_variables) || max_variables < 1)
    stop("Invalid witness flag or variable cap.", call. = FALSE)
  model <- region$template
  lp_backend <- .poc_lp_backend(lp_backend, model)
  map <- poc_margin_map(model)
  cells <- region$cells
  if (!is.data.frame(cells) || !all(c(names(map), "lower", "upper") %in% names(cells)) ||
      !identical(cells[, names(map), drop = FALSE], map) ||
      !is.numeric(cells$lower) || !is.numeric(cells$upper) ||
      any(!is.finite(cells$lower)) || any(!is.finite(cells$upper)) ||
      any(cells$lower < 0 | cells$upper > 1 | cells$lower > cells$upper))
    stop("Confidence cells no longer match the structural template.", call. = FALSE)
  longitudinal <- inherits(model, "pnsbib_longitudinal_model")
  if (longitudinal) {
    event <- .poc_long_event(model, query)
    system <- event$system
    objective <- event$event
    condition <- event$condition
    conditional <- event$conditional
  } else {
    idx <- .poc_indices(model, query)
    system <- .poc_lp_system(model, max_variables, lp_backend)
    objective <- .poc_event_indicator(system, idx, lp_backend)
    condition_indices <- idx
    condition_indices$x <- condition_indices$y <- integer()
    condition <- .poc_event_indicator(system, condition_indices, lp_backend)
    conditional <- query$conditional
  }
  stopifnot(nrow(system$A) == nrow(cells) + 1L)
  margin_A <- system$A[-1L, , drop = FALSE]
  A <- rbind(rep(1, system$nvar), margin_A, margin_A)
  directions <- c("=", rep(">=", nrow(cells)), rep("<=", nrow(cells)))
  rhs <- c(1, cells$lower, cells$upper)
  solve <- function(sense, objective) .poc_lp(sense, objective, A, directions, rhs, lp_backend)
  result <- function(lower, upper, status, cmin = NA_real_, cmax = NA_real_, endpoints = NULL) {
    denominator <- if (!conditional) 1 else
      if (is.finite(cmin) && is.finite(cmax) && abs(cmax - cmin) <= model$tolerance)
        (cmin + cmax) / 2 else NA_real_
    if (longitudinal) {
      out <- structure(data.frame(estimand = toupper(query$kind), event = event$label,
        regime = query$regime, comparison = query$reference, horizon = query$horizon,
        lower = lower, upper = upper, denominator = denominator,
        condition_on = event$condition_on, condition_horizon = event$condition_horizon,
        regime_outcome = event$regime_outcome, reference_outcome = event$reference_outcome,
        outcome_order = event$outcome_order, method = "confidence_region_projection",
        sharpness = "sharp_projection_of_region", status = status),
        class = c("pnsbib_longitudinal_result", "data.frame"))
    } else out <- .poc_result(lower, upper, "confidence_region_projection",
      "sharp_projection_of_region", status, denominator, query)
    out$conf_level <- region$conf_level
    out$region_method <- paste0(region$method, "_bonferroni")
    out$family_cells <- region$family_cells
    out$sampling_unit <- region$sampling_unit
    out$target <- region$target
    out$uncertainty_kind <- "simultaneous_identified_set_confidence_region"
    out$condition_probability_min <- if (conditional) cmin else 1
    out$condition_probability_max <- if (conditional) cmax else 1
    out$condition_may_be_zero <- if (!conditional) FALSE else
      if (!is.finite(cmin)) NA else cmin <= model$tolerance
    out$assumptions <- paste0(region$assumptions, if (longitudinal)
      paste0(";absorbing=", model$absorbing, ";no_anticipation=", model$no_anticipation) else "")
    if (witnesses && !is.null(endpoints)) attr(out, "witnesses") <- endpoints
    out
  }
  feasible <- solve("min", rep(0, system$nvar))
  if (feasible$status != 0L) return(result(NA_real_, NA_real_,
    if (feasible$status == 2L) "empty_confidence_region" else "solver_failure"))
  cmin <- cmax <- 1
  if (conditional) {
    low_c <- solve("min", condition)
    high_c <- solve("max", condition)
    if (low_c$status != 0L || high_c$status != 0L)
      return(result(NA_real_, NA_real_, "solver_failure"))
    cmin <- low_c$objval; cmax <- high_c$objval
    if (cmax <= model$tolerance)
      return(result(NA_real_, NA_real_, "undefined_condition", cmin, cmax))
    # z = q / d(q), s = 1 / d(q). The normalization row becomes sum(z)=s.
    A <- rbind(cbind(A, -rhs), c(condition, 0))
    directions <- c(directions, "=")
    rhs <- c(rep(0, length(rhs)), 1)
    objective <- c(objective, 0)
  }
  lo <- solve("min", objective)
  hi <- solve("max", objective)
  if (lo$status != 0L || hi$status != 0L)
    return(result(NA_real_, NA_real_, "solver_failure", cmin, cmax))
  unpack <- function(fit) {
    q <- fit$solution[seq_len(system$nvar)]
    if (conditional) q <- q / fit$solution[system$nvar + 1L]
    list(probability = q, response = system$response,
         factual = if (longitudinal) system$factual else system$factual_x)
  }
  result(lo$objval, hi$objval, "ok", cmin, cmax,
    if (witnesses) list(lower = unpack(lo), upper = unpack(hi)) else NULL)
}
