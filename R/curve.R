#' Tabulate longitudinal bounds over a declared horizon/query grid
#'
#' Every row is the result of an existing fixed-regime query on the same model.
#' The helper preserves requested order and undefined statuses; it does not
#' enforce monotone bounds, identify margins or supply simultaneous confidence
#' intervals. Horizons are model-period positions, not automatically time since
#' indexing. Exact path-pair events use [poc_longitudinal_query()] separately.
#'
#' @param model A [poc_longitudinal_model()].
#' @param regime,reference Active and comparison regime names.
#' @param horizons Unique positive integer model-period positions.
#' @param kind Unique query kinds; default pns, pn, ps. Also accepts persistent,
#'   event_time and lagged. Every requested kind/horizon pair must be valid.
#' @param condition_on PN/PS factual conditioning scope: history or prefix.
#' @param regime_outcome,reference_outcome Character outcome sets, as in
#'   [poc_longitudinal_query()].
#' @param start,lag Persistent starting period and lagged-query lag.
#' @param bands NULL for exact inputs, or a named list of arguments to
#'   [poc_longitudinal_sensitivity()]: observed_delta, path_delta, time_delta,
#'   minimum_condition_probability. The conditioning floor applies only to
#'   PN/PS rows and requires at least one such query in the grid.
#' @inheritParams poc_exact
#' @return A multirow pnsbib_longitudinal_result, ordered by the supplied horizon
#'   vector and then kind vector, with a curve_row index. Compatible with
#'   [poc_analysis()], summary and [poc_export_summary()].
#' @examples
#' regimes <- rbind(off = c(0, 0), on = c(1, 1))
#' observed <- expand.grid(history = c("off", "on"), path = c("00", "01", "11"),
#'                         stringsAsFactors = FALSE)
#' observed$prob <- c(.2, .1, .1, .1, .2, .3)
#' model <- poc_longitudinal_model(observed, regimes, absorbing = TRUE)
#' poc_curve(model, "on", "off")
#' poc_curve(model, "on", "off", kind = "pns", bands = list(observed_delta = .05))
#' @export
poc_curve <- function(model, regime, reference,
                       horizons = seq_len(ncol(model$regimes)),
                       kind = c("pns", "pn", "ps"),
                       condition_on = c("history", "prefix"),
                       regime_outcome = "1", reference_outcome = "0",
                       start = 1L, lag = 1L, bands = NULL, lp_backend = NULL) {
  lp_backend <- .poc_lp_backend(lp_backend, model)
  if (!inherits(model, "pnsbib_longitudinal_model"))
    stop("Expected a longitudinal model.", call. = FALSE)
  condition_on <- match.arg(condition_on)
  if (!is.numeric(horizons) || !length(horizons) || any(!is.finite(horizons)) ||
      any(horizons < 1 | horizons > ncol(model$regimes) | horizons != floor(horizons)) ||
      anyDuplicated(horizons) || !.poc_long_labels(kind) ||
      any(!kind %in% c("pns", "pn", "ps", "persistent", "event_time", "lagged")))
    stop("Invalid curve horizons or kinds.", call. = FALSE)
  if (!is.null(bands)) {
    permitted <- c("observed_delta", "path_delta", "time_delta", "minimum_condition_probability")
    if (!is.list(bands) || is.data.frame(bands) ||
        (length(bands) && (!.poc_long_labels(names(bands)) || any(!names(bands) %in% permitted))))
      stop("bands must be a named list of longitudinal sensitivity settings.", call. = FALSE)
    if (!is.null(bands$minimum_condition_probability) && !any(kind %in% c("pn", "ps")))
      stop("A conditioning floor requires PN or PS in the curve.", call. = FALSE)
  }
  queries <- unlist(lapply(horizons, function(t) lapply(kind, function(k)
    poc_longitudinal_query(k, regime, reference, t, start, lag,
      condition_on = if (k %in% c("pn", "ps")) condition_on else "history",
      regime_outcome = regime_outcome, reference_outcome = reference_outcome))), recursive = FALSE)
  result <- lapply(queries, function(q) {
    if (is.null(bands)) return(poc_longitudinal_bounds(model, q, lp_backend = lp_backend))
    settings <- bands
    if (!q$kind %in% c("pn", "ps")) settings$minimum_condition_probability <- NULL
    do.call(poc_longitudinal_sensitivity,
      c(list(model = model, query = q, lp_backend = lp_backend), settings))
  })
  out <- do.call(rbind, result)
  rownames(out) <- NULL
  out$curve_row <- seq_len(nrow(out))
  out
}
