#' Specify a probability-of-causation query
#'
#' @param counterfactual Named character vector mapping treatment labels to
#'   required potential-outcome labels.
#' @param observed_x Optional observed treatment label.
#' @param observed_y Optional observed outcome label.
#' @param conditional Divide by the probability of the specified factual
#'   condition. A zero-probability condition yields an undefined result.
#' @return A pnsbib_query object.
#' @export
poc_query <- function(counterfactual, observed_x = NULL, observed_y = NULL,
                      conditional = FALSE) {
  if (!is.character(counterfactual) || length(counterfactual) < 1L ||
      is.null(names(counterfactual)) || anyNA(counterfactual) ||
      anyNA(names(counterfactual)) ||
      any(!nzchar(counterfactual)) || any(!nzchar(names(counterfactual))) ||
      anyDuplicated(names(counterfactual))) {
    stop("counterfactual must be a named character vector with distinct treatments.", call. = FALSE)
  }
  valid_label <- function(z) is.null(z) ||
    (is.character(z) && length(z) == 1L && !is.na(z) && nzchar(z))
  if (!valid_label(observed_x) || !valid_label(observed_y) ||
      !is.logical(conditional) || length(conditional) != 1L ||
      is.na(conditional) ||
      (conditional && is.null(observed_x) && is.null(observed_y))) {
    stop("Invalid factual condition or conditional flag.", call. = FALSE)
  }
  structure(list(counterfactual = counterfactual, observed_x = observed_x,
                 observed_y = observed_y, conditional = conditional),
            class = "pnsbib_query")
}

.poc_indices <- function(model, query) {
  if (!inherits(model, "pnsbib_model") || !inherits(query, "pnsbib_query")) {
    stop("Expected a pnsbib_model and pnsbib_query.", call. = FALSE)
  }
  .poc_query_indices(rownames(model$o), colnames(model$o), query)
}

.poc_query_indices <- function(treatments, outcomes, query) {
  if (!inherits(query, "pnsbib_query"))
    stop("Expected a pnsbib_query.", call. = FALSE)
  x <- match(names(query$counterfactual), treatments)
  y <- match(unname(query$counterfactual), outcomes)
  ox <- if (is.null(query$observed_x)) NA_integer_ else
    match(query$observed_x, treatments)
  oy <- if (is.null(query$observed_y)) NA_integer_ else
    match(query$observed_y, outcomes)
  if (anyNA(x) || anyNA(y) || (is.na(ox) && !is.null(query$observed_x)) ||
      (is.na(oy) && !is.null(query$observed_y))) {
    stop("The query contains a treatment or outcome absent from the model.", call. = FALSE)
  }
  list(x = x, y = y, ox = ox, oy = oy)
}

.poc_event_indicator <- function(system, indices, lp_backend = "reference") {
  if (.poc_lp_backend(lp_backend) == "reference")
    return(.poc_event_indicator_reference(system, indices))
  poc_static_event_zig(system$nvar, ncol(system$response), as.integer(system$response),
    as.integer(system$factual_x), as.integer(system$factual_y),
    as.integer(indices$x), as.integer(indices$y),
    if (is.na(indices$ox)) 0L else as.integer(indices$ox),
    if (is.na(indices$oy)) 0L else as.integer(indices$oy))
}

.poc_event_indicator_reference <- function(system, indices) {
  event <- rep(TRUE, system$nvar)
  for (i in seq_along(indices$x))
    event <- event & system$response[, indices$x[i]] == indices$y[i]
  if (!is.na(indices$ox)) event <- event & system$factual_x == indices$ox
  if (!is.na(indices$oy)) event <- event & system$factual_y == indices$oy
  as.numeric(event)
}

.poc_denominator <- function(model, q, indices) {
  if (!q$conditional) return(1)
  if (!is.na(indices$ox) && !is.na(indices$oy)) {
    return(model$o[indices$ox, indices$oy])
  }
  if (!is.na(indices$ox)) return(model$px[indices$ox])
  model$py[indices$oy]
}

.poc_theorem <- function(indices) {
  k <- length(indices$x)
  if (k == 1L) {
    if (!is.na(indices$ox) && indices$ox == indices$x[1L]) return("consistency")
    if (is.na(indices$ox) && is.na(indices$oy)) return("margin")
    if (is.na(indices$ox)) return(if (indices$y[1L] == indices$oy) "4" else "5")
    if (is.na(indices$oy)) return("6")
    return("7")
  }
  if (!is.na(indices$ox) && indices$ox %in% indices$x) return("reduced by consistency")
  if (is.na(indices$ox) && is.na(indices$oy)) return("8")
  if (!is.na(indices$ox) && is.na(indices$oy)) return("9")
  if (is.na(indices$ox)) return("10")
  "11"
}
