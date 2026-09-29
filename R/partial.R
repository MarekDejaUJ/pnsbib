#' Specify a partial probability constraint
#'
#' A constraint is a declared input, not an identification estimator or a
#' confidence interval. Missing probabilities remain unknown. For example,
#' P(Y_0=1 | X=1) constrains a treated-subgroup counterfactual risk, not the
#' whole-population P(Y_0=1).
#'
#' @param query A [poc_query()] event, optionally conditional on factual X/Y.
#' @param lower,upper Probability endpoints in `[0,1]`; upper defaults to lower
#'   for an equality.
#' @param stratum Optional stratum label. The probability is conditional on
#'   that stratum, and additionally on factual X/Y if query is conditional.
#'   NULL refers to the whole target population.
#' @param label Optional unique, nonempty constraint label for diagnostics.
#' @return A pnsbib_constraint for [poc_partial_model()].
#' @export
poc_constraint <- function(query, lower, upper = lower, stratum = NULL,
                           label = NULL) {
  probability <- function(value) is.numeric(value) && length(value) == 1L &&
    is.finite(value) && value >= 0 && value <= 1
  text <- function(value) is.null(value) || (is.character(value) &&
    length(value) == 1L && !is.na(value) && nzchar(value))
  if (!inherits(query, "pnsbib_query") || !probability(lower) ||
      !probability(upper) || lower > upper || !text(stratum) || !text(label))
    stop("Expected a query, ordered probability endpoints and optional labels.",
         call. = FALSE)
  structure(list(query = query, lower = lower, upper = upper,
                 stratum = stratum, label = label), class = "pnsbib_constraint")
}

#' Construct a finite model with partial or subgroup probability constraints
#'
#' Supply fixed observational probabilities and only the additional probabilities
#' actually justified for the declared target population. Unspecified intervention
#' margins remain free. All constraints are tested jointly for feasibility.
#' This class is separate from the complete-margin [poc_model()] used by the
#' Li-Pearl analytical formulas. Use [poc_partial_bounds()] for its queries.
#'
#' Strata are a disjoint, exhaustive, prespecified population partition. Their
#' observed cells are joint P(Z=s,X=x,Y=y), not separately normalized conditional
#' tables. The package does not establish causal identification, estimate strata
#' weights, impose exchangeability or supply sampling uncertainty.
#'
#' @param observed_joint A named numeric X-by-Y joint probability matrix, or
#'   a named list of conformable matrices whose cells sum to one across strata.
#'   Include every factual cell, including zeros; missing cells are invalid.
#' @param constraints A list of [poc_constraint()] objects; an empty list leaves
#'   all additional potential-outcome probabilities unconstrained.
#' @param tolerance Positive numerical tolerance for normalization and conditions.
#' @param max_variables Limit on stratum-by-factual-X-by-response-type variables.
#' @inheritParams validate_poc_inputs
#' @return A pnsbib_partial_model, with the supplied constraints and a
#'   constraint_table showing their probabilities, denominators and labels.
#' @examples
#' observed <- matrix(c(.3, .2, .1, .4), 2, 2, byrow = TRUE,
#'                    dimnames = list(c("0", "1"), c("0", "1")))
#' treated_risk <- poc_constraint(
#'   poc_query(c("0" = "1"), observed_x = "1", conditional = TRUE),
#'   lower = .4, label = "hypothetical_treated_counterfactual")
#' model <- poc_partial_model(observed, list(treated_risk))
#' poc_partial_bounds(model, poc_query(c("0" = "0", "1" = "1")))
#' @export
poc_partial_model <- function(observed_joint, constraints = list(),
                               tolerance = 1e-9, max_variables = 10000L,
                               lp_backend = "reference") {
  lp_backend <- .poc_lp_backend(lp_backend)
  valid_names <- function(x) !is.null(x) && !anyNA(x) &&
    all(nzchar(x)) && !anyDuplicated(x)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0 ||
      !is.numeric(max_variables) || length(max_variables) != 1L ||
      !is.finite(max_variables) || max_variables < 1)
    stop("Invalid tolerance or max_variables.", call. = FALSE)
  strata <- if (is.matrix(observed_joint)) list(population = observed_joint) else
    observed_joint
  if (!is.list(strata) || is.data.frame(strata) || !length(strata) ||
      !valid_names(names(strata)))
    stop("observed_joint needs a matrix or uniquely named list of stratum matrices.",
         call. = FALSE)
  good_matrix <- function(x) is.matrix(x) && is.numeric(x) &&
    nrow(x) >= 2L && ncol(x) >= 2L && valid_names(rownames(x)) &&
    valid_names(colnames(x)) && all(is.finite(x)) && all(x >= 0 & x <= 1)
  if (!all(vapply(strata, good_matrix, logical(1L))) ||
      !all(vapply(strata, function(x) identical(dimnames(x), dimnames(strata[[1L]])),
                  logical(1L))))
    stop("Strata require conformable named probability matrices with matching labels.",
         call. = FALSE)
  if (abs(sum(vapply(strata, sum, 0)) - 1) > tolerance)
    stop("Joint probabilities must sum to one across all strata.", call. = FALSE)
  if (!is.list(constraints) || !all(vapply(constraints, inherits, logical(1L),
                                           "pnsbib_constraint")))
    stop("constraints must be a list of poc_constraint objects.", call. = FALSE)
  for (i in seq_along(constraints)) {
    z <- constraints[[i]]
    # Revalidate a deserialized or manually modified constraint.
    constraints[[i]] <- poc_constraint(z$query, z$lower, z$upper, z$stratum,
      if (is.null(z$label)) paste0("constraint_", i) else z$label)
  }
  labels <- vapply(constraints, function(z) z$label, "")
  if (anyDuplicated(labels)) stop("Constraint labels must be unique.", call. = FALSE)
  model <- structure(list(strata = strata, o = Reduce(`+`, strata),
                           constraints = constraints, tolerance = tolerance,
                           max_variables = max_variables, lp_backend = lp_backend),
                      class = "pnsbib_partial_model")
  system <- .poc_partial_system(model)
  fit <- .poc_lp("min", rep(0, system$nvar), system$A, system$direction, system$b, lp_backend)
  if (fit$status != 0L)
    stop(if (fit$status == 2L) "No response-type distribution satisfies the partial constraints." else
      "LP solver failed while checking partial constraints.", call. = FALSE)
  model$n_variables <- system$nvar
  model$n_response_types <- system$ntypes
  model$constraint_table <- system$constraint_table
  model
}

.poc_partial_event <- function(model, query, stratum, system) {
  idx <- .poc_query_indices(rownames(model$o), colnames(model$o), query)
  if (!is.null(stratum) && (!is.character(stratum) || length(stratum) != 1L ||
      is.na(stratum) || !stratum %in% names(model$strata)))
    stop("Query stratum absent from model.", call. = FALSE)
  selected <- if (is.null(stratum)) rep(TRUE, system$nvar) else
    system$stratum == match(stratum, names(model$strata))
  o <- if (is.null(stratum)) model$o else model$strata[[stratum]]
  denom <- if (!query$conditional) sum(o) else {
    x <- if (is.na(idx$ox)) seq_len(nrow(o)) else idx$ox
    y <- if (is.na(idx$oy)) seq_len(ncol(o)) else idx$oy
    sum(o[x, y, drop = FALSE])
  }
  list(event = .poc_event_indicator(system, idx) * selected,
       denominator = denom)
}

.poc_partial_system <- function(model) {
  m <- nrow(model$o)
  n <- ncol(model$o)
  S <- length(model$strata)
  ntypes <- n^m
  nvar <- S * m * ntypes
  if (!is.finite(nvar) || nvar > model$max_variables)
    stop(sprintf("Partial LP requires %.0f variables, above max_variables.", nvar),
         call. = FALSE)
  types <- as.matrix(do.call(expand.grid, rep(list(seq_len(n)), m)))
  storage.mode(types) <- "integer"
  response <- types[rep(seq_len(ntypes), times = S * m), , drop = FALSE]
  factual_x <- rep(rep(seq_len(m), each = ntypes), times = S)
  factual_y <- response[cbind(seq_len(nvar), factual_x)]
  stratum <- rep(seq_len(S), each = m * ntypes)
  system <- list(response = response, factual_x = factual_x, factual_y = factual_y,
                  stratum = stratum, nvar = nvar, ntypes = ntypes)
  A <- matrix(0, 1L + S * m * n, nvar)
  A[1L, ] <- 1
  b <- c(1, numeric(S * m * n))
  direction <- rep("=", nrow(A))
  row <- 1L
  for (s in seq_len(S)) for (x in seq_len(m)) for (y in seq_len(n)) {
    row <- row + 1L
    A[row, ] <- as.numeric(stratum == s & factual_x == x & factual_y == y)
    b[row] <- model$strata[[s]][x, y]
  }
  table <- data.frame(label = character(), stratum = character(),
    lower = numeric(), upper = numeric(), denominator = numeric(),
    event = character(), observed_x = character(), observed_y = character(),
    conditional = logical(), stringsAsFactors = FALSE)
  for (constraint in model$constraints) {
    query <- constraint$query
    e <- .poc_partial_event(model, query, constraint$stratum, system)
    if (e$denominator <= model$tolerance)
      stop(sprintf("Constraint '%s' has a zero-probability condition.", constraint$label),
           call. = FALSE)
    if (constraint$lower == constraint$upper) {
      A <- rbind(A, e$event)
      b <- c(b, constraint$lower * e$denominator)
      direction <- c(direction, "=")
    } else {
      A <- rbind(A, e$event, e$event)
      b <- c(b, c(constraint$lower, constraint$upper) * e$denominator)
      direction <- c(direction, ">=", "<=")
    }
    table <- rbind(table, data.frame(label = constraint$label,
      stratum = if (is.null(constraint$stratum)) NA_character_ else constraint$stratum,
      lower = constraint$lower, upper = constraint$upper, denominator = e$denominator,
      event = paste(paste0("Y_", names(query$counterfactual), "=",
                           unname(query$counterfactual)), collapse = " & "),
      observed_x = if (is.null(query$observed_x)) NA_character_ else query$observed_x,
      observed_y = if (is.null(query$observed_y)) NA_character_ else query$observed_y,
      conditional = query$conditional, stringsAsFactors = FALSE))
  }
  c(system, list(A = A, b = b, direction = direction, constraint_table = table))
}

#' Exact bounds from partial and subgroup probability information
#'
#' Optimizes the complete finite response-type distribution subject to fixed
#' observed cells and the model's declared constraints. Missing margins remain
#' free. The result is sharp for those inputs, not evidence of their causal
#' validity or a statistical confidence interval. This function does not use
#' the complete-margin Li-Pearl analytical formulas.
#'
#' @param model A [poc_partial_model()].
#' @param query A [poc_query()].
#' @param stratum Optional stratum to condition on; NULL selects the whole
#'   population. A conditional query also conditions on its factual X/Y event.
#' @param return_witness Include lower/upper endpoint distributions in the
#'   witnesses attribute. Each row gives stratum, factual treatment, response
#'   type, potential outcomes and joint probability.
#' @inheritParams poc_exact
#' @return A one-row pnsbib_result with conditioning stratum and input labels.
#' @export
poc_partial_bounds <- function(model, query, stratum = NULL, return_witness = FALSE,
                                lp_backend = NULL) {
  lp_backend <- .poc_lp_backend(lp_backend, model)
  if (!inherits(model, "pnsbib_partial_model") ||
      !is.logical(return_witness) || length(return_witness) != 1L ||
      is.na(return_witness))
    stop("Expected a partial model and a logical return_witness flag.", call. = FALSE)
  system <- .poc_partial_system(model)
  event <- .poc_partial_event(model, query, stratum, system)
  denominator <- event$denominator
  result <- function(lower, upper, status) {
    out <- .poc_result(lower, upper, "partial_exact_lp", "sharp_given_constraints",
                       status, denominator, query)
    out$stratum <- if (is.null(stratum)) NA_character_ else stratum
    out$input_information <- "partial_probability_constraints"
    out$n_constraints <- length(model$constraints)
    out$constraint_labels <- paste(model$constraint_table$label, collapse = ";")
    out$assumptions <- "consistency;fixed_observed_joint;declared_probability_constraints"
    out
  }
  if (denominator <= model$tolerance)
    return(result(NA_real_, NA_real_, "undefined_condition"))
  lo <- .poc_lp("min", event$event, system$A, system$direction, system$b, lp_backend)
  hi <- .poc_lp("max", event$event, system$A, system$direction, system$b, lp_backend)
  if (lo$status != 0L || hi$status != 0L)
    return(result(NA_real_, NA_real_, if (lo$status == 2L || hi$status == 2L)
      "incompatible_constraints" else "solver_failure"))
  out <- result(max(0, lo$objval / denominator), min(1, hi$objval / denominator), "ok")
  if (return_witness) {
    base <- data.frame(stratum = names(model$strata)[system$stratum],
      factual_x = rownames(model$o)[system$factual_x],
      response_type = rep(seq_len(system$ntypes), times = length(model$strata) * nrow(model$o)),
      stringsAsFactors = FALSE)
    for (x in seq_len(nrow(model$o)))
      base[[paste0("Y_", rownames(model$o)[x])]] <- colnames(model$o)[system$response[, x]]
    witness <- function(fit) {
      table <- base
      table$prob <- fit$solution
      table
    }
    attr(out, "witnesses") <- list(lower = witness(lo), upper = witness(hi))
  }
  out
}
