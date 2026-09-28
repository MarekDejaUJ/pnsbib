# Finite longitudinal response-type oracle. Time is 1,...,T; Y_t follows A_t.
.poc_binary_paths <- function(periods, absorbing) {
  if (absorbing) {
    out <- t(vapply(0:periods, function(k)
      c(rep(0L, k), rep(1L, periods - k)), integer(periods)))
  } else {
    out <- as.matrix(do.call(expand.grid, rep(list(0:1), periods)))
    storage.mode(out) <- "integer"
  }
  rownames(out) <- apply(out, 1L, paste0, collapse = "")
  colnames(out) <- paste0("t", seq_len(periods))
  out
}

.poc_long_system <- function(model) {
  G <- nrow(model$regimes)
  H <- nrow(model$observed)
  P <- nrow(model$paths)
  ntypes <- P^G
  nvar <- H * ntypes
  if (!is.finite(nvar) || nvar > model$max_variables) {
    stop(sprintf("Longitudinal LP requires %.0f variables, above max_variables.",
                 nvar), call. = FALSE)
  }
  types <- as.matrix(do.call(expand.grid, rep(list(seq_len(P)), G)))
  storage.mode(types) <- "integer"
  if (model$no_anticipation && G > 1L) {
    keep <- rep(TRUE, nrow(types))
    for (g in seq_len(G - 1L)) for (h in (g + 1L):G) {
      for (t in seq_len(ncol(model$regimes))) {
        if (identical(unname(model$regimes[g, seq_len(t)]),
                      unname(model$regimes[h, seq_len(t)]))) {
          keep <- keep & model$paths[types[, g], t] ==
            model$paths[types[, h], t]
        }
      }
    }
    types <- types[keep, , drop = FALSE]
  }
  ntypes <- nrow(types)
  nvar <- H * ntypes
  if (nvar == 0L) stop("No response types satisfy temporal restrictions.", call. = FALSE)
  response <- types[rep(seq_len(ntypes), times = H), , drop = FALSE]
  factual <- rep(seq_len(H), each = ntypes)
  factual_regime <- match(rownames(model$observed), rownames(model$regimes))
  factual_path <- response[cbind(seq_len(nvar), factual_regime[factual])]
  A <- matrix(0, 1L, nvar)
  A[1L, ] <- 1
  b <- 1
  for (h in seq_len(H)) for (p in seq_len(P)) {
    A <- rbind(A, as.numeric(factual == h & factual_path == p))
    b <- c(b, model$observed[h, p])
  }
  if (nrow(model$interventional_paths)) for (i in seq_len(nrow(model$interventional_paths))) {
    g <- match(model$interventional_paths$regime[i], rownames(model$regimes))
    p <- match(model$interventional_paths$path[i], rownames(model$paths))
    A <- rbind(A, as.numeric(response[, g] == p))
    b <- c(b, model$interventional_paths$prob[i])
  }
  info <- .poc_long_time_data(model$time_margins)
  if (nrow(model$time_margins)) for (i in seq_len(nrow(model$time_margins))) {
    g <- match(model$time_margins$regime[i], rownames(model$regimes))
    t <- model$time_margins$horizon[i]
    A <- rbind(A, as.numeric(model$paths[response[, g], t] == info$outcome[i]))
    b <- c(b, info$prob[i])
  }
  list(A = A, b = b, response = response, factual = factual,
       factual_path = factual_path, nvar = nvar, ntypes = ntypes)
}

#' Construct a finite longitudinal probability-of-causation model
#'
#' Observed full outcome-path cells are joint probabilities P(H=h,Ybar=p).
#' Regime-specific intervention margins are optional and must be independently
#' justified. An absent margin is free, not estimated from observed risk.
#'
#' @param observed Data frame with history, path and prob columns; include all
#'   history/path cells, including zeros.
#' @param regimes Named matrix of categorical treatment states, one row per
#'   fixed history and one column per period. For two overlapping databases,
#'   use the disjoint joint-state labels `00`, `01`, `10`, `11`. Observed
#'   histories must be named regime rows.
#' @param interventional_paths Optional data frame regime, path, prob.
#' @param time_margins Optional data frame regime, horizon, prob1 for P(Y_t^g=1),
#'   or regime, horizon, outcome, prob for a named categorical state. Do not mix
#'   these schemas. Missing state/horizon margins remain unconstrained.
#' @param absorbing Restrict default binary paths to nondecreasing 0/1. Explicit
#'   categorical paths must be nondecreasing under outcome_order when TRUE.
#' @param allowed_paths Optional vector of admissible path row identifiers,
#'   for restrictions implied by the outcome definition itself.
#' @param no_anticipation Equal potential outcomes when treatment prefixes agree.
#' @param tolerance Probability tolerance.
#' @param max_variables Maximum LP variables.
#' @param outcome_paths Optional matrix of categorical states, with unique named
#'   rows identifying distinct complete paths and one column per model period.
#'   NULL retains automatic binary paths. This declares the potential path domain,
#'   not just the observed paths. Explicit matrices avoid ambiguous string parsing.
#' @param outcome_order Optional unique character state labels in increasing
#'   order. Required for explicit absorbing paths; never inferred from labels.
#' @inheritParams validate_poc_inputs
#' @return A pnsbib_longitudinal_model.
#' @export
poc_longitudinal_model <- function(observed, regimes,
                                   interventional_paths = NULL,
                                   time_margins = NULL,
                                   absorbing = FALSE,
                                   allowed_paths = NULL,
                                   no_anticipation = TRUE,
                                   tolerance = 1e-9,
                                   max_variables = 10000L,
                                   outcome_paths = NULL,
                                   outcome_order = NULL, lp_backend = "reference") {
  lp_backend <- .poc_lp_backend(lp_backend)
  if (!is.matrix(regimes) ||
      !(is.numeric(regimes) || is.character(regimes)) ||
      nrow(regimes) < 2L || ncol(regimes) < 1L ||
      is.null(rownames(regimes)) || anyNA(rownames(regimes)) ||
      any(!nzchar(rownames(regimes))) || anyDuplicated(rownames(regimes)) ||
      anyNA(regimes) ||
      (is.numeric(regimes) && any(!is.finite(regimes) |
                                  regimes != floor(regimes) | regimes < 0)) ||
      (is.character(regimes) && any(!nzchar(regimes))) ||
      anyDuplicated(apply(regimes, 1L, paste, collapse = "\r"))) {
    stop("regimes must be a named matrix of distinct categorical treatment histories.",
         call. = FALSE)
  }
  if (!is.logical(absorbing) || length(absorbing) != 1L || is.na(absorbing) ||
      !is.logical(no_anticipation) || length(no_anticipation) != 1L ||
      is.na(no_anticipation) || !is.numeric(tolerance) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      !is.numeric(max_variables) || length(max_variables) != 1L ||
      !is.finite(max_variables) || max_variables < 1) {
    stop("Invalid restrictions, tolerance or max_variables.", call. = FALSE)
  }
  domain <- .poc_long_paths(ncol(regimes), absorbing, outcome_paths,
                            outcome_order, max_variables)
  paths <- domain$paths
  if (!is.null(allowed_paths)) {
    if (!is.character(allowed_paths) || !length(allowed_paths) ||
        anyNA(allowed_paths) || anyDuplicated(allowed_paths) ||
        any(!allowed_paths %in% rownames(paths)))
      stop("allowed_paths must be unique admissible outcome paths.", call. = FALSE)
    paths <- paths[rownames(paths) %in% allowed_paths, , drop = FALSE]
  }
  if (!is.data.frame(observed) ||
      !all(c("history", "path", "prob") %in% names(observed)) ||
      nrow(observed) == 0L || anyNA(observed[, c("history", "path", "prob")]) ||
      !is.numeric(observed$prob) || any(!is.finite(observed$prob)) ||
      any(observed$prob < -tolerance | observed$prob > 1 + tolerance) ||
      any(!observed$history %in% rownames(regimes)) ||
      any(!observed$path %in% rownames(paths)) ||
      anyDuplicated(paste(observed$history, observed$path, sep = "\r"))) {
    stop("Invalid observed history/path probability table.", call. = FALSE)
  }
  histories <- rownames(regimes)[rownames(regimes) %in% observed$history]
  if (nrow(observed) != length(histories) * nrow(paths) ||
      abs(sum(observed$prob) - 1) > tolerance) {
    stop("Observed table must have every path cell for each history and sum to one.",
         call. = FALSE)
  }
  o <- matrix(NA_real_, length(histories), nrow(paths),
              dimnames = list(histories, rownames(paths)))
  for (i in seq_len(nrow(observed)))
    o[as.character(observed$history[i]), as.character(observed$path[i])] <-
      observed$prob[i]
  if (anyNA(o)) stop("Missing observed history/path cell.", call. = FALSE)
  if (is.null(interventional_paths))
    interventional_paths <- data.frame(regime = character(), path = character(),
                                      prob = numeric())
  if (is.null(time_margins))
    time_margins <- data.frame(regime = character(), horizon = integer(),
                               prob1 = numeric())
  a <- interventional_paths
  t <- time_margins
  if (!is.data.frame(a) || !all(c("regime", "path", "prob") %in% names(a)) ||
      anyNA(a) || !is.numeric(a$prob) || any(!is.finite(a$prob)) ||
      any(a$prob < -tolerance | a$prob > 1 + tolerance) ||
      any(!a$regime %in% rownames(regimes)) ||
      any(!a$path %in% rownames(paths)) ||
      anyDuplicated(paste(a$regime, a$path, sep = "\r"))) {
    stop("Invalid interventional path margins.", call. = FALSE)
  }
  time_info <- .poc_long_time_data(t)
  if (!is.data.frame(t) || !all(c("regime", "horizon") %in% names(t)) ||
      anyNA(t) || !is.numeric(t$horizon) || any(!is.finite(t$horizon)) ||
      any(t$horizon != floor(t$horizon) | t$horizon < 1L |
          t$horizon > ncol(regimes)) ||
      !is.numeric(time_info$prob) || any(!is.finite(time_info$prob)) ||
      any(time_info$prob < -tolerance | time_info$prob > 1 + tolerance) ||
      any(!time_info$outcome %in% domain$levels) ||
      any(!t$regime %in% rownames(regimes)) ||
      anyDuplicated(data.frame(regime = t$regime, horizon = t$horizon,
                               outcome = time_info$outcome))) {
    stop("Invalid interventional time margins.", call. = FALSE)
  }
  model <- structure(list(regimes = regimes, paths = paths, observed = o,
                          interventional_paths = a, time_margins = t,
                          absorbing = absorbing, no_anticipation = no_anticipation,
                          outcome_order = domain$order,
                          outcome_levels = domain$levels,
                          tolerance = tolerance, max_variables = max_variables,
                          lp_backend = lp_backend),
                     class = "pnsbib_longitudinal_model")
  system <- .poc_long_system(model)
  fit <- .poc_lp("min", rep(0, system$nvar), system$A,
                  rep("=", nrow(system$A)), system$b, lp_backend)
  if (fit$status != 0L)
    stop("No finite longitudinal response-type distribution fits supplied margins.",
         call. = FALSE)
  model$n_response_types <- system$ntypes
  model$n_variables <- system$nvar
  model
}

#' Specify a fixed-regime longitudinal probability-of-causation event
#'
#' @param kind One of pns, pn, ps, trajectory, persistent, event_time, lagged.
#' @param regime Active regime name.
#' @param reference Comparison regime name.
#' @param horizon Evaluation time 1,...,T.
#' @param start Start time for persistent event (default 1).
#' @param lag Lag for a lagged event (default 1).
#' @param regime_path,reference_path Exact path row identifiers for trajectory.
#' @param regime_outcome,reference_outcome Nonempty character sets of outcome
#'   labels for the active and reference events. Defaults preserve binary 1/0
#'   queries. PN conditions on the factual active set; PS on the factual reference
#'   set. Event-time means first entry into the active set, with reference-set
#'   membership at that horizon. It does not generally mean reference nonentry.
#'   Customized outcome sets are inapplicable to exact trajectory queries.
#' @param condition_on For PN/PS, condition on the named complete factual
#'   history (default "history") or all observed histories with its treatment
#'   prefix through the queried horizon ("prefix"). PN uses the active prefix
#'   and factual active-set membership; PS uses the reference prefix and factual
#'   reference-set membership (binary Y=1/0 by default). Only PN/PS
#'   accept prefix mode. With no anticipation, a prefix query is invariant to
#'   the named regimes' later continuations. Otherwise the counterfactual still
#'   refers to its named complete regime.
#' @return A pnsbib_longitudinal_query.
#' @examples
#' regimes <- rbind(off = c(0, 0), later = c(0, 1),
#'                  early = c(1, 0), always = c(1, 1))
#' observed <- expand.grid(history = rownames(regimes), path = c("00", "11"),
#'                         stringsAsFactors = FALSE)
#' observed$prob <- c(.1, .2, 0, 0, 0, 0, .3, .4)
#' margins <- data.frame(regime = rownames(regimes), horizon = 1,
#'                       prob1 = c(.3, .3, .8, .8))
#' model <- poc_longitudinal_model(observed, regimes, time_margins = margins,
#'                                 allowed_paths = c("00", "11"))
#' # Pool early and always for the factual first-period treatment condition.
#' query <- poc_longitudinal_query("pn", "early", "off", 1,
#'                                 condition_on = "prefix")
#' poc_longitudinal_bounds(model, query)
#' @export
poc_longitudinal_query <- function(kind = c("pns", "pn", "ps", "trajectory",
                                            "persistent", "event_time", "lagged"),
                                   regime, reference, horizon,
                                   start = 1L, lag = 1L,
                                   regime_path = NULL, reference_path = NULL,
                                   condition_on = c("history", "prefix"),
                                   regime_outcome = "1", reference_outcome = "0") {
  kind <- match.arg(kind)
  condition_on <- match.arg(condition_on)
  if (condition_on == "prefix" && !kind %in% c("pn", "ps"))
    stop("Prefix conditioning requires PN or PS.", call. = FALSE)
  if (!.poc_long_labels(regime_outcome) || !.poc_long_labels(reference_outcome))
    stop("Outcome sets must contain unique nonempty character labels.", call. = FALSE)
  if (kind == "trajectory" && (!identical(regime_outcome, "1") ||
                                !identical(reference_outcome, "0")))
    stop("Outcome sets are inapplicable to exact trajectory queries.", call. = FALSE)
  if (!is.character(regime) || length(regime) != 1L || is.na(regime) ||
      !nzchar(regime) || !is.character(reference) || length(reference) != 1L ||
      is.na(reference) || !nzchar(reference) || regime == reference ||
      !is.numeric(horizon) || length(horizon) != 1L || !is.finite(horizon) ||
      horizon < 1 || horizon > .Machine$integer.max || horizon != floor(horizon) ||
      !is.numeric(start) || length(start) != 1L || !is.finite(start) ||
      start < 1 || start > horizon || start != floor(start) ||
      !is.numeric(lag) || length(lag) != 1L || !is.finite(lag) ||
      lag < 1 || lag > .Machine$integer.max || lag != floor(lag)) {
    stop("Invalid longitudinal query labels or time indices.", call. = FALSE)
  }
  if (kind == "lagged" && lag >= horizon)
    stop("lag must be smaller than horizon.", call. = FALSE)
  if (kind == "trajectory" &&
      (!is.character(regime_path) || length(regime_path) != 1L ||
       !is.character(reference_path) || length(reference_path) != 1L ||
       anyNA(c(regime_path, reference_path)))) {
    stop("A trajectory query needs two exact path labels.", call. = FALSE)
  }
  structure(list(kind = kind, regime = regime, reference = reference,
                 horizon = as.integer(horizon), start = as.integer(start),
                 lag = as.integer(lag), regime_path = regime_path,
                 reference_path = reference_path, condition_on = condition_on,
                 regime_outcome = regime_outcome, reference_outcome = reference_outcome),
            class = "pnsbib_longitudinal_query")
}

.poc_long_event <- function(model, query, system = NULL) {
  if (!inherits(model, "pnsbib_longitudinal_model") ||
      !inherits(query, "pnsbib_longitudinal_query"))
    stop("Expected a longitudinal model and query.", call. = FALSE)
  g <- match(query$regime, rownames(model$regimes))
  r <- match(query$reference, rownames(model$regimes))
  T <- ncol(model$regimes)
  t <- query$horizon
  if (anyNA(c(g, r)) || t > T ||
      (query$kind == "trajectory" &&
       (is.na(match(query$regime_path, rownames(model$paths))) ||
        is.na(match(query$reference_path, rownames(model$paths))))))
    stop("Query regime, horizon or path absent from model.", call. = FALSE)
  if (is.null(system)) system <- .poc_long_system(model)
  path_g <- system$response[, g]
  path_r <- system$response[, r]
  Yg <- model$paths[path_g, , drop = FALSE]
  Yr <- model$paths[path_r, , drop = FALSE]
  kind <- query$kind
  # Saved <=0.0.5 queries retain their binary meaning.
  active <- if (is.null(query$regime_outcome)) "1" else query$regime_outcome
  reference <- if (is.null(query$reference_outcome)) "0" else query$reference_outcome
  levels <- if (is.null(model$outcome_levels)) c("0", "1") else model$outcome_levels
  if (kind != "trajectory" &&
      (any(!c(active, reference) %in% levels) ||
       !.poc_long_labels(active) || !.poc_long_labels(reference)))
    stop("Query outcome set absent from model path domain.", call. = FALSE)
  in_active <- matrix(Yg %in% active, nrow(Yg), ncol(Yg))
  in_reference <- matrix(Yr %in% reference, nrow(Yr), ncol(Yr))
  event <- switch(kind,
    pns = in_active[, t] & in_reference[, t],
    pn = in_reference[, t],
    ps = in_active[, t],
    trajectory = path_g == match(query$regime_path, rownames(model$paths)) &
      path_r == match(query$reference_path, rownames(model$paths)),
    persistent = rowSums(in_active[, query$start:t, drop = FALSE] &
                         in_reference[, query$start:t, drop = FALSE]) ==
      t - query$start + 1L,
    event_time = in_active[, t] &
      (if (t == 1L) TRUE else rowSums(in_active[, seq_len(t - 1L), drop = FALSE]) == 0L) &
      in_reference[, t],
    lagged = in_active[, t] & in_reference[, t - query$lag])
  denom <- 1
  condition <- rep(TRUE, system$nvar)
  condition_on <- "none"
  condition_horizon <- NA_integer_
  condition_label <- ""
  if (kind %in% c("pn", "ps")) {
    h <- if (kind == "pn") query$regime else query$reference
    wanted <- if (kind == "pn") active else reference
    # Serialized queries from <= 0.0.3 retain their full-history meaning.
    condition_on <- if (is.null(query$condition_on)) "history" else
      query$condition_on
    histories <- rownames(model$observed)
    if (condition_on == "prefix") {
      prefix <- model$regimes[h, seq_len(t)]
      matching <- vapply(histories, function(name)
        all(model$regimes[name, seq_len(t)] == prefix), logical(1L))
      condition_horizon <- as.integer(t)
      condition_label <- sprintf("A[1:%d]=(%s)", t,
                                 paste(prefix, collapse = ","))
    } else {
      matching <- histories == h
      condition_horizon <- as.integer(T)
      condition_label <- paste0("H=", h)
    }
    if (!any(matching))
      stop("Conditional query requires observed factual history matching its condition.",
           call. = FALSE)
    denom <- sum(model$observed[matching, model$paths[, t] %in% wanted, drop = FALSE])
    condition <- system$factual %in% which(matching) &
      model$paths[system$factual_path, t] %in% wanted
    event <- event & condition
  }
  ag <- .poc_long_set_label(active)
  ar <- .poc_long_set_label(reference)
  label <- switch(kind,
    pns = sprintf("Y_%s(%d)%s & Y_%s(%d)%s", query$regime, t, ag, query$reference, t, ar),
    pn = sprintf("Y_%s(%d)%s | %s,Y(%d)%s", query$reference, t, ar,
                 condition_label, t, ag),
    ps = sprintf("Y_%s(%d)%s | %s,Y(%d)%s", query$regime, t, ag,
                 condition_label, t, ar),
    trajectory = sprintf("Y_%s=%s & Y_%s=%s", query$regime,
                         query$regime_path, query$reference,
                         query$reference_path),
    persistent = sprintf("Y_%s%s & Y_%s%s for t=%d..%d", query$regime, ag,
                         query$reference, ar, query$start, t),
    event_time = sprintf("first Y_%s%s at %d; Y_%s(%d)%s", query$regime, ag,
                         t, query$reference, t, ar),
    lagged = sprintf("Y_%s(%d)%s & Y_%s(%d)%s", query$regime, t, ag,
                     query$reference, t - query$lag, ar))
  list(event = as.numeric(event), condition = as.numeric(condition),
       conditional = kind %in% c("pn", "ps"), denominator = denom,
       condition_on = condition_on, condition_horizon = condition_horizon,
       regime_outcome = if (kind == "trajectory") NA_character_ else paste(active, collapse = ";"),
       reference_outcome = if (kind == "trajectory") NA_character_ else paste(reference, collapse = ";"),
       outcome_order = if (is.null(model$outcome_order)) NA_character_ else
         paste(model$outcome_order, collapse = ";"),
       label = label, system = system)
}

#' Exact longitudinal response-type LP bounds
#'
#' Sharp conditional on supplied margins and restrictions; no sampling
#' uncertainty or causal identification is implied. PN and PS condition on
#' a complete factual history or its prefix, as specified by the query.
#'
#' @param model A pnsbib_longitudinal_model.
#' @param query A pnsbib_longitudinal_query.
#' @return A one-row pnsbib_longitudinal_result data frame.
#' @inheritParams poc_exact
#' @export
poc_longitudinal_bounds <- function(model, query, lp_backend = NULL) {
  lp_backend <- .poc_lp_backend(lp_backend, model)
  specification <- .poc_long_event(model, query)
  system <- specification$system
  event <- specification$event
  denom <- specification$denominator
  label <- specification$label
  kind <- query$kind
  t <- query$horizon
  if (kind %in% c("pn", "ps") && denom <= model$tolerance) {
    return(structure(data.frame(estimand = toupper(kind), event = label,
      regime = query$regime, comparison = query$reference, horizon = t,
      lower = NA_real_, upper = NA_real_, denominator = denom,
      condition_on = specification$condition_on,
      condition_horizon = specification$condition_horizon,
      regime_outcome = specification$regime_outcome,
      reference_outcome = specification$reference_outcome,
      outcome_order = specification$outcome_order,
      method = "longitudinal_exact_lp", sharpness = "sharp_given_inputs",
      status = "undefined_condition", assumptions = paste0(
        "consistency;absorbing=", model$absorbing,
        ";no_anticipation=", model$no_anticipation),
      stringsAsFactors = FALSE), class = c("pnsbib_longitudinal_result", "data.frame")))
  }
  lo <- .poc_lp("min", as.numeric(event), system$A,
                 rep("=", nrow(system$A)), system$b, lp_backend)
  hi <- .poc_lp("max", as.numeric(event), system$A,
                 rep("=", nrow(system$A)), system$b, lp_backend)
  status <- if (lo$status == 0L && hi$status == 0L) "ok" else "lp_failed"
  structure(data.frame(estimand = toupper(kind), event = label,
    regime = query$regime, comparison = query$reference, horizon = t,
    lower = if (status == "ok") max(0, lo$objval / denom) else NA_real_,
    upper = if (status == "ok") min(1, hi$objval / denom) else NA_real_,
    denominator = denom, method = "longitudinal_exact_lp",
    condition_on = specification$condition_on,
    condition_horizon = specification$condition_horizon,
    regime_outcome = specification$regime_outcome,
    reference_outcome = specification$reference_outcome,
    outcome_order = specification$outcome_order,
    sharpness = "sharp_given_inputs", status = status,
    assumptions = paste0("consistency;absorbing=", model$absorbing,
                         ";no_anticipation=", model$no_anticipation),
    stringsAsFactors = FALSE),
    class = c("pnsbib_longitudinal_result", "data.frame"))
}
