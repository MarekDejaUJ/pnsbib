#' Map structural template margins to confidence-region cell identifiers
#'
#' The returned order matches the finite response-type system. Identifiers are
#' local to this template: retain the full map when preparing sample vectors.
#' All numerical template probabilities are ignored by the confidence API.
#' Template support and temporal restrictions remain substantive assumptions.
#'
#' @param template A static [poc_model()] or [poc_longitudinal_model()].
#'   Partial-model objects are not supported by this interface.
#' @return A data frame with cell_id, type and explicit event labels.
#' @export
poc_margin_map <- function(template) {
  rows <- list()
  add <- function(type, treatment = NA_character_, outcome = NA_character_,
                  history = NA_character_, path = NA_character_,
                  regime = NA_character_, horizon = NA_integer_) {
    rows[[length(rows) + 1L]] <<- data.frame(type, treatment, outcome, history,
      path, regime, horizon, stringsAsFactors = FALSE)
  }
  if (inherits(template, "pnsbib_model")) {
    for (type in c("observed", "intervention"))
      for (x in rownames(template$o)) for (y in colnames(template$o))
        add(type, treatment = x, outcome = y)
  } else if (inherits(template, "pnsbib_longitudinal_model")) {
    for (h in rownames(template$observed)) for (p in rownames(template$paths))
      add("observed_path", history = h, path = p)
    a <- template$interventional_paths
    if (nrow(a)) for (i in seq_len(nrow(a)))
      add("intervention_path", regime = as.character(a$regime[i]), path = as.character(a$path[i]))
    t <- template$time_margins
    info <- .poc_long_time_data(t)
    if (nrow(t)) for (i in seq_len(nrow(t)))
      add("intervention_time", regime = as.character(t$regime[i]),
          horizon = t$horizon[i], outcome = info$outcome[i])
  } else stop("Expected a static or longitudinal structural template.", call. = FALSE)
  out <- do.call(rbind, rows)
  out$cell_id <- paste0("cell_", seq_len(nrow(out)))
  out[, c("cell_id", setdiff(names(out), "cell_id"))]
}

#' Construct a simultaneous sampling confidence region for model margins
#'
#' The template supplies labels and response support only. Its probability
#' values are ignored, including margins without samples. Missing sample cells
#' and numeric(0) leave probabilities free in `[0,1]`. Thus a region can be used
#' even when the empirical cell centers are mutually incompatible.
#'
#' Hoeffding uses independent `[0,1]` unit contributions and fixed nonnegative
#' weights (equal weights by default). Its radius is
#' sqrt(sum(w^2) * log(2 * J / alpha) / 2), where J is the complete margin-map
#' row count and alpha = 1 - conf_level. Clopper-Pearson uses unweighted iid
#' Bernoulli contributions and per-cell error alpha/J. A union bound supplies
#' at least conf_level simultaneous coverage, without independence across cells.
#'
#' Cluster vectors contain one bounded contribution per independent cluster,
#' not one value per dependent paper. Equal cluster weighting generally targets
#' a different population from equal paper weighting. Weights must be fixed
#' externally; estimated propensity weights and outcome-dependent selection or
#' stopping are not covered. Sampling and identification assumptions require
#' external justification. This does not supply g-formula standard errors.
#'
#' @param template A static or longitudinal structural template.
#' @param samples Named list of numeric `[0,1]` contribution vectors, keyed by
#'   [poc_margin_map()] cell_id. No missing values or automatic deletion.
#' @param target Nonempty description of the common target population.
#' @param sampling_unit Explicitly "individual" or "cluster".
#' @param method "hoeffding" or "clopper_pearson". The latter accepts only
#'   unweighted binary individual contributions.
#' @param conf_level Simultaneous confidence level strictly between zero and one.
#' @param weights Optional named list of fixed nonnegative weight vectors for
#'   supplied cells; only for Hoeffding. Omitted weights are equal.
#' @param assume_sampling Must be TRUE to declare the documented sampling and
#'   common-target assumptions. This records a declaration, not a verification.
#' @return A pnsbib_confidence_region with a rectangular cells table, sampling
#'   metadata and structural template. Raw sample vectors are not retained.
#' @references Hoeffding (1963), doi:10.1080/01621459.1963.10500830;
#'   [stats::binom.test()] for Clopper-Pearson intervals.
#' @examples
#' o <- matrix(.25, 2, 2, dimnames = list(c("0", "1"), c("0", "1")))
#' template <- poc_model(o, o * 2)
#' map <- poc_margin_map(template)
#' samples <- setNames(lapply(1:4, function(j) as.numeric(rep(1:4, 25) == j)),
#'                     map$cell_id[map$type == "observed"])
#' region <- poc_confidence_region(template, samples,
#'   target = "Illustrative iid article population", sampling_unit = "individual",
#'   method = "clopper_pearson", assume_sampling = TRUE)
#' region$cells
#' poc_confidence_bounds(region, poc_query(c("0" = "0", "1" = "1")))
#' @export
poc_confidence_region <- function(template, samples, target, sampling_unit,
                                  method = c("hoeffding", "clopper_pearson"),
                                  conf_level = .95, weights = NULL,
                                  assume_sampling = FALSE) {
  method <- match.arg(method)
  if (missing(sampling_unit) || !is.character(sampling_unit) ||
      length(sampling_unit) != 1L || is.na(sampling_unit) ||
      !sampling_unit %in% c("individual", "cluster"))
    stop("Declare sampling_unit as individual or cluster.", call. = FALSE)
  if (!identical(assume_sampling, TRUE))
    stop("Declare the sampling assumptions with assume_sampling = TRUE.", call. = FALSE)
  if (missing(target) || !is.character(target) || length(target) != 1L ||
      is.na(target) || !nzchar(trimws(target)) || !is.numeric(conf_level) ||
      length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("Supply a target description and confidence level strictly between zero and one.", call. = FALSE)
  cells <- poc_margin_map(template)
  valid_list <- function(x) is.list(x) && !is.data.frame(x) &&
    (length(x) == 0L || .poc_long_labels(names(x)))
  if (!valid_list(samples) || any(!names(samples) %in% cells$cell_id))
    stop("samples must be a named list keyed by template cell identifiers.", call. = FALSE)
  if (!is.null(weights) && (!valid_list(weights) || any(!names(weights) %in% names(samples))))
    stop("weights must be a named list for supplied sample cells.", call. = FALSE)
  if (method == "clopper_pearson" && (!is.null(weights) || sampling_unit != "individual"))
    stop("Clopper-Pearson requires unweighted binary individual contributions.", call. = FALSE)
  cells$n_units <- 0L
  cells$effective_units <- 0
  cells$estimate <- NA_real_
  cells$lower <- 0
  cells$upper <- 1
  J <- nrow(cells)
  alpha <- 1 - conf_level
  for (name in names(samples)) {
    z <- samples[[name]]
    if (!is.numeric(z) || !is.null(dim(z)) || any(!is.finite(z)) || any(z < 0 | z > 1))
      stop("Each sample must be a complete numeric vector in [0,1].", call. = FALSE)
    i <- match(name, cells$cell_id)
    n <- length(z)
    w <- weights[[name]]
    if (!is.null(w) && (!is.numeric(w) || !is.null(dim(w)) || length(w) != n ||
        !length(w) || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0 || !is.finite(sum(w))))
      stop("Weights must match a nonempty sample and have finite positive sum.", call. = FALSE)
    if (!n) next
    if (is.null(w)) w <- rep(1 / n, n) else w <- w / sum(w)
    cells$n_units[i] <- n
    cells$effective_units[i] <- 1 / sum(w^2)
    cells$estimate[i] <- sum(w * z)
    if (method == "hoeffding") {
      radius <- sqrt(sum(w^2) * (log(2 * J) - log(alpha)) / 2)
      bounds <- c(max(0, cells$estimate[i] - radius), min(1, cells$estimate[i] + radius))
    } else {
      if (any(!z %in% c(0, 1)))
        stop("Clopper-Pearson requires binary contributions.", call. = FALSE)
      successes <- sum(z)
      # Tail probabilities avoid rounding a near-one confidence level to one.
      tail <- alpha / (2 * J)
      bounds <- c(if (successes == 0) 0 else stats::qbeta(tail, successes, n - successes + 1),
                  if (successes == n) 1 else stats::qbeta(tail, successes + 1, n - successes,
                                                        lower.tail = FALSE))
    }
    cells$lower[i] <- bounds[1L]
    cells$upper[i] <- bounds[2L]
  }
  structure(list(template = template, cells = cells, target = target,
    sampling_unit = sampling_unit, method = method, conf_level = conf_level,
    family_cells = J, assume_sampling = TRUE,
    assumptions = paste0("independent_", sampling_unit, "_contributions;fixed_sampling_design;",
      if (method == "hoeffding") "bounded_contributions;fixed_weights;" else "iid_bernoulli;equal_weights;",
      "same_target_margins;valid_causal_model;template_probabilities_ignored")),
    class = "pnsbib_confidence_region")
}
