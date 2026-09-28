# Explicit categorical path domains and marginal schemas. No coarsening of
# categorical states is used in consistency or shared-prefix restrictions.
.poc_long_labels <- function(x) is.character(x) && length(x) > 0L &&
  !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)

.poc_long_paths <- function(periods, absorbing, outcome_paths, outcome_order,
                            max_variables) {
  if (is.null(outcome_paths)) {
    if (!is.null(outcome_order) && !identical(outcome_order, c("0", "1")))
      stop("Default binary paths use outcome_order c('0', '1').", call. = FALSE)
    size <- if (absorbing) periods + 1 else 2^periods
    if (!is.finite(size) || size > max_variables)
      stop("Binary path enumeration exceeds max_variables; supply explicit outcome_paths.", call. = FALSE)
    return(list(paths = .poc_binary_paths(periods, absorbing), levels = c("0", "1"),
                order = if (absorbing) c("0", "1") else NULL))
  }
  p <- outcome_paths
  if (!is.matrix(p) || !(is.character(p) || is.numeric(p)) ||
      nrow(p) < 1L || ncol(p) != periods ||
      !.poc_long_labels(rownames(p)) || anyNA(p) ||
      (is.numeric(p) && any(!is.finite(p))) ||
      (is.character(p) && any(!nzchar(p))))
    stop("outcome_paths must be a row-named, complete categorical path matrix.", call. = FALSE)
  p <- matrix(as.character(p), nrow(p), ncol(p), dimnames = dimnames(p))
  if (anyDuplicated(as.data.frame(p)))
    stop("outcome_paths must represent distinct state sequences.", call. = FALSE)
  if (!is.null(outcome_order) &&
      (!.poc_long_labels(outcome_order) || any(!p %in% outcome_order)))
    stop("outcome_order must list every path state once, in the intended order.", call. = FALSE)
  if (absorbing) {
    if (is.null(outcome_order))
      stop("Explicit absorbing paths require an outcome_order.", call. = FALSE)
    ranks <- matrix(match(p, outcome_order), nrow(p), ncol(p))
    if (ncol(p) > 1L && any(ranks[, -1L, drop = FALSE] <
                            ranks[, -ncol(p), drop = FALSE]))
      stop("outcome_paths violate the absorbing outcome_order.", call. = FALSE)
  }
  list(paths = p, order = outcome_order,
       levels = if (is.null(outcome_order)) unique(as.character(p)) else outcome_order)
}

.poc_long_time_data <- function(t) {
  if (!is.data.frame(t) || !all(c("regime", "horizon") %in% names(t)))
    stop("Invalid interventional time margins.", call. = FALSE)
  legacy <- "prob1" %in% names(t)
  categorical <- any(c("outcome", "prob") %in% names(t))
  if (legacy == categorical ||
      (categorical && !all(c("outcome", "prob") %in% names(t))))
    stop("Time margins need either prob1 or outcome/prob, not mixed schemas.", call. = FALSE)
  list(outcome = if (legacy) rep("1", nrow(t)) else as.character(t$outcome),
       prob = if (legacy) t$prob1 else t$prob)
}

.poc_long_set_label <- function(x) {
  if (length(x) == 1L) paste0("=", x) else paste0(" in {", paste(x, collapse = ","), "}")
}
