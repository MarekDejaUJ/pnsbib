#' Validate static probabilities-of-causation inputs
#'
#' The observational matrix contains joint probabilities P(X = x, Y = y).
#' Each row of the interventional matrix contains P(Y_x = y). Row and column
#' names must agree. The optional feasibility check tests whether one finite
#' response-type distribution can generate both matrices.
#'
#' @param observed_joint A named numeric matrix, or a long data frame with
#'   treatment, outcome, observed_joint, and interventional columns.
#' @param interventional A named numeric matrix. Omit when observed_joint is a
#'   long data frame.
#' @param tolerance Probability tolerance.
#' @param check_lp Check full response-type feasibility when tractable.
#' @param max_variables Maximum number of response-type LP variables.
#' @param lp_backend Continuous LP solver: "reference" (lpSolve, the current
#'   default) or "zig" (native rzig solver). Stored in the model and inherited
#'   by LP bounds and sensitivity calls. No automatic fallback is performed.
#' @return A pnsbib_model object.
#' @export
validate_poc_inputs <- function(observed_joint, interventional = NULL,
                                tolerance = 1e-9, check_lp = TRUE,
                                max_variables = 10000L, lp_backend = "reference") {
  lp_backend <- .poc_lp_backend(lp_backend)
  if (is.data.frame(observed_joint)) {
    d <- observed_joint
    required <- c("treatment", "outcome", "observed_joint", "interventional")
    if (!all(required %in% names(d)) || !is.null(interventional)) {
      stop("The long table needs treatment, outcome, observed_joint and interventional columns.", call. = FALSE)
    }
    xs <- unique(as.character(d$treatment))
    ys <- unique(as.character(d$outcome))
    if (anyNA(xs) || anyNA(ys) || any(!nzchar(xs)) || any(!nzchar(ys)) ||
        nrow(d) != length(xs) * length(ys) ||
        anyDuplicated(paste(d$treatment, d$outcome, sep = "\r"))) {
      stop("The long table must contain one row for every treatment/outcome pair.", call. = FALSE)
    }
    observed <- matrix(NA_real_, length(xs), length(ys), dimnames = list(xs, ys))
    intervention <- observed
    for (i in seq_len(nrow(d))) {
      x <- as.character(d$treatment[i])
      y <- as.character(d$outcome[i])
      observed[x, y] <- d$observed_joint[i]
      intervention[x, y] <- d$interventional[i]
    }
  } else {
    observed <- observed_joint
    intervention <- interventional
  }
  if (!is.matrix(observed) || !is.matrix(intervention) ||
      !is.numeric(observed) || !is.numeric(intervention) ||
      !identical(dim(observed), dim(intervention)) ||
      nrow(observed) < 2L || ncol(observed) < 2L ||
      is.null(rownames(observed)) || is.null(colnames(observed)) ||
      !identical(dimnames(observed), dimnames(intervention))) {
    stop("Inputs must be conformable named numeric matrices with at least two levels each.", call. = FALSE)
  }
  if (anyNA(dimnames(observed)[[1L]]) || anyNA(dimnames(observed)[[2L]]) ||
      any(!nzchar(rownames(observed))) || any(!nzchar(colnames(observed))) ||
      anyDuplicated(rownames(observed)) || anyDuplicated(colnames(observed))) {
    stop("Treatment and outcome labels must be unique and nonempty.", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      length(max_variables) != 1L || !is.finite(max_variables) ||
      max_variables < 1L) {
    stop("Invalid tolerance or max_variables.", call. = FALSE)
  }
  if (any(!is.finite(observed)) || any(!is.finite(intervention)) ||
      any(observed < -tolerance | observed > 1 + tolerance) ||
      any(intervention < -tolerance | intervention > 1 + tolerance)) {
    stop("Probabilities must be finite and in [0, 1].", call. = FALSE)
  }
  if (abs(sum(observed) - 1) > tolerance ||
      any(abs(rowSums(intervention) - 1) > tolerance)) {
    stop("The observational distribution and every interventional row must sum to one.", call. = FALSE)
  }
  px <- rowSums(observed)
  if (any(observed - intervention > tolerance) ||
      any(intervention - observed - (1 - px) > tolerance)) {
    stop("The observational and interventional margins violate consistency.", call. = FALSE)
  }
  model <- structure(
    list(o = observed, a = intervention, px = px,
         py = colSums(observed), tolerance = tolerance,
         checked_lp = FALSE, lp_backend = lp_backend),
    class = "pnsbib_model"
  )
  nvar <- nrow(observed) * ncol(observed)^nrow(observed)
  if (isTRUE(check_lp) && is.finite(nvar) && nvar <= max_variables) {
    system <- .poc_lp_system(model, max_variables = max_variables)
    fit <- .poc_lp("min", rep(0, system$nvar), system$A,
                    rep("=", nrow(system$A)), system$b, lp_backend)
    if (fit$status != 0L) {
      stop("No finite response-type distribution satisfies both input matrices.", call. = FALSE)
    }
    model$checked_lp <- TRUE
  }
  model
}

#' @rdname validate_poc_inputs
#' @export
poc_model <- validate_poc_inputs
