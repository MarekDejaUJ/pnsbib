# Explicit solver selection, with compatibility for pre-backend model objects.
# Constructors keep the reference default until the full native adoption gate.
.poc_lp_backend <- function(lp_backend = NULL, model = NULL) {
  if (is.null(lp_backend) && is.list(model)) lp_backend <- model$lp_backend
  if (is.null(lp_backend)) lp_backend <- "reference"
  if (!is.character(lp_backend) || length(lp_backend) != 1L ||
      is.na(lp_backend) || !lp_backend %in% c("reference", "zig"))
    stop('lp_backend must be "reference" or "zig" (NULL inherits the model).', call. = FALSE)
  lp_backend
}

# Both adapters solve continuous nonnegative-variable programs. Integer
# allocation references belong to analysis code and never enter this adapter.
.poc_lp <- function(direction, objective, A, directions, rhs, lp_backend) {
  lp_backend <- .poc_lp_backend(lp_backend)
  if (lp_backend == "reference") {
    fit <- lpSolve::lp(direction, objective, A, directions, rhs)
  } else {
    if (!is.character(direction) || length(direction) != 1L ||
        is.na(direction) || !direction %in% c("min", "max") ||
        !is.matrix(A) || !is.numeric(A) ||
        length(directions) != nrow(A) || anyNA(directions) ||
        any(!directions %in% c("<=", "=", ">=")))
      stop("Malformed internal LP specification.", call. = FALSE)
    m <- nrow(A); n <- ncol(A)
    z <- poc_lp_zig(m, n, as.numeric(A), as.numeric(rhs),
      as.integer(match(directions, c("<=", "=", ">=")) - 2L),
      as.numeric(objective), direction == "max", 200000L)
    fit <- list(status = as.integer(z[1L]), objval = z[2L],
      solution = z[6L + seq_len(n)], dual = z[6L + n + seq_len(m)],
      ray = z[6L + n + m + seq_len(n)], iterations = z[3L],
      primal_violation = z[4L], dual_violation = z[5L], duality_gap = z[6L])
  }
  # Every public probability objective is bounded. An unexpected unbounded
  # status or another solver failure must not become "incompatible margins".
  if (!fit$status %in% c(0L, 2L))
    stop(sprintf("LP solver (%s) returned computational status %s.", lp_backend, fit$status), call. = FALSE)
  fit
}
