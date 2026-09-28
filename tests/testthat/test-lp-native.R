native_lp_fixture <- function(direction, objective, A, dirs, b, iterations = 200000L) {
  raw <- poc_lp_zig(nrow(A), ncol(A), as.numeric(A), as.numeric(b),
    as.integer(match(dirs, c("<=", "=", ">=")) - 2L), as.numeric(objective),
    direction == "max", iterations)
  n <- ncol(A); m <- nrow(A)
  list(status = raw[1L], objval = raw[2L], iterations = raw[3L],
    primal_violation = raw[4L], dual_violation = raw[5L], gap = raw[6L],
    solution = raw[6L + seq_len(n)], dual = raw[6L + n + seq_len(m)],
    ray = raw[6L + n + m + seq_len(n)])
}

check_native_lp_witness <- function(fit, objective, A, dirs, b, direction) {
  expect_equal(fit$status, 0)
  expect_true(all(is.finite(fit$solution)))
  expect_true(all(fit$solution >= -1e-8))
  value <- drop(A %*% fit$solution)
  for (i in seq_along(b)) {
    if (dirs[i] == "=") expect_equal(value[i], b[i], tolerance = 1e-7)
    if (dirs[i] == "<=") expect_lte(value[i], b[i] + 1e-7)
    if (dirs[i] == ">=") expect_gte(value[i], b[i] - 1e-7)
  }
  sign <- if (direction == "max") 1 else -1
  expect_true(all(drop(crossprod(A, sign * fit$dual)) >= sign * objective - 1e-7))
  expect_true(all(sign * fit$dual[dirs == "<="] >= -1e-7))
  expect_true(all(sign * fit$dual[dirs == ">="] <= 1e-7))
  expect_equal(sum(b * fit$dual), fit$objval, tolerance = 1e-7)
  expect_equal(sum(objective * fit$solution), fit$objval, tolerance = 1e-7)
}

test_that("native continuous LP matches independent seeded generic programs", {
  set.seed(20261002)
  for (n in c(1L, 2L, 3L, 5L, 8L)) for (replicate in 1:8) {
    x <- if (replicate %% 2L) runif(n) else sample(c(0, 0, 1, 2), n, replace = TRUE)
    A <- matrix(sample(-3:3, 4L * n, replace = TRUE), 4L, n)
    dirs <- c("=", "=", "<=", ">=")
    b <- drop(A %*% x) + c(0, 0, runif(1), -runif(1))
    A <- rbind(A, diag(n), -A[1L, , drop = FALSE])
    b <- c(b, x + 1, -b[1L]); dirs <- c(dirs, rep("<=", n), "=")
    objective <- runif(n, -1, 1)
    saved <- list(objective, A, dirs, b)
    for (direction in c("min", "max")) {
      fit <- native_lp_fixture(direction, objective, A, dirs, b)
      reference <- lpSolve::lp(direction, objective, A, dirs, b)
      expect_equal(reference$status, 0L)
      expect_equal(fit$objval, reference$objval, tolerance = 1e-7)
      check_native_lp_witness(fit, objective, A, dirs, b, direction)
    }
    expect_identical(list(objective, A, dirs, b), saved)
  }
})

test_that("native LP preserves complete response-type systems", {
  for (shape in list(c(2L, 2L), c(3L, 3L), c(4L, 3L))) {
    f <- fixture_from_response_types(shape[1L], shape[2L], seed = 89L)
    model <- poc_model(f$o, f$a, check_lp = FALSE)
    system <- pnsbib:::.poc_lp_system(model)
    dirs <- rep("=", nrow(system$A))
    objective <- as.numeric(system$response[, 1L] == 1L & system$response[, 2L] == 2L)
    for (direction in c("min", "max")) {
      fit <- native_lp_fixture(direction, objective, system$A, dirs, system$b)
      reference <- lpSolve::lp(direction, objective, system$A, dirs, system$b)
      expect_equal(fit$objval, reference$objval, tolerance = 1e-8)
      check_native_lp_witness(fit, objective, system$A, dirs, system$b, direction)
    }
  }
})

test_that("native LP failure states carry independently checked certificates", {
  fit <- native_lp_fixture("max", 0, matrix(c(1, 1), 2L), c("<=", ">="), c(0, 1))
  expect_equal(fit$status, 2)
  expect_true(is.nan(fit$objval))
  expect_gte(sum(fit$dual), -1e-8)
  expect_lt(fit$dual[2L], -1e-8)
  fit <- native_lp_fixture("max", c(1, 0), matrix(c(1, -1), 1L), ">=", 0)
  expect_equal(fit$status, 3)
  expect_true(all(fit$ray >= -1e-8))
  expect_gt(fit$ray[1L], 0)
  expect_gte(fit$ray[1L] - fit$ray[2L], -1e-8)
})

test_that("native LP rejects malformed or exhausted calls without poisoning R", {
  run <- function(m = 1L, n = 1L, A = 1, b = 1, dirs = 0L, objective = 1, iterations = 200000L)
    poc_lp_zig(m, n, A, b, dirs, objective, TRUE, iterations)
  expect_error(run(n = 2L), "InvalidDimensions")
  expect_error(run(A = Inf), "(InvalidInput|finite)")
  expect_error(run(dirs = 2L), "InvalidInput")
  expect_error(run(iterations = 0L), "InvalidOptions")
  expect_error(run(m = 10000000L), "WorkspaceLimit")
  expect_error(native_lp_fixture("max", c(3, 2), matrix(c(1, 1, 0, 1, 0, 1), 3L),
    rep("<=", 3L), c(4, 2, 3), iterations = 1L), "IterationLimit")
  expect_equal(run()[1:2], c(0, 1), tolerance = 1e-9)
})
