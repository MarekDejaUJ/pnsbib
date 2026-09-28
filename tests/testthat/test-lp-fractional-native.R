test_that("native LP resolves the retained degenerate prefix-conditioning programs", {
  # Reconstruct the independent seed-20260927 population used by the original
  # prefix test, without relying on a local RDS failure artifact.
  set.seed(20260927)
  regimes <- rbind(off = c("00", "00"), later = c("00", "11"),
    early = c("11", "00"), always = c("11", "11"))
  n <- 80L
  potential <- array(0L, c(n, 4L, 2L))
  initial <- matrix(sample(0:1, n * 2L, replace = TRUE), n, 2L)
  for (g in 1:4) {
    potential[, g, 1L] <- initial[, if (g <= 2L) 1L else 2L]
    potential[, g, 2L] <- sample(0:1, n, replace = TRUE)
  }
  factual <- rep(1:4, length.out = n)
  weights <- runif(n); weights <- weights / sum(weights)
  labels <- c("00", "10", "01", "11")
  paths <- sapply(1:4, function(g) paste0(potential[, g, 1L], potential[, g, 2L]))
  actual <- paths[cbind(seq_len(n), factual)]
  observed <- expand.grid(history = rownames(regimes), path = labels, stringsAsFactors = FALSE)
  observed$prob <- vapply(seq_len(nrow(observed)), function(i)
    sum(weights[factual == match(observed$history[i], rownames(regimes)) & actual == observed$path[i]]), 0)
  margins <- expand.grid(regime = rownames(regimes), path = labels, stringsAsFactors = FALSE)
  margins$prob <- vapply(seq_len(nrow(margins)), function(i)
    sum(weights[paths[, match(margins$regime[i], rownames(regimes))] == margins$path[i]]), 0)
  model <- poc_longitudinal_model(observed, regimes, margins)
  for (horizon in 1:2) for (kind in c("pn", "ps")) {
    query <- poc_longitudinal_query(kind, "early", "off", horizon, condition_on = "prefix")
    event <- pnsbib:::.poc_long_event(model, query)
    cells <- event$system$A[-1L, , drop = FALSE]
    centers <- event$system$b[-1L]
    widths <- c(rep(.01, length(model$observed)), rep(.02, nrow(model$interventional_paths)))
    rhs <- c(1, pmax(centers - widths, 0), pmin(centers + widths, 1))
    A <- rbind(event$system$A[1L, , drop = FALSE], cells, cells)
    A <- rbind(cbind(A, -rhs), c(event$condition, 0))
    dirs <- c("=", rep(">=", nrow(cells)), rep("<=", nrow(cells)), "=")
    b <- c(rep(0, length(rhs)), 1)
    objective <- c(event$event, 0)
    expect_identical(dim(A), c(66L, 257L))
    for (direction in c("min", "max")) {
      reference <- lpSolve::lp(direction, objective, A, dirs, b)
      z <- poc_lp_zig(nrow(A), ncol(A), as.numeric(A), as.numeric(b),
        as.integer(match(dirs, c("<=", "=", ">=")) - 2L),
        as.numeric(objective), direction == "max", 200000L)
      expect_equal(z[1L], 0)
      expect_equal(z[2L], reference$objval, tolerance = 1e-8)
      x <- z[6L + seq_len(ncol(A))]
      expect_true(all(x >= -1e-8))
      residual <- drop(A %*% x) - b
      expect_true(all(abs(residual[dirs == "="]) < 1e-8))
      expect_true(all(residual[dirs == "<="] < 1e-8))
      expect_true(all(residual[dirs == ">="] > -1e-8))
      expect_equal(sum(objective * x), z[2L], tolerance = 1e-8)
    }
  }
})
