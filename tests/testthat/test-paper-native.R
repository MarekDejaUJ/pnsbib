test_that("the full native recursion agrees with the R oracle on every query shape", {
  for (shape in list(c(2L, 2L), c(2L, 3L), c(3L, 2L), c(3L, 3L), c(4L, 2L))) {
    m <- shape[1L]
    n <- shape[2L]
    f <- fixture_from_response_types(m, n, seed = 47L + 10L * m + n)
    model <- poc_model(f$o, f$a, check_lp = FALSE)
    grid <- as.matrix(expand.grid(rep(list(0:n), m)))
    max_error <- 0
    checked <- 0L
    for (i in seq_len(nrow(grid))) {
      x <- which(grid[i, ] != 0L)
      y <- as.integer(grid[i, x])
      for (ox in 0:m) for (oy in 0:n) {
        got <- poc_paper_zig(m, n, as.numeric(f$o), as.numeric(f$a),
                             as.integer(x), y, ox, oy, 1e-9)
        want <- pnsbib:::.poc_paper_interval(model, x, y,
          if (ox) ox else NA_integer_, if (oy) oy else NA_integer_, use_zig = FALSE)
        expect_length(got, 2L)
        expect_true(all(is.finite(got)))
        max_error <- max(max_error, abs(got - want))
        checked <- checked + 1L
      }
    }
    expect_equal(checked, (n + 1)^m * (m + 1L) * (n + 1L))
    expect_lte(max_error, 1e-12)
  }
})

test_that("native analytical errors are recoverable and borrowed inputs are unchanged", {
  o <- c(.1, .2, .3, .4)
  a <- c(.5, .3, .5, .7)
  saved_o <- o
  saved_a <- a
  run <- function(m = 2L, n = 2L, obs = o, inter = a, x = 1L, y = 2L,
                  ox = 0L, oy = 0L, tolerance = 1e-9)
    poc_paper_zig(m, n, obs, inter, x, y, ox, oy, tolerance)
  expect_error(run(m = 1L), "InvalidDimensions")
  expect_error(run(obs = o[-1L]), "InvalidDimensions")
  expect_error(run(x = c(1L, 1L), y = c(1L, 2L)), "InvalidQuery")
  expect_error(run(x = 0L), "InvalidQuery")
  expect_error(run(x = 3L), "InvalidQuery")
  expect_error(run(y = 3L), "InvalidQuery")
  expect_error(run(y = integer()), "InvalidQuery")
  expect_error(run(ox = -1L), "InvalidQuery")
  expect_error(run(oy = 3L), "InvalidQuery")
  expect_error(run(x = NA_integer_), "(InvalidQuery|missing|NA)")
  expect_error(run(obs = c(NA_real_, o[-1L])), "(InvalidProbability|missing|NA)")
  expect_error(run(inter = c(Inf, a[-1L])), "(InvalidProbability|finite)")
  expect_error(run(obs = c(-.1, .4, .3, .4)), "InvalidProbability")
  expect_error(run(obs = o / 2), "InvalidNormalization")
  expect_error(run(inter = rep(.1, 4)), "(InvalidNormalization|InconsistentMargins)")
  expect_error(run(inter = c(0, .3, 1, .7)), "InconsistentMargins")
  expect_error(run(tolerance = 0), "InvalidTolerance")
  expect_error(run(tolerance = NaN), "(InvalidTolerance|missing|finite|NA)")
  expect_error(run(m = 100L, n = 100L, x = 1:20, y = rep(1L, 20)), "WorkspaceLimit")
  expect_equal(run(), c(.5, .5), tolerance = 1e-12)
  expect_identical(o, saved_o)
  expect_identical(a, saved_a)
})

test_that("public native bounds retain conditional and contradiction semantics", {
  f <- fixture_from_response_types(seed = 29L)
  model <- poc_model(f$o, f$a)
  for (condition in list(list(), list(observed_x = "x3"), list(observed_y = "y1"),
                        list(observed_x = "x3", observed_y = "y1"),
                        list(observed_x = "x1", observed_y = "y2"))) {
    for (conditional in c(FALSE, TRUE)) {
      if (conditional && !length(condition)) next
      query <- do.call(poc_query, c(list(counterfactual = c(x1 = "y1", x2 = "y3"),
                                        conditional = conditional), condition))
      expect_equal(poc_bounds(model, query, use_zig = TRUE),
                   poc_bounds(model, query, use_zig = FALSE), tolerance = 1e-12)
    }
  }
})
