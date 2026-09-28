test_that("native static matrices and every event match the preserved R oracle", {
  for (shape in list(c(2L, 2L), c(2L, 3L), c(3L, 2L), c(3L, 3L), c(4L, 2L), c(4L, 3L))) {
    m <- shape[1L]; n <- shape[2L]
    f <- fixture_from_response_types(m, n, seed = 47L + 10L * m + n)
    model <- poc_model(f$o, f$a, check_lp = FALSE)
    before <- serialize(model, NULL)
    native <- pnsbib:::.poc_lp_system(model, lp_backend = "zig")
    oracle <- pnsbib:::.poc_lp_system_reference(model)
    expect_identical(native, oracle)
    expect_identical(serialize(model, NULL), before)
    grid <- as.matrix(expand.grid(rep(list(0:n), m)))
    count <- 0L
    for (i in seq_len(nrow(grid))) for (ox in 0:m) for (oy in 0:n) {
      x <- which(grid[i, ] != 0L)
      idx <- list(x = x, y = as.integer(grid[i, x]),
        ox = if (ox) ox else NA_integer_, oy = if (oy) oy else NA_integer_)
      want <- pnsbib:::.poc_event_indicator_reference(oracle, idx)
      got <- pnsbib:::.poc_event_indicator(native, idx, lp_backend = "zig")
      expect_identical(got, want)
      count <- count + 1L
    }
    expect_equal(count, (n + 1)^m * (m + 1L) * (n + 1L))
  }
})

test_that("public static native calls do not invoke numerical R constructors", {
  f <- fixture_from_response_types(seed = 135L)
  model <- poc_model(f$o, f$a, lp_backend = "reference")
  q <- poc_query(c(x1 = "y1", x2 = "y2"))
  conditional <- poc_query(c(x1 = "y1"), observed_x = "x2", observed_y = "y2", conditional = TRUE)
  expected <- poc_exact(model, q)
  band <- poc_input_sensitivity(model, conditional, .01, .01)
  testthat::local_mocked_bindings(.poc_lp_system_reference = function(...) stop("R construction"),
    .poc_event_indicator_reference = function(...) stop("R event mask"), .package = "pnsbib")
  testthat::local_mocked_bindings(lp = function(...) stop("R LP"), .package = "lpSolve")
  native <- poc_model(f$o, f$a, lp_backend = "zig")
  expect_equal(poc_exact(native, q), expected, tolerance = 1e-8)
  expect_equal(poc_exact(model, q, lp_backend = "zig"), expected, tolerance = 1e-8)
  expect_equal(poc_input_sensitivity(native, conditional, .01, .01), band, tolerance = 1e-8)
  expect_error(poc_exact(native, q, max_variables = 10L), "(variables|VariableLimit)")
})

test_that("native constructor and event input failures are recoverable", {
  build <- function(m = 2L, n = 2L, o = c(.3, .1, .2, .4), a = c(.6, .2, .4, .8), cap = 10000L)
    poc_static_system_zig(m, n, o, a, cap)
  expect_error(build(m = 1L), "InvalidDimensions")
  expect_error(build(o = 1), "InvalidDimensions")
  expect_error(build(a = rep(Inf, 4)), "(InvalidInput|finite)")
  expect_error(build(cap = 7L), "VariableLimit")
  expect_error(build(m = 2L, n = 1000L, cap = 10000000L), "WorkspaceLimit")
  expect_length(build(), 116L)
  r <- c(1L, 2L, 2L, 1L); fx <- c(1L, 2L); fy <- c(1L, 1L)
  before <- serialize(list(r, fx, fy), NULL)
  event <- function(x = 1L, y = 2L, ox = 0L, oy = 0L, response = r, factual_x = fx, factual_y = fy)
    poc_static_event_zig(2L, 2L, response, factual_x, factual_y, x, y, ox, oy)
  expect_identical(event(), c(0, 1))
  expect_error(event(x = c(1L, 1L), y = c(1L, 2L)), "InvalidQuery")
  expect_error(event(x = 0L), "InvalidQuery")
  expect_error(event(x = NA_integer_), "(InvalidQuery|missing|NA)")
  expect_error(event(ox = 3L), "InvalidQuery")
  expect_error(event(oy = -1L), "InvalidQuery")
  expect_error(event(response = r[-1L]), "InvalidDimensions")
  expect_error(event(factual_y = c(2L, 1L)), "InvalidInput")
  expect_identical(serialize(list(r, fx, fy), NULL), before)
})
