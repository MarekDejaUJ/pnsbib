backend_fixture <- function(backend) {
  o <- rbind(`0` = c(`0` = .3, `1` = .2), `1` = c(`0` = .1, `1` = .4))
  observed <- expand.grid(history = c("off", "on"), path = c("00", "01", "11"),
    stringsAsFactors = FALSE)
  observed$prob <- c(.2, .1, .1, .1, .2, .3)
  list(static = poc_model(o, 2 * o, lp_backend = backend),
    partial = poc_partial_model(o, lp_backend = backend),
    longitudinal = poc_longitudinal_model(observed,
      rbind(off = c(0, 0), on = c(1, 1)), absorbing = TRUE, lp_backend = backend))
}

backend_results <- function(models) {
  query <- poc_query(c(`0` = "0", `1` = "1"))
  conditional <- poc_query(c(`0` = "0"), "1", "1", conditional = TRUE)
  long_query <- poc_longitudinal_query("pn", "on", "off", 1, condition_on = "prefix")
  map <- poc_margin_map(models$static)
  samples <- setNames(rep(list(c(0, 1, 0, 1)), nrow(map)), map$cell_id)
  region <- poc_confidence_region(models$static, samples, target = "test population",
    sampling_unit = "individual", assume_sampling = TRUE)
  list(exact = poc_exact(models$static, query),
    exact_conditional = poc_exact(models$static, conditional),
    partial = poc_partial_bounds(models$partial, query),
    partial_conditional = poc_partial_bounds(models$partial, conditional),
    sensitivity = poc_sensitivity(models$static, conditional),
    joint_bands = poc_input_sensitivity(models$static, conditional, .02, .03),
    longitudinal = poc_longitudinal_bounds(models$longitudinal, long_query),
    longitudinal_bands = poc_longitudinal_sensitivity(models$longitudinal, long_query, .02),
    curve = poc_curve(models$longitudinal, "on", "off"),
    curve_bands = poc_curve(models$longitudinal, "on", "off", bands = list(observed_delta = .02)),
    confidence = poc_confidence_bounds(region, query),
    confidence_conditional = poc_confidence_bounds(region, conditional))
}

test_that("all public continuous LP routes can use Zig without the R solver", {
  reference <- backend_results(backend_fixture("reference"))
  calls <- 0L
  kernel <- poc_lp_zig
  testthat::local_mocked_bindings(lp = function(...) stop("R solver must not be called"),
    .package = "lpSolve")
  testthat::local_mocked_bindings(poc_lp_zig = function(...) {
    calls <<- calls + 1L
    kernel(...)
  }, .package = "pnsbib")
  models <- backend_fixture("zig")
  expect_true(all(vapply(models, function(m) identical(m$lp_backend, "zig"), logical(1))))
  native <- backend_results(models)
  for (name in names(reference)) expect_equal(native[[name]], reference[[name]], tolerance = 1e-8)
  expect_gte(calls, 50L)
})

test_that("LP choice is explicit, inherited and compatible with old models", {
  models <- backend_fixture("reference")
  q <- poc_query(c(`0` = "0", `1` = "1"))
  expect_equal(poc_exact(models$static, q, lp_backend = "zig"), poc_exact(models$static, q), tolerance = 1e-8)
  expect_equal(poc_partial_bounds(models$partial, q, lp_backend = "zig"),
    poc_partial_bounds(models$partial, q), tolerance = 1e-8)
  expect_equal(poc_curve(models$longitudinal, "on", "off", lp_backend = "zig"),
    poc_curve(models$longitudinal, "on", "off"), tolerance = 1e-8)
  models$static$lp_backend <- NULL
  expect_identical(pnsbib:::.poc_lp_backend(NULL, models$static), "reference")
  expect_equal(poc_exact(models$static, q), poc_exact(models$static, q, lp_backend = "zig"), tolerance = 1e-8)
  for (bad in list("auto", "ZIG", TRUE, NA_character_, character(), c("zig", "reference"))) {
    expect_error(poc_exact(models$static, q, lp_backend = bad), "lp_backend")
    expect_error(poc_partial_bounds(models$partial, q, lp_backend = bad), "lp_backend")
    expect_error(poc_curve(models$longitudinal, "on", "off", lp_backend = bad), "lp_backend")
  }
})

test_that("native endpoint witnesses satisfy the supplied constraints", {
  model <- backend_fixture("zig")$partial
  q <- poc_query(c(`0` = "0", `1` = "1"))
  fit <- poc_partial_bounds(model, q, return_witness = TRUE)
  expect_identical(fit$status, "ok")
  witness <- attr(fit, "witnesses")
  expect_length(witness, 2L)
  system <- pnsbib:::.poc_partial_system(model)
  event <- pnsbib:::.poc_partial_event(model, q, NULL, system)
  for (i in 1:2) {
    p <- witness[[i]]$prob
    expect_equal(drop(system$A %*% p), system$b, tolerance = 1e-8)
    expect_gte(min(p), -1e-8)
    expect_equal(sum(p * event$event), c(fit$lower, fit$upper)[i], tolerance = 1e-8)
  }
})

test_that("LP computational failures cannot masquerade as incompatible margins", {
  lp <- pnsbib:::.poc_lp
  expect_error(lp("max", 1, matrix(numeric(), 0, 1), character(), numeric(), "zig"),
    "LP solver.*status 3")
  impossible <- lp("min", 0, matrix(1, 1, 1), "=", -1, "zig")
  expect_identical(impossible$status, 2L)
  testthat::local_mocked_bindings(poc_lp_zig = function(...) stop("IterationLimit"),
    .package = "pnsbib")
  expect_error(lp("min", 0, matrix(1, 1, 1), "=", 1, "zig"), "IterationLimit")
  testthat::local_mocked_bindings(lp = function(...) list(status = 5L), .package = "lpSolve")
  expect_error(lp("min", 0, matrix(1, 1, 1), "=", 1, "reference"), "LP solver.*status 5")
})
