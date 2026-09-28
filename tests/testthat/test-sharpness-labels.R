single_term_counterexample <- function() {
  o <- rbind(c(.1, .1, .3), c(.1, .125, .025), c(.1, .125, .025))
  a <- rbind(c(.55, .125, .325), c(.35, .375, .275), c(.35, .375, .275))
  dimnames(o) <- dimnames(a) <- list(paste0("x", 1:3), paste0("y", 1:3))
  poc_model(o, a)
}

test_that("printed preservation and replacement can be non-sharp", {
  model <- single_term_counterexample()
  for (k in 1:2) for (conditional in c(FALSE, TRUE)) {
    q <- poc_query(c(x1 = "y1"), observed_y = paste0("y", k), conditional = conditional)
    paper <- poc_bounds(model, q, use_zig = FALSE)
    native <- poc_bounds(model, q, use_zig = TRUE)
    exact <- poc_exact(model, q)
    den <- if (conditional) unname(model$py[k]) else 1
    expect_identical(paper$sharpness, "valid_analytical")
    expect_identical(native$sharpness, "valid_analytical")
    expect_identical(exact$sharpness, "sharp")
    expect_equal(paper$lower, c(.1, .15)[k] / den)
    expect_equal(paper$upper, c(.3, .25)[k] / den)
    expect_equal(exact$lower, c(.25, .2)[k] / den)
    expect_equal(exact$upper, c(.3, .25)[k] / den)
    expect_equal(native, paper)
    expect_gt(exact$lower, paper$lower + .01)
  }
})

test_that("independent outside-factual-stratum algebra agrees with the LP", {
  model <- single_term_counterexample()
  for (x in 1:3) for (i in 1:3) for (k in 1:3) {
    # On X=x consistency fixes Y_x. On X!=x two events with masses
    # d,c lie in a universe of mass w. Their overlap obeys Frechet bounds.
    d <- model$a[x, i] - model$o[x, i]
    c <- model$py[k] - model$o[x, k]
    w <- 1 - model$px[x]
    fixed <- if (i == k) model$o[x, i] else 0
    expected <- c(fixed + max(0, d + c - w), fixed + min(d, c))
    q <- poc_query(setNames(paste0("y", i), paste0("x", x)), observed_y = paste0("y", k))
    exact <- poc_exact(model, q)
    paper <- poc_bounds(model, q)
    expect_equal(c(exact$lower, exact$upper), expected, tolerance = 1e-9)
    expect_lte(paper$lower, exact$lower + 1e-9)
    expect_gte(paper$upper, exact$upper - 1e-9)
    expect_identical(paper$sharpness, "valid_analytical")
  }
})

test_that("sharpness labels survive summary exports and zero conditions", {
  model <- single_term_counterexample()
  q <- poc_query(c(x1 = "y1"), observed_y = "y1")
  result <- poc_bounds(model, q)
  directory <- tempfile("pnsbib-sharpness-")
  poc_export_summary(poc_analysis(list(result, poc_exact(model, q)),
    provenance = list(input = "synthetic_counterexample")), directory)
  exported <- read.delim(file.path(directory, "bounds.tsv"))
  expect_identical(exported$sharpness, c("valid_analytical", "sharp"))
  expect_equal(exported$lower, c(.1, .25))
  o <- a <- matrix(c(1, 0, 0, 0), 2, byrow = TRUE,
                    dimnames = list(c("x1", "x2"), c("y1", "y2")))
  a[2, 1] <- 1
  zero <- poc_model(o, a)
  for (required in c("y1", "y2")) {
    result <- poc_bounds(zero, poc_query(c(x1 = required), observed_y = "y2", conditional = TRUE))
    expect_identical(result$status, "undefined_condition")
    expect_identical(result$sharpness, "valid_analytical")
    expect_true(is.na(result$lower))
  }
  for (q in list(poc_query(c(x1 = "y1")),
                 poc_query(c(x1 = "y1"), observed_x = "x1"),
                 poc_query(c(x1 = "y1"), observed_x = "x2"),
                 poc_query(c(x1 = "y1"), observed_x = "x2", observed_y = "y2"))) {
    expect_identical(poc_bounds(model, q)$sharpness, "sharp")
  }
})
