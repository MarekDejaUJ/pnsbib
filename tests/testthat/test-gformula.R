baseline_fixture <- function() {
  cells <- data.frame(
    L = c(0, 0, 0, 0, 1, 1, 1, 1),
    X = c("x0", "x0", "x1", "x1", "x0", "x0", "x1", "x1"),
    Y = c("y0", "y1", "y0", "y1", "y0", "y1", "y0", "y1"),
    n = c(18, 2, 7, 3, 5, 5, 2, 18)
  )
  cells[rep(seq_len(nrow(cells)), cells$n), c("L", "X", "Y")]
}

test_that("baseline standardization uses one target covariate distribution", {
  d <- baseline_fixture()
  fit <- poc_gformula(d, "X", "Y", "L", assume_identification = TRUE)
  expected_o <- matrix(c(23, 9, 7, 21) / 60, 2, 2,
                       dimnames = list(c("x0", "x1"), c("y0", "y1")))
  expected_a <- matrix(c(0.7, 0.4, 0.3, 0.6), 2, 2,
                       dimnames = dimnames(expected_o))
  expect_equal(fit$model$o, expected_o)
  expect_equal(fit$model$a, expected_a)
  expect_true(fit$model$checked_lp)
  expect_equal(fit$diagnostics$minimum_propensity, 1 / 3)
  expect_equal(fit$diagnostics$effective_n, 60)
  expect_equal(nrow(fit$support), 4L)
})

test_that("weights change the target population without hiding support", {
  d <- baseline_fixture()
  d$w <- ifelse(d$L == 1, 2, 1)
  fit <- poc_gformula(d, "X", "Y", "L", weights = "w",
                      assume_identification = TRUE)
  expect_equal(fit$model$a["x0", "y1"], 11 / 30)
  expect_equal(fit$model$a["x1", "y1"], 7 / 10)
  expect_equal(fit$model$o["x1", "y1"], 39 / 90)
  expect_equal(fit$diagnostics$strata, 2L)
})

test_that("an unsupported treatment and undeclared design stop estimation", {
  d <- baseline_fixture()
  expect_error(poc_gformula(d, "X", "Y", "L"),
               "assume_identification")
  d <- d[!(d$L == 1 & d$X == "x1"), ]
  expect_error(poc_gformula(d, "X", "Y", "L",
                            assume_identification = TRUE), "Positivity fails")
  d <- baseline_fixture()
  d$X <- factor(d$X, levels = c("x0", "x1", "x2"))
  expect_error(poc_gformula(d, "X", "Y", "L",
                            assume_identification = TRUE), "Positivity fails")
})
