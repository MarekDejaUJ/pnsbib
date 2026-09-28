test_that("the model rejects malformed and inconsistent margins", {
  f <- fixture_from_response_types()
  expect_s3_class(poc_model(f$o, f$a), "pnsbib_model")
  bad <- f$a
  bad[1, 1] <- bad[1, 1] + 0.1
  expect_error(poc_model(f$o, bad), "sum to one")
  bad <- f$a
  bad[1, ] <- c(1, 0, 0)
  expect_error(poc_model(f$o, bad), "consistency")
  expect_error(poc_query(c(x1 = "y1", x1 = "y2")), "distinct")
})

test_that("the long table and matrix contracts agree", {
  f <- fixture_from_response_types()
  d <- expand.grid(treatment = rownames(f$o), outcome = colnames(f$o),
                   stringsAsFactors = FALSE)
  d$observed_joint <- mapply(function(x, y) f$o[x, y],
                             d$treatment, d$outcome)
  d$interventional <- mapply(function(x, y) f$a[x, y],
                             d$treatment, d$outcome)
  model <- poc_model(d)
  expect_equal(model$o[rownames(f$o), colnames(f$o)], f$o)
  expect_equal(model$a[rownames(f$a), colnames(f$a)], f$a)
})

test_that("zero-mass factual conditioning is undefined", {
  o <- matrix(c(0.5, 0, 0.5, 0), nrow = 2,
              dimnames = list(c("x0", "x1"), c("y0", "y1")))
  a <- matrix(0.5, 2, 2, dimnames = dimnames(o))
  model <- poc_model(o, a)
  q <- poc_query(c(x0 = "y0"), observed_x = "x1",
                 observed_y = "y0", conditional = TRUE)
  result <- poc_bounds(model, q)
  expect_identical(result$status, "undefined_condition")
  expect_true(is.na(result$lower))
})
