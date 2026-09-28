test_that("margin bands contain the fixed-margin exact interval", {
  f <- fixture_from_response_types(seed = 93L)
  model <- poc_model(f$o, f$a)
  query <- poc_query(c(x1 = "y1", x2 = "y3"))
  exact <- poc_exact(model, query)
  sensitivity <- poc_sensitivity(model, query, delta = c(0, 0.03, 0.1))
  expect_true(all(sensitivity$status == "ok"))
  expect_equal(sensitivity$lower[1], exact$lower, tolerance = 1e-8)
  expect_equal(sensitivity$upper[1], exact$upper, tolerance = 1e-8)
  expect_true(all(diff(sensitivity$lower) <= 1e-8))
  expect_true(all(diff(sensitivity$upper) >= -1e-8))
  expect_error(poc_sensitivity(model, query, delta = -0.01), "delta")
})

test_that("joint input bands agree with exact LP at zero width and widen", {
  f <- fixture_from_response_types(seed = 107L)
  model <- poc_model(f$o, f$a)
  query <- poc_query(c(x1 = "y1", x2 = "y3"))
  exact <- poc_exact(model, query)
  fixed <- poc_input_sensitivity(model, query)
  wider <- poc_input_sensitivity(model, query, observed_delta = 0.01,
                                 interventional_delta = 0.02)
  expect_equal(fixed$status, "ok")
  expect_equal(fixed$lower, exact$lower, tolerance = 1e-8)
  expect_equal(fixed$upper, exact$upper, tolerance = 1e-8)
  expect_equal(wider$status, "ok")
  expect_lte(wider$lower, fixed$lower + 1e-8)
  expect_gte(wider$upper, fixed$upper - 1e-8)
  expect_error(poc_input_sensitivity(model, query, observed_delta = -0.1),
               "observed_delta")
  expect_error(poc_input_sensitivity(model, query,
                                      interventional_delta = Inf),
               "interventional_delta")
})

test_that("conditional input bands use the varying factual denominator", {
  f <- fixture_from_response_types(seed = 108L)
  model <- poc_model(f$o, f$a)
  queries <- list(
    poc_query(c(x2 = "y3"), observed_x = "x1", observed_y = "y1",
              conditional = TRUE),
    poc_query(c(x2 = "y3"), observed_x = "x1", conditional = TRUE),
    poc_query(c(x2 = "y3"), observed_y = "y1", conditional = TRUE)
  )
  for (query in queries) {
    exact <- poc_exact(model, query)
    fixed <- poc_input_sensitivity(model, query)
    expect_equal(fixed$status, "ok")
    expect_equal(fixed$lower, exact$lower, tolerance = 1e-8)
    expect_equal(fixed$upper, exact$upper, tolerance = 1e-8)
    expect_equal(fixed$condition_probability_min,
                 fixed$condition_probability_max, tolerance = 1e-8)
  }
  query <- queries[[1L]]
  fixed <- poc_input_sensitivity(model, query)
  wider <- poc_input_sensitivity(model, query, observed_delta = 0.005,
                                 interventional_delta = 0.01)
  expect_equal(wider$status, "ok")
  expect_lte(wider$lower, fixed$lower + 1e-8)
  expect_gte(wider$upper, fixed$upper - 1e-8)
  expect_equal(poc_input_sensitivity(
    model, query, minimum_condition_probability = 0.9)$status,
    "undefined_condition")
  expect_error(poc_input_sensitivity(
    model, queries[[2L]], minimum_condition_probability = 0),
    "minimum_condition_probability")
})
