band_fixture <- function() {
  labels <- c("0", "1")
  observed <- data.frame(history = rep(labels, each = 2),
                         path = rep(labels, 2), prob = c(.3, .2, .1, .4))
  intervention <- data.frame(regime = observed$history, path = observed$path,
                             prob = c(.7, .3, .4, .6))
  regimes <- matrix(0:1, ncol = 1, dimnames = list(labels, "t1"))
  long <- poc_longitudinal_model(observed, regimes, intervention)
  o <- matrix(observed$prob, 2, 2, byrow = TRUE, dimnames = list(labels, labels))
  a <- matrix(intervention$prob, 2, 2, byrow = TRUE, dimnames = dimnames(o))
  list(model = long, static = poc_model(o, a), observed = observed,
       regimes = regimes, intervention = intervention)
}

test_that("longitudinal bands reduce to the static cell-band LP at one period", {
  f <- band_fixture()
  queries <- list(pns = poc_query(c("0" = "0", "1" = "1")),
    pn = poc_query(c("0" = "0"), observed_x = "1", observed_y = "1", conditional = TRUE),
    ps = poc_query(c("1" = "1"), observed_x = "0", observed_y = "0", conditional = TRUE))
  for (kind in names(queries)) for (od in c(0, .1, 1)) for (ad in c(0, .1, 1)) {
    q <- poc_longitudinal_query(kind, "1", "0", 1)
    long <- poc_longitudinal_sensitivity(f$model, q, od, ad)
    static <- poc_input_sensitivity(f$static, queries[[kind]], od, ad)
    expect_identical(long$status, static$status)
    expect_equal(c(long$lower, long$upper), c(static$lower, static$upper), tolerance = 1e-8)
    if (kind != "pns")
      expect_equal(c(long$condition_probability_min, long$condition_probability_max),
                   c(static$condition_probability_min, static$condition_probability_max),
                   tolerance = 1e-8)
  }
  for (kind in c("pn", "ps")) {
    q <- poc_longitudinal_query(kind, "1", "0", 1)
    long <- poc_longitudinal_sensitivity(f$model, q, .3, .1,
                                         minimum_condition_probability = .4)
    static <- poc_input_sensitivity(f$static, queries[[kind]], .3, .1,
                                      minimum_condition_probability = .4)
    expect_identical(long$status, "ok")
    expect_equal(c(long$lower, long$upper), c(static$lower, static$upper), tolerance = 1e-8)
    expect_true(is.na(long$denominator))
  }
})

test_that("a fully relaxed intervention band equals absent intervention information", {
  f <- band_fixture()
  unknown <- poc_longitudinal_model(f$observed, f$regimes)
  for (kind in c("pns", "pn", "ps")) {
    q <- poc_longitudinal_query(kind, "1", "0", 1)
    wide <- poc_longitudinal_sensitivity(f$model, q, path_delta = 1)
    free <- poc_longitudinal_bounds(unknown, q)
    exact <- poc_longitudinal_bounds(f$model, q)
    zero <- poc_longitudinal_sensitivity(f$model, q)
    expect_equal(c(wide$lower, wide$upper), c(free$lower, free$upper), tolerance = 1e-9)
    expect_equal(c(zero$lower, zero$upper), c(exact$lower, exact$upper), tolerance = 1e-9)
    expect_equal(zero$denominator, exact$denominator, tolerance = 1e-9)
    no_information <- poc_longitudinal_sensitivity(unknown, q, path_delta = .6, time_delta = .8)
    expect_equal(c(no_information$lower, no_information$upper), c(free$lower, free$upper))
  }
})

test_that("time and path bands retain their independent information", {
  f <- band_fixture()
  times <- data.frame(regime = c("0", "1"), horizon = 1, prob1 = c(.3, .6))
  timed <- poc_longitudinal_model(f$observed, f$regimes, time_margins = times)
  both <- poc_longitudinal_model(f$observed, f$regimes, f$intervention, times)
  q <- poc_longitudinal_query("pns", "1", "0", 1)
  for (d in c(0, .1, 1)) {
    a <- poc_longitudinal_sensitivity(f$model, q, path_delta = d)
    b <- poc_longitudinal_sensitivity(timed, q, time_delta = d)
    expect_equal(c(a$lower, a$upper), c(b$lower, b$upper), tolerance = 1e-9)
  }
  exact <- poc_longitudinal_bounds(f$model, q)
  for (pair in list(c(1, 0), c(0, 1))) {
    bound <- poc_longitudinal_sensitivity(both, q, path_delta = pair[1], time_delta = pair[2])
    expect_equal(c(bound$lower, bound$upper), c(exact$lower, exact$upper), tolerance = 1e-9)
  }
  partial <- poc_longitudinal_model(f$observed, f$regimes,
                                    interventional_paths = f$intervention[1, ],
                                    time_margins = times[2, ])
  bound <- poc_longitudinal_sensitivity(partial, q)
  expect_equal(c(bound$lower, bound$upper), c(exact$lower, exact$upper), tolerance = 1e-9)
})

test_that("conditional zero cells can become possible only through observed bands", {
  f <- band_fixture()
  observed <- f$observed
  observed$prob <- c(.5, 0, .5, 0)
  model <- poc_longitudinal_model(observed, f$regimes)
  q <- poc_longitudinal_query("pn", "1", "0", 1)
  zero <- poc_longitudinal_sensitivity(model, q, path_delta = 1)
  expect_identical(zero$status, "undefined_condition")
  expect_equal(zero$denominator, 0)
  possible <- poc_longitudinal_sensitivity(model, q, observed_delta = .1)
  expect_identical(possible$status, "ok")
  expect_equal(c(possible$condition_probability_min, possible$condition_probability_max), c(0, .1))
  expect_equal(c(possible$lower, possible$upper), c(0, 1))
  expect_true(is.na(possible$denominator))
  impossible <- poc_longitudinal_sensitivity(model, q, observed_delta = .1,
                                              minimum_condition_probability = .2)
  expect_identical(impossible$status, "undefined_condition")
  expect_true(is.na(impossible$lower))
})

test_that("path restrictions and no anticipation survive even maximum cell bands", {
  regimes <- rbind(early = c(0, 0), delayed = c(0, 1))
  observed <- expand.grid(history = rownames(regimes), path = c("00", "01"),
                          stringsAsFactors = FALSE)
  observed$prob <- .25
  model <- poc_longitudinal_model(observed, regimes, absorbing = TRUE,
                                  allowed_paths = c("00", "01"))
  q <- poc_longitudinal_query("pns", "delayed", "early", 1)
  bound <- poc_longitudinal_sensitivity(model, q, 1, 1, 1)
  expect_equal(c(bound$lower, bound$upper), c(0, 0))
  observed <- expand.grid(history = rownames(regimes), path = c("00", "10", "01", "11"),
                          stringsAsFactors = FALSE)
  observed$prob <- 1 / 8
  shared <- poc_longitudinal_model(observed, regimes)
  free <- poc_longitudinal_model(observed, regimes, no_anticipation = FALSE)
  bound <- poc_longitudinal_sensitivity(shared, q, 1, 1, 1)
  expect_equal(c(bound$lower, bound$upper), c(0, 0))
  expect_equal(poc_longitudinal_sensitivity(free, q, 1, 1, 1)$upper, 1)
})

test_that("mixed latent populations give nested intervals containing known query truth", {
  set.seed(20260927)
  for (absorbing in c(FALSE, TRUE)) {
    labels <- if (absorbing) c("00", "01", "11") else c("00", "10", "01", "11")
    latent <- expand.grid(h = c("no", "yes"), p0 = labels, p1 = labels,
                           stringsAsFactors = FALSE)
    latent$mass <- runif(nrow(latent))
    latent$mass <- latent$mass / sum(latent$mass)
    factual <- ifelse(latent$h == "yes", latent$p1, latent$p0)
    observed <- expand.grid(history = c("no", "yes"), path = labels, stringsAsFactors = FALSE)
    observed$prob <- vapply(seq_len(nrow(observed)), function(i)
      sum(latent$mass[latent$h == observed$history[i] & factual == observed$path[i]]), 0)
    margins <- data.frame(regime = observed$history, path = observed$path)
    margins$prob <- vapply(seq_len(nrow(margins)), function(i)
      sum(latent$mass[(if (margins$regime[i] == "yes") latent$p1 else latent$p0) == margins$path[i]]), 0)
    model <- poc_longitudinal_model(observed, rbind(no = c(0, 0), yes = c(1, 1)),
                                    margins, absorbing = absorbing)
    y1 <- substr(latent$p1, 2, 2) == "1"
    y0 <- substr(latent$p0, 2, 2) == "1"
    truth <- c(pns = sum(latent$mass[y1 & !y0]),
      pn = sum(latent$mass[latent$h == "yes" & y1 & !y0]) / sum(latent$mass[latent$h == "yes" & y1]),
      ps = sum(latent$mass[latent$h == "no" & !y0 & y1]) / sum(latent$mass[latent$h == "no" & !y0]),
      trajectory = sum(latent$mass[latent$p1 == "11" & latent$p0 == "00"]),
      persistent = sum(latent$mass[latent$p1 == "11" & latent$p0 == "00"]),
      event_time = sum(latent$mass[latent$p1 == "01" & !y0]),
      lagged = sum(latent$mass[y1 & substr(latent$p0, 1, 1) == "0"]))
    for (kind in names(truth)) {
      q <- poc_longitudinal_query(kind, "yes", "no", 2,
                                  regime_path = "11", reference_path = "00")
      bounds <- lapply(c(0, .05, .2, 1), function(d)
        poc_longitudinal_sensitivity(model, q, d, d, d))
      exact <- poc_longitudinal_bounds(model, q)
      expect_equal(c(bounds[[1]]$lower, bounds[[1]]$upper), c(exact$lower, exact$upper), tolerance = 1e-8)
      for (b in bounds) {
        expect_identical(b$status, "ok")
        expect_lte(b$lower, unname(truth[kind]) + 1e-8)
        expect_gte(b$upper, unname(truth[kind]) - 1e-8)
      }
      expect_true(all(diff(vapply(bounds, function(b) b$lower, 0)) <= 1e-8))
      expect_true(all(diff(vapply(bounds, function(b) b$upper, 0)) >= -1e-8))
    }
  }
})

test_that("bands validate inputs and keep uncertainty labels through export", {
  f <- band_fixture()
  q <- poc_longitudinal_query("pn", "1", "0", 1)
  for (bad in list(-.1, 1.1, Inf, NA_real_, numeric(), c(0, .1), "0")) {
    expect_error(poc_longitudinal_sensitivity(f$model, q, observed_delta = bad), "observed_delta")
    expect_error(poc_longitudinal_sensitivity(f$model, q, path_delta = bad), "path_delta")
    expect_error(poc_longitudinal_sensitivity(f$model, q, time_delta = bad), "time_delta")
  }
  for (bad in list(0, -.1, 1.1, Inf, NA_real_, c(.1, .2)))
    expect_error(poc_longitudinal_sensitivity(f$model, q, minimum_condition_probability = bad),
                 "minimum_condition_probability")
  expect_error(poc_longitudinal_sensitivity(f$model, poc_longitudinal_query("pns", "1", "0", 1),
                                           minimum_condition_probability = .1), "PN or PS")
  result <- poc_longitudinal_sensitivity(f$model, q, .1, .1)
  object <- poc_analysis(result, provenance = list(source = "synthetic_band_fixture", seed = 20260927))
  directory <- tempfile("long-band-export-")
  poc_export_summary(object, directory)
  exported <- utils::read.delim(file.path(directory, "bounds.tsv"))
  expect_equal(exported$lower, result$lower)
  expect_true(is.na(exported$denominator))
  expect_equal(exported$condition_probability_max, result$condition_probability_max)
  expect_identical(exported$sharpness, "sharp_given_bands")
  expect_match(exported$assumptions, "joint_absolute_longitudinal_cell_bands")
})
