test_that("one-period longitudinal LP agrees with binary static exact LP", {
  regimes <- matrix(0:1, ncol = 1, dimnames = list(c("0", "1"), "t1"))
  observed <- data.frame(history = rep(c("0", "1"), each = 2),
                         path = rep(c("0", "1"), 2),
                         prob = c(.3, .2, .1, .4))
  intervention <- data.frame(regime = rep(c("0", "1"), each = 2),
                             path = rep(c("0", "1"), 2),
                             prob = c(.7, .3, .4, .6))
  long <- poc_longitudinal_model(observed, regimes, intervention)
  o <- matrix(c(.3, .2, .1, .4), 2, 2, byrow = TRUE,
              dimnames = list(c("0", "1"), c("0", "1")))
  a <- matrix(c(.7, .3, .4, .6), 2, 2, byrow = TRUE,
              dimnames = dimnames(o))
  static <- poc_model(o, a)
  queries <- list(
    list(kind = "pns", q = poc_query(c("0" = "0", "1" = "1"))),
    list(kind = "pn", q = poc_query(c("0" = "0"), observed_x = "1",
                                    observed_y = "1", conditional = TRUE)),
    list(kind = "ps", q = poc_query(c("1" = "1"), observed_x = "0",
                                    observed_y = "0", conditional = TRUE)))
  for (item in queries) {
    l <- poc_longitudinal_bounds(long,
      poc_longitudinal_query(item$kind, "1", "0", 1))
    s <- poc_exact(static, item$q)
    expect_equal(l$lower, s$lower, tolerance = 1e-9)
    expect_equal(l$upper, s$upper, tolerance = 1e-9)
    expect_identical(l$status, "ok")
  }
})

test_that("shared treatment prefixes enforce no anticipation", {
  regimes <- rbind(early = c(0, 0), delayed = c(0, 1))
  observed <- expand.grid(history = rownames(regimes),
                          path = c("00", "10", "01", "11"),
                          stringsAsFactors = FALSE)
  observed$prob <- rep(1 / 8, nrow(observed))
  model <- poc_longitudinal_model(observed, regimes, no_anticipation = TRUE)
  result <- poc_longitudinal_bounds(model,
    poc_longitudinal_query("pns", "delayed", "early", 1))
  expect_equal(result$lower, 0)
  expect_equal(result$upper, 0)
  free <- poc_longitudinal_model(observed, regimes, no_anticipation = FALSE)
  expect_gt(poc_longitudinal_bounds(free,
    poc_longitudinal_query("pns", "delayed", "early", 1))$upper, 0)
})

test_that("joint database states can label multivalued time histories", {
  regimes <- matrix(c("00", "01", "10", "11"), ncol = 1,
                    dimnames = list(c("neither", "second", "first", "both"),
                                    "t1"))
  observed <- expand.grid(history = rownames(regimes), path = c("0", "1"),
                          stringsAsFactors = FALSE)
  observed$prob <- rep(1 / nrow(observed), nrow(observed))
  model <- poc_longitudinal_model(observed, regimes)
  expect_equal(model$n_response_types, 16)
  query <- poc_longitudinal_query("pns", "both", "neither", 1)
  bounds <- poc_longitudinal_bounds(model, query)
  expect_identical(bounds$status, "ok")
  expect_equal(bounds$lower, 0)
  expect_gt(bounds$upper, 0)
  expect_identical(as.character(model$regimes[, 1]),
                   c("00", "01", "10", "11"))
})

test_that("absorbing response types and known truth are respected", {
  regimes <- rbind(no = c(0, 0), yes = c(1, 1))
  observed <- expand.grid(history = c("no", "yes"),
                          path = c("00", "01", "11"),
                          stringsAsFactors = FALSE)
  intervention <- expand.grid(regime = c("no", "yes"),
                              path = c("00", "01", "11"),
                              stringsAsFactors = FALSE)
  # expand.grid ordering is history-fast; assign by named cells.
  observed$prob <- 0
  observed$prob[observed$history == "no" & observed$path == "00"] <- .25
  observed$prob[observed$history == "no" & observed$path == "01"] <- .25
  observed$prob[observed$history == "yes" & observed$path == "01"] <- .5
  intervention$prob <- 0
  intervention$prob[intervention$regime == "no" & intervention$path == "00"] <- .5
  intervention$prob[intervention$regime == "no" & intervention$path == "01"] <- .5
  intervention$prob[intervention$regime == "yes" & intervention$path == "01"] <- 1
  bad <- intervention
  bad$prob[bad$regime == "yes" & bad$path == "01"] <- 0
  bad$prob[bad$regime == "yes" & bad$path == "00"] <- 1
  expect_error(poc_longitudinal_model(observed, regimes, bad, absorbing = TRUE),
               "No finite")
  model <- poc_longitudinal_model(observed, regimes, intervention,
                                  absorbing = TRUE)
  expect_equal(model$n_response_types, 9)
  result <- poc_longitudinal_bounds(model,
    poc_longitudinal_query("pns", "yes", "no", 2))
  expect_lte(result$lower, .5 + 1e-9)
  expect_gte(result$upper, .5 - 1e-9)
  expect_equal(result$lower, .5, tolerance = 1e-9)
  expect_equal(result$upper, .5, tolerance = 1e-9)
  persistent <- poc_longitudinal_bounds(model,
    poc_longitudinal_query("persistent", "yes", "no", 2))
  expect_lte(persistent$upper, result$upper)
})

test_that("outcome-definition path restrictions apply to every potential outcome", {
  regimes <- rbind(no = c(0, 0), yes = c(1, 1))
  observed <- expand.grid(history = c("no", "yes"),
                          path = c("00", "01"), stringsAsFactors = FALSE)
  observed$prob <- rep(.25, 4)
  model <- poc_longitudinal_model(observed, regimes, absorbing = TRUE,
                                  allowed_paths = c("00", "01"))
  expect_equal(model$n_response_types, 4)
  first <- poc_longitudinal_bounds(model,
    poc_longitudinal_query("pns", "yes", "no", 1))
  expect_equal(c(first$lower, first$upper), c(0, 0))
  expect_error(poc_longitudinal_model(observed, regimes, absorbing = TRUE,
                                      allowed_paths = c("00", "02")),
               "allowed_paths")
})

test_that("small potential-path populations are enclosed by sharp bounds", {
  for (absorbing in c(FALSE, TRUE)) {
    labels <- if (absorbing) c("00", "01", "11") else
      c("00", "10", "01", "11")
    regimes <- rbind(no = c(0, 0), yes = c(1, 1))
    for (p0 in labels) for (p1 in labels) {
      observed <- expand.grid(history = c("no", "yes"), path = labels,
                              stringsAsFactors = FALSE)
      observed$prob <- 0
      observed$prob[observed$history == "no" & observed$path == p0] <- .5
      observed$prob[observed$history == "yes" & observed$path == p1] <- .5
      intervention <- expand.grid(regime = c("no", "yes"), path = labels,
                                  stringsAsFactors = FALSE)
      intervention$prob <- as.numeric(
        (intervention$regime == "no" & intervention$path == p0) |
        (intervention$regime == "yes" & intervention$path == p1))
      identified <- poc_longitudinal_model(observed, regimes, intervention,
                                            absorbing = absorbing)
      unknown <- poc_longitudinal_model(observed, regimes,
                                       absorbing = absorbing)
      for (t in 1:2) {
        truth <- as.numeric(substr(p1, t, t) == "1" &
                            substr(p0, t, t) == "0")
        q <- poc_longitudinal_query("pns", "yes", "no", t)
        exact <- poc_longitudinal_bounds(identified, q)
        wide <- poc_longitudinal_bounds(unknown, q)
        expect_equal(c(exact$lower, exact$upper), c(truth, truth),
                     tolerance = 1e-9)
        expect_lte(wide$lower, truth + 1e-9)
        expect_gte(wide$upper, truth - 1e-9)
      }
      extra <- list(
        list(kind = "lagged", truth = as.numeric(substr(p1, 2, 2) == "1" &
                                                   substr(p0, 1, 1) == "0")),
        list(kind = "event_time", truth = as.numeric(substr(p1, 1, 1) == "0" &
                                                      substr(p1, 2, 2) == "1" &
                                                      substr(p0, 2, 2) == "0")),
        list(kind = "persistent", truth = as.numeric(p1 == "11" & p0 == "00")))
      for (item in extra) {
        b <- poc_longitudinal_bounds(identified,
          poc_longitudinal_query(item$kind, "yes", "no", 2))
        expect_equal(c(b$lower, b$upper), rep(item$truth, 2),
                     tolerance = 1e-9)
      }
      pair <- poc_longitudinal_bounds(identified,
        poc_longitudinal_query("trajectory", "yes", "no", 2,
                               regime_path = p1, reference_path = p0))
      expect_equal(c(pair$lower, pair$upper), c(1, 1), tolerance = 1e-9)
    }
  }
})

test_that("undefined denominators and plain summaries remain explicit", {
  regimes <- rbind(no = 0, yes = 1)
  observed <- data.frame(history = c("no", "no", "yes", "yes"),
                         path = c("0", "1", "0", "1"),
                         prob = c(.5, 0, .5, 0))
  model <- poc_longitudinal_model(observed, regimes)
  result <- poc_longitudinal_bounds(model,
    poc_longitudinal_query("pn", "yes", "no", 1))
  expect_identical(result$status, "undefined_condition")
  expect_true(is.na(result$lower))
  expect_equal(result$denominator, 0)
  analysis <- poc_analysis(result, provenance = list(source = "fixture", seed = 17))
  summed <- summary(analysis)
  expect_identical(summed$bounds$status, "undefined_condition")
  log1 <- capture.output(print(summed))
  expect_identical(log1, capture.output(print(summed)))
  expect_true(any(grepl("undefined_condition", log1, fixed = TRUE)))
  expect_true(any(grepl("fixture", log1, fixed = TRUE)))
  directory <- tempfile("pnsbib-export-")
  paths <- poc_export_summary(analysis, directory, "tsv")
  expect_true(all(file.exists(paths)))
  from_disk <- utils::read.delim(file.path(directory, "bounds.tsv"),
                                 check.names = FALSE)
  expect_identical(from_disk$status, "undefined_condition")
  expect_true(is.na(from_disk$lower))
  expect_identical(readLines(file.path(directory, "analysis.log")), log1)
  static <- poc_model(matrix(c(.3, .1, .2, .4), 2, 2,
    dimnames = list(c("0", "1"), c("0", "1"))),
    matrix(c(.7, .3, .4, .6), 2, 2, byrow = TRUE,
      dimnames = list(c("0", "1"), c("0", "1"))))
  single <- poc_exact(static, poc_query(c("0" = "0", "1" = "1")))
  expect_true(inherits(summary(single), "summary.pnsbib_analysis"))
  expect_true(inherits(summary(result), "summary.pnsbib_analysis"))
  mixed <- poc_analysis(list(static = single, longitudinal = result),
                        provenance = list(source = "mixed_fixture"))
  expect_equal(nrow(summary(mixed)$bounds), 2)
  expect_true(all(c("observed_x", "regime", "horizon") %in%
                  names(summary(mixed)$bounds)))
  expect_true(any(grepl("mixed_fixture", capture.output(print(summary(mixed))),
                        fixed = TRUE)))
})
