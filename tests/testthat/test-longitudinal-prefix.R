prefix_fixture <- function() {
  regimes <- rbind(off = c(0, 0), later = c(0, 1),
                   early = c(1, 0), always = c(1, 1))
  observed <- expand.grid(history = rownames(regimes), path = c("00", "11"),
                          stringsAsFactors = FALSE)
  observed$prob <- c(.1, .2, 0, 0, 0, 0, .3, .4)
  margins <- data.frame(regime = rownames(regimes), horizon = 1,
                         prob1 = c(.3, .3, .8, .8))
  poc_longitudinal_model(observed, regimes, time_margins = margins,
                         allowed_paths = c("00", "11"))
}

test_that("prefix PN and PS pool continuations before optimizing", {
  model <- prefix_fixture()
  for (kind in c("pn", "ps")) {
    query <- poc_longitudinal_query(kind, "early", "off", 1,
                                    condition_on = "prefix")
    result <- poc_longitudinal_bounds(model, query)
    expected <- if (kind == "pn") 4 / 7 else 1 / 3
    denominator <- if (kind == "pn") .7 else .3
    expect_identical(result$status, "ok")
    expect_equal(c(result$lower, result$upper), rep(expected, 2), tolerance = 1e-8)
    expect_equal(result$denominator, denominator)
    expect_identical(result$condition_on, "prefix")
    expect_identical(result$condition_horizon, 1L)
    expect_match(result$event, "A[1:1]", fixed = TRUE)
    history <- poc_longitudinal_bounds(model,
      poc_longitudinal_query(kind, "early", "off", 1))
    expect_equal(c(history$lower, history$upper), c(0, 1), tolerance = 1e-8)
    expect_equal(history$denominator, if (kind == "pn") .3 else .1)
    expect_identical(history$condition_on, "history")
    expect_identical(history$condition_horizon, 2L)
    # Selecting the other continuation leaves an early-horizon prefix query unchanged.
    other <- poc_longitudinal_bounds(model,
      poc_longitudinal_query(kind, "always", "later", 1, condition_on = "prefix"))
    expect_equal(other$lower, result$lower, tolerance = 1e-8)
    expect_equal(other$upper, result$upper, tolerance = 1e-8)
    expect_equal(other$denominator, result$denominator)
  }
})

test_that("terminal, singleton-prefix and legacy queries preserve history results", {
  model <- prefix_fixture()
  for (kind in c("pn", "ps")) {
    for (t in 1:2) {
      # The restricted observed support has one history per first-period state.
      observed <- data.frame(history = c("off", "off", "early", "early"),
                              path = c("00", "11", "00", "11"),
                              prob = c(.25, 0, 0, .75))
      single <- poc_longitudinal_model(observed, model$regimes,
                                       allowed_paths = c("00", "11"))
      for (candidate in if (t == 2L) list(model, single) else list(single)) {
        query <- poc_longitudinal_query(kind, "early", "off", t)
        history <- poc_longitudinal_bounds(candidate, query)
        legacy <- query
        legacy$condition_on <- NULL
        expect_identical(poc_longitudinal_bounds(candidate, legacy), history)
        query$condition_on <- "prefix"
        prefix <- poc_longitudinal_bounds(candidate, query)
        expect_equal(prefix[, c("lower", "upper", "denominator", "status")],
                     history[, c("lower", "upper", "denominator", "status")])
      }
    }
  }
})

test_that("prefix T=1 reduces to static conditional PoC", {
  regimes <- rbind(off = 0, on = 1)
  observed <- data.frame(history = c("off", "off", "on", "on"),
                          path = c("0", "1", "0", "1"), prob = c(.3, .2, .1, .4))
  intervention <- data.frame(regime = observed$history, path = observed$path,
                              prob = c(.7, .3, .4, .6))
  model <- poc_longitudinal_model(observed, regimes, intervention)
  o <- matrix(observed$prob, 2, 2, byrow = TRUE,
               dimnames = list(c("off", "on"), c("0", "1")))
  a <- matrix(intervention$prob, 2, 2, byrow = TRUE, dimnames = dimnames(o))
  static <- poc_model(o, a)
  for (kind in c("pn", "ps")) {
    long <- poc_longitudinal_bounds(model,
      poc_longitudinal_query(kind, "on", "off", 1, condition_on = "prefix"))
    q <- if (kind == "pn") poc_query(c(off = "0"), "on", "1", TRUE) else
      poc_query(c(on = "1"), "off", "0", TRUE)
    exact <- poc_exact(static, q)
    expect_equal(c(long$lower, long$upper, long$denominator),
                 c(exact$lower, exact$upper, exact$denominator), tolerance = 1e-8)
  }
})

test_that("prefix conditions use actual Y even without no anticipation", {
  regimes <- rbind(named = c(1, 0), other = c(1, 1), ref = c(0, 0))
  for (kind in c("pn", "ps")) {
    observed <- data.frame(history = "other", path = c("00", "10", "01", "11"),
                            prob = if (kind == "pn") c(0, 0, 0, 1) else c(1, 0, 0, 0))
    margins <- data.frame(regime = c("named", "ref"), horizon = 1,
                           prob1 = if (kind == "pn") c(0, 0) else c(1, 1))
    model <- poc_longitudinal_model(observed, regimes, time_margins = margins,
                                     no_anticipation = FALSE)
    active <- if (kind == "pn") "named" else "ref"
    reference <- if (kind == "pn") "ref" else "named"
    query <- poc_longitudinal_query(kind, active, reference, 1,
                                    condition_on = "prefix")
    result <- poc_longitudinal_bounds(model, query)
    expect_equal(c(result$lower, result$upper, result$denominator), c(1, 1, 1))
    expect_equal(poc_longitudinal_sensitivity(model, query)$lower, 1)
    expect_error(poc_longitudinal_bounds(model,
      poc_longitudinal_query(kind, active, reference, 1)), "observed factual history")
    expect_error(poc_longitudinal_model(observed, regimes, time_margins = margins),
                 "No finite")
    query$horizon <- 2L
    expect_error(poc_longitudinal_bounds(model, query), "observed factual history")
  }
})

test_that("prefix bands retain exact limits, floors, and reporting", {
  model <- prefix_fixture()
  rows <- list()
  for (kind in c("pn", "ps")) {
    query <- poc_longitudinal_query(kind, "early", "off", 1,
                                    condition_on = "prefix")
    exact <- poc_longitudinal_bounds(model, query)
    zero <- poc_longitudinal_sensitivity(model, query)
    expect_equal(zero[, c("lower", "upper", "denominator")],
                 exact[, c("lower", "upper", "denominator")], tolerance = 1e-8)
    narrow <- poc_longitudinal_sensitivity(model, query, observed_delta = .01,
                                            time_delta = .02)
    wide <- poc_longitudinal_sensitivity(model, query, observed_delta = .03,
                                          time_delta = .05)
    expect_lte(narrow$lower, exact$lower + 1e-8)
    expect_gte(narrow$upper, exact$upper - 1e-8)
    expect_lte(wide$lower, narrow$lower + 1e-8)
    expect_gte(wide$upper, narrow$upper - 1e-8)
    expect_true(is.na(wide$denominator))
    expect_lt(wide$condition_probability_min, wide$condition_probability_max)
    floor <- poc_longitudinal_sensitivity(model, query, observed_delta = .03,
      time_delta = .05, minimum_condition_probability = exact$denominator)
    expect_lte(wide$lower, floor$lower + 1e-8)
    expect_gte(wide$upper, floor$upper - 1e-8)
    expect_equal(floor$condition_probability_min, wide$condition_probability_min)
    expect_identical(poc_longitudinal_sensitivity(model, query,
      minimum_condition_probability = .9)$status, "undefined_condition")
    rows[[kind]] <- wide
  }
  analysis <- poc_analysis(rows, provenance = list(fixture = "prefix_scope"))
  output <- tempfile("prefix-export-")
  poc_export_summary(analysis, output)
  back <- utils::read.delim(file.path(output, "bounds.tsv"))
  expect_identical(back$condition_on, rep("prefix", 2))
  expect_identical(back$condition_horizon, rep(1L, 2))
  expect_true(all(grepl("A[1:1]", back$event, fixed = TRUE)))
})

test_that("shared-prefix no anticipation and zero masses remain explicit", {
  model <- prefix_fixture()
  for (kind in c("pn", "ps")) {
    same <- if (kind == "pn") c("early", "always") else c("later", "off")
    query <- poc_longitudinal_query(kind, same[1], same[2], 1,
                                    condition_on = "prefix")
    result <- poc_longitudinal_bounds(model, query)
    expect_equal(c(result$lower, result$upper), c(0, 0))
    relaxed <- poc_longitudinal_sensitivity(model, query, observed_delta = .1,
                                             time_delta = .1)
    expect_equal(c(relaxed$lower, relaxed$upper), c(0, 0))
    observed <- expand.grid(history = rownames(model$regimes), path = c("00", "11"),
                            stringsAsFactors = FALSE)
    observed$prob <- if (kind == "pn") c(rep(.25, 4), rep(0, 4)) else
      c(rep(0, 4), rep(.25, 4))
    empty <- poc_longitudinal_model(observed, model$regimes, allowed_paths = c("00", "11"))
    result <- poc_longitudinal_bounds(empty, query)
    expect_identical(result$status, "undefined_condition")
    expect_true(is.na(result$lower))
    expect_equal(result$denominator, 0)
    expect_identical(result$condition_on, "prefix")
    expect_identical(poc_longitudinal_sensitivity(empty, query)$status, "undefined_condition")
    expect_identical(poc_longitudinal_sensitivity(empty, query,
      observed_delta = .05)$status, "ok")
  }
  expect_error(poc_longitudinal_query("pns", "early", "off", 1,
    condition_on = "prefix"), "requires PN or PS")
  expect_error(poc_longitudinal_query("pn", "early", "off", 1,
    condition_on = "typo"), "arg")
  ordinary <- poc_longitudinal_bounds(model,
    poc_longitudinal_query("pns", "early", "off", 1))
  expect_identical(ordinary$condition_on, "none")
  expect_true(is.na(ordinary$condition_horizon))
})

test_that("mixed synthetic prefix truth lies within exact and band bounds", {
  set.seed(20260927)
  regimes <- rbind(off = c("00", "00"), later = c("00", "11"),
                   early = c("11", "00"), always = c("11", "11"))
  n <- 80L
  for (absorbing in c(FALSE, TRUE)) {
    # Generate potential trajectories independently of the engine's enumeration.
    potentials <- array(0L, c(n, 4L, 2L))
    initial <- matrix(sample(0:1, n * 2, replace = TRUE), n, 2)
    for (g in 1:4) {
      potentials[, g, 1] <- initial[, if (g <= 2) 1 else 2]
      terminal <- sample(0:1, n, replace = TRUE)
      potentials[, g, 2] <- if (absorbing) pmax(potentials[, g, 1], terminal) else terminal
    }
    factual <- rep(1:4, length.out = n)
    weights <- runif(n)
    weights <- weights / sum(weights)
    labels <- if (absorbing) c("00", "01", "11") else c("00", "10", "01", "11")
    path_labels <- sapply(1:4, function(g)
      paste0(potentials[, g, 1], potentials[, g, 2]))
    actual_labels <- path_labels[cbind(seq_len(n), factual)]
    observed <- expand.grid(history = rownames(regimes), path = labels,
                            stringsAsFactors = FALSE)
    observed$prob <- vapply(seq_len(nrow(observed)), function(i)
      sum(weights[factual == match(observed$history[i], rownames(regimes)) &
                    actual_labels == observed$path[i]]), 0)
    margins <- expand.grid(regime = rownames(regimes), path = labels,
                           stringsAsFactors = FALSE)
    margins$prob <- vapply(seq_len(nrow(margins)), function(i)
      sum(weights[path_labels[, match(margins$regime[i], rownames(regimes))] ==
                    margins$path[i]]), 0)
    model <- poc_longitudinal_model(observed, regimes, margins, absorbing = absorbing)
    for (t in 1:2) for (kind in c("pn", "ps")) {
      query <- poc_longitudinal_query(kind, "early", "off", t, condition_on = "prefix")
      h <- if (kind == "pn") 3L else 1L
      matching <- vapply(factual, function(i)
        all(regimes[i, seq_len(t)] == regimes[h, seq_len(t)]), TRUE)
      actual_y <- potentials[cbind(seq_len(n), factual, rep(t, n))]
      condition <- matching & actual_y == if (kind == "pn") 1L else 0L
      event <- condition & if (kind == "pn") potentials[, 1, t] == 0L else
        potentials[, 3, t] == 1L
      denominator <- sum(weights[condition])
      truth <- sum(weights[event]) / denominator
      exact <- poc_longitudinal_bounds(model, query)
      band <- poc_longitudinal_sensitivity(model, query, observed_delta = .01,
                                            path_delta = .02)
      expect_identical(exact$status, "ok")
      expect_identical(band$status, "ok")
      expect_equal(exact$denominator, denominator, tolerance = 1e-8)
      expect_lte(exact$lower, truth + 1e-8)
      expect_gte(exact$upper, truth - 1e-8)
      expect_lte(band$lower, exact$lower + 1e-8)
      expect_gte(band$upper, exact$upper - 1e-8)
    }
  }
})
