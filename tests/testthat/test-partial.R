partial_observed <- function() matrix(c(.3, .2, .1, .4), 2, 2, byrow = TRUE,
  dimnames = list(c("0", "1"), c("0", "1")))

partial_treated <- function(lower = .4, upper = lower, stratum = NULL) {
  poc_constraint(poc_query(c("0" = "1"), observed_x = "1", conditional = TRUE),
                  lower, upper, stratum = stratum, label = "treated_counterfactual")
}

partial_witness_value <- function(witness, query, stratum = NULL) {
  selected <- if (is.null(stratum)) rep(TRUE, nrow(witness)) else witness$stratum == stratum
  event <- selected
  factual <- selected
  for (x in names(query$counterfactual))
    event <- event & witness[[paste0("Y_", x)]] == query$counterfactual[[x]]
  if (!is.null(query$observed_x)) factual <- factual & witness$factual_x == query$observed_x
  if (!is.null(query$observed_y)) {
    actual <- vapply(seq_len(nrow(witness)), function(i)
      witness[[paste0("Y_", witness$factual_x[i])]][i], "")
    factual <- factual & actual == query$observed_y
  }
  denominator <- sum(witness$prob[if (query$conditional) factual else selected])
  sum(witness$prob[event & factual]) / denominator
}

test_that("partial constraints preserve treated and whole-population distinctions", {
  o <- partial_observed()
  pns <- poc_query(c("0" = "0", "1" = "1"))
  pn <- poc_query(c("0" = "0"), "1", "1", TRUE)
  ps <- poc_query(c("1" = "1"), "0", "0", TRUE)
  unknown <- poc_partial_model(o)
  limited <- poc_partial_model(o, list(partial_treated()))
  wrong_population <- poc_partial_model(o, list(poc_constraint(
    poc_query(c("0" = "1")), .4, label = "whole_population")))
  expect_equal(unlist(poc_partial_bounds(unknown, pns)[c("lower", "upper")]),
                 c(lower = 0, upper = .7), tolerance = 1e-9)
  for (item in list(list(q = pns, limits = c(.2, .6)),
                   list(q = pn, limits = c(.5, .75)),
                   list(q = ps, limits = c(0, 1)))) {
    result <- poc_partial_bounds(limited, item$q, return_witness = TRUE)
    expect_equal(c(result$lower, result$upper), item$limits, tolerance = 1e-8)
    expect_identical(result$status, "ok")
    for (end in c("lower", "upper")) {
      witness <- attr(result, "witnesses")[[end]]
      expect_equal(partial_witness_value(witness, item$q), result[[end]], tolerance = 1e-8)
      expect_equal(partial_witness_value(witness, partial_treated()$query), .4, tolerance = 1e-8)
    }
  }
  # Here both lead to the same global Y0 risk only because observed untreated
  # Y0 mass is .2 and the treated mass is .5. Change the subgroup risk to expose
  # the distinction, rather than concluding that its value is globally portable.
  subgroup <- poc_partial_model(o, list(partial_treated(.6)))
  marginal <- poc_partial_model(o, list(poc_constraint(poc_query(c("0" = "1")), .6)))
  a <- poc_partial_bounds(subgroup, poc_query(c("0" = "1")))
  b <- poc_partial_bounds(marginal, poc_query(c("0" = "1")))
  expect_equal(c(a$lower, a$upper), c(.5, .5), tolerance = 1e-8)
  expect_equal(c(b$lower, b$upper), c(.6, .6), tolerance = 1e-8)
  expect_equal(limited$constraint_table$denominator, .5)
  expect_equal(wrong_population$constraint_table$denominator, 1)
  wide <- poc_partial_model(o, list(partial_treated(.2, .6)))
  expect_lte(poc_partial_bounds(wide, pn)$lower, poc_partial_bounds(limited, pn)$lower)
  expect_gte(poc_partial_bounds(wide, pn)$upper, poc_partial_bounds(limited, pn)$upper)
})

test_that("strata use joint weights and combine pooled and local information", {
  o <- partial_observed()
  strata <- list(a = .25 * o, b = .75 * o)
  partial <- poc_partial_model(strata, list(partial_treated(stratum = "a")))
  query <- poc_query(c("0" = "0"), "1", "1", TRUE)
  global <- poc_partial_bounds(partial, query)
  local <- poc_partial_bounds(partial, query, stratum = "a")
  expect_equal(c(global$lower, global$upper, global$denominator), c(.125, .9375, .4), tolerance = 1e-8)
  expect_equal(c(local$lower, local$upper, local$denominator), c(.5, .75, .1), tolerance = 1e-8)
  expect_identical(local$stratum, "a")
  expect_equal(partial$constraint_table$denominator, .125)
  pooled <- poc_constraint(poc_query(c("0" = "1")), .4, label = "pooled_risk")
  combined <- poc_partial_model(strata, list(partial_treated(stratum = "a"), pooled))
  total <- poc_partial_bounds(combined, query, return_witness = TRUE)
  expect_equal(c(total$lower, total$upper), c(.5, .75), tolerance = 1e-8)
  for (w in attr(total, "witnesses")) {
    expect_equal(sum(w$prob[w$stratum == "a"]), .25, tolerance = 1e-8)
    expect_equal(sum(w$prob[w$stratum == "b"]), .75, tolerance = 1e-8)
    expect_equal(partial_witness_value(w, pooled$query), .4, tolerance = 1e-8)
    expect_equal(partial_witness_value(w, partial_treated()$query, "a"), .4, tolerance = 1e-8)
  }
  # Stratum conditioning applies even for a nonconditional query.
  risk_a <- poc_partial_bounds(partial, poc_query(c("0" = "1")), "a")
  expect_equal(c(risk_a$lower, risk_a$upper, risk_a$denominator), c(.4, .4, .25), tolerance = 1e-8)
  expect_error(poc_partial_model(list(a = o, b = o)), "sum to one across")
})

test_that("complete supplied margins reduce to the existing multivalued exact LP", {
  set.seed(20260928)
  for (n in 2:3) {
    m <- 3L
    treatments <- c("none", "one", "both")
    outcomes <- as.character(seq_len(n))
    N <- 180L
    potential <- matrix(sample(outcomes, N * m, replace = TRUE), N, m,
                          dimnames = list(NULL, treatments))
    factual <- sample(treatments, N, replace = TRUE)
    actual <- potential[cbind(seq_len(N), match(factual, treatments))]
    w <- runif(N)
    w <- w / sum(w)
    o <- a <- matrix(0, m, n, dimnames = list(treatments, outcomes))
    for (x in treatments) for (y in outcomes) {
      o[x, y] <- sum(w[factual == x & actual == y])
      a[x, y] <- sum(w[potential[, x] == y])
    }
    constraints <- list()
    for (x in treatments) for (y in outcomes)
      constraints[[length(constraints) + 1L]] <- poc_constraint(
        poc_query(stats::setNames(y, x)), a[x, y], label = paste(x, y))
    full <- poc_model(o, a)
    partial <- poc_partial_model(o, constraints)
    queries <- list(poc_query(c(none = "1")),
      poc_query(c(none = "1"), observed_y = "1"),
      poc_query(c(none = "1"), observed_y = "2"),
      poc_query(c(none = "1"), observed_x = "one"),
      poc_query(c(none = "1"), "one", "2"),
      poc_query(c(none = "1", one = "2")),
      poc_query(c(none = "1", one = "2"), observed_x = "both"),
      poc_query(c(none = "1", one = "2"), observed_y = "2"),
      poc_query(c(none = "1", one = "2"), "both", "2"),
      poc_query(c(none = "1", one = "2", both = "1")),
      poc_query(c(none = "1"), "none", "2", TRUE))
    for (q in queries) for (conditional in c(FALSE, TRUE)) {
      if (conditional && is.null(q$observed_x) && is.null(q$observed_y)) next
      q$conditional <- conditional
      reference <- poc_exact(full, q)
      result <- poc_partial_bounds(partial, q, return_witness = TRUE)
      expect_identical(result$status, "ok")
      expect_equal(c(result$lower, result$upper, result$denominator),
                   c(reference$lower, reference$upper, reference$denominator), tolerance = 1e-8)
      for (end in c("lower", "upper")) {
        witness <- attr(result, "witnesses")[[end]]
        expect_equal(sum(witness$prob), 1, tolerance = 1e-8)
        expect_gte(min(witness$prob), -1e-8)
        expect_equal(partial_witness_value(witness, q), result[[end]], tolerance = 1e-8)
        for (x in treatments) for (y in outcomes) {
          expect_equal(sum(witness$prob[witness$factual_x == x &
                             witness[[paste0("Y_", x)]] == y]), o[x, y], tolerance = 1e-8)
          expect_equal(sum(witness$prob[witness[[paste0("Y_", x)]] == y]), a[x, y], tolerance = 1e-8)
        }
      }
    }
    # Half the margins with nonzero widths must contain this independently
    # generated population's multi-hypothetical event probability.
    selected <- constraints[seq(1, length(constraints), by = 2)]
    selected <- lapply(selected, function(c) poc_constraint(c$query,
      max(0, c$lower - .02), min(1, c$upper + .02), label = c$label))
    limited <- poc_partial_model(o, selected)
    q <- poc_query(c(none = "1", one = "2", both = "1"))
    truth <- sum(w[potential[, "none"] == "1" & potential[, "one"] == "2" &
                     potential[, "both"] == "1"])
    limits <- poc_partial_bounds(limited, q)
    expect_lte(limits$lower, truth + 1e-8)
    expect_gte(limits$upper, truth - 1e-8)
  }
})

test_that("invalid constraints, zero populations and size limits are explicit", {
  o <- partial_observed()
  q <- poc_query(c("0" = "0"))
  for (bounds in list(c(-.1, .5), c(.3, 1.1), c(.8, .2), c(NA, .2), c(0, Inf)))
    expect_error(poc_constraint(q, bounds[1], bounds[2]), "ordered probability")
  expect_error(poc_constraint(q, numeric()), "ordered probability")
  expect_error(poc_constraint(q, .2, label = ""), "optional labels")
  expect_error(poc_constraint(q, .2, stratum = c("a", "b")), "optional labels")
  expect_error(poc_constraint("bad", .2), "query")
  expect_error(poc_partial_model(o, list(partial_treated(), partial_treated())), "unique")
  expect_error(poc_partial_model(o, list(q)), "poc_constraint")
  expect_error(poc_partial_model(o, list(partial_treated(stratum = "missing"))), "stratum absent")
  expect_error(poc_partial_model(o, list(poc_constraint(poc_query(c(bad = "0")), .2))), "absent")
  expect_error(poc_partial_model(o, list(poc_constraint(poc_query(c("0" = "1")), .1))), "No response-type")
  expect_error(poc_partial_model(o, max_variables = 7), "above max_variables")
  expect_error(poc_partial_model(list(a = o / 2, b = o / 2), max_variables = 15), "above max_variables")
  expect_error(poc_partial_model(o, tolerance = NA_real_), "Invalid tolerance")
  expect_error(poc_partial_model(o, max_variables = Inf), "max_variables")
  expect_error(poc_partial_model(unname(o)), "named probability")
  expect_error(poc_partial_model(list(o)), "uniquely named")
  expect_error(poc_partial_model(list(a = o / 2, b = o[2:1, ] / 2)), "matching labels")
  invalid <- o
  invalid[1, 1] <- NA_real_
  expect_error(poc_partial_model(invalid), "probability matrices")
  invalid[1, 1] <- -.1
  expect_error(poc_partial_model(invalid), "probability matrices")
  model <- poc_partial_model(list(empty = 0 * o, full = o))
  result <- poc_partial_bounds(model, q, "empty", return_witness = TRUE)
  expect_identical(result$status, "undefined_condition")
  expect_true(is.na(result$lower))
  expect_equal(result$denominator, 0)
  expect_null(attr(result, "witnesses"))
  expect_error(poc_partial_model(list(empty = 0 * o, full = o),
    list(poc_constraint(q, 0, 1, stratum = "empty"))), "zero-probability")
  one <- o
  one[1, ] <- 0
  one[2, ] <- one[2, ] / sum(one[2, ])
  factual <- poc_query(c("1" = "1"), observed_x = "0", conditional = TRUE)
  expect_identical(poc_partial_bounds(poc_partial_model(one), factual)$status, "undefined_condition")
  expect_error(poc_partial_model(one, list(poc_constraint(factual, 0, 1))), "zero-probability")
  expect_error(poc_partial_bounds(model, q, stratum = "unknown"), "stratum absent")
  expect_error(poc_partial_bounds(model, q, return_witness = NA), "logical")
  expect_error(poc_exact(model, q), "pnsbib_model")
  expect_error(poc_bounds(model, q), "pnsbib_model")
})

test_that("partial information labels survive analysis summaries and exports", {
  model <- poc_partial_model(partial_observed(), list(partial_treated()))
  result <- poc_partial_bounds(model, poc_query(c("0" = "0"), "1", "1", TRUE))
  analysis <- poc_analysis(result, intervention = model$constraint_table,
                            provenance = list(input_type = "hypothetical_subgroup_fixture"))
  output <- tempfile("partial-export-")
  poc_export_summary(analysis, output)
  bounds <- utils::read.delim(file.path(output, "bounds.tsv"))
  constraints <- utils::read.delim(file.path(output, "intervention.tsv"))
  expect_identical(bounds$method, "partial_exact_lp")
  expect_identical(bounds$sharpness, "sharp_given_constraints")
  expect_identical(bounds$constraint_labels, "treated_counterfactual")
  expect_equal(bounds$denominator, .4)
  expect_equal(constraints$denominator, .5)
  expect_true(constraints$conditional)
  expect_identical(constraints$observed_x, 1L)
  expect_equal(bounds$lower, .5, tolerance = 1e-8)
})
