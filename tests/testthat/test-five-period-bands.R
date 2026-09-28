five_period_fixture <- function(absorbing = FALSE) {
  regimes <- rbind(reference = rep(0, 5), active = rep(1, 5))
  paths <- if (absorbing) t(vapply(0:5, function(k)
    c(rep(0, k), rep(1, 5-k)), numeric(5))) else
    as.matrix(expand.grid(rep(list(0:1), 5)))
  labels <- apply(paths, 1, paste0, collapse = "")
  observed <- expand.grid(history = rownames(regimes), path = labels,
                          stringsAsFactors = FALSE)
  probabilities <- if (absorbing)
    list(reference = c(.1, .1, .1, .1, .1, .5),
         active = c(.2, .2, .2, .1, .1, .2)) else
    lapply(c(reference = .4, active = .6), function(p)
      apply(paths, 1, function(y) prod(ifelse(y == 1, p, 1-p))))
  observed$prob <- vapply(seq_len(nrow(observed)), function(i)
    .5 * probabilities[[observed$history[i]]][match(observed$path[i], labels)], numeric(1))
  margins <- expand.grid(regime = rownames(regimes), horizon = 1:5,
                         stringsAsFactors = FALSE)
  margins$prob1 <- vapply(seq_len(nrow(margins)), function(i)
    sum(probabilities[[margins$regime[i]]] * paths[, margins$horizon[i]]), numeric(1))
  poc_longitudinal_model(observed, regimes, time_margins = margins,
                        absorbing = absorbing, lp_backend = "zig")
}

test_that("five-period annual sufficiency bands retain feasible endpoints", {
  model <- five_period_fixture()
  query <- poc_longitudinal_query("ps", "active", "reference", 4)
  result <- poc_longitudinal_sensitivity(model, query, time_delta = .05)
  reference <- poc_longitudinal_sensitivity(model, query, time_delta = .05,
                                           lp_backend = "reference")
  expect_identical(result$status, "ok")
  expect_equal(c(result$lower, result$upper), c(1/6, 1), tolerance = 1e-8)
  expect_equal(c(result$lower, result$upper),
               c(reference$lower, reference$upper), tolerance = 1e-8)
})

test_that("five-period annual and cumulative curves agree across LP backends", {
  for (absorbing in c(FALSE, TRUE)) {
    model <- five_period_fixture(absorbing)
    for (bands in list(NULL, list(time_delta = .05))) {
      native <- poc_curve(model, "active", "reference", bands = bands)
      reference <- poc_curve(model, "active", "reference", bands = bands,
                             lp_backend = "reference")
      expect_equal(nrow(native), 15L)
      expect_true(all(native$status == "ok"))
      expect_equal(native$lower, reference$lower, tolerance = 1e-8)
      expect_equal(native$upper, reference$upper, tolerance = 1e-8)
    }
  }
})
