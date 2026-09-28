# Independent eight-cell population, not constructed with package internals.
nested_population <- function(q) {
  response <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))[rep(1:4, 2), ]
  assignment <- rep(0:1, each = 4)
  factual <- response[cbind(seq_len(8), assignment + 1L)]
  observed <- intervention <- matrix(0, 2, 2,
    dimnames = list(c("0", "1"), c("0", "1")))
  for (x in 0:1) for (y in 0:1) {
    observed[x + 1L, y + 1L] <- sum(q[assignment == x & factual == y])
    intervention[x + 1L, y + 1L] <- sum(q[response[, x + 1L] == y])
  }
  benefit <- response[, 1L] == 0 & response[, 2L] == 1
  list(o = observed, a = intervention, q = q, r = response, x = assignment,
    y = factual, truth = c(PNS = sum(q[benefit]),
      PN = sum(q[benefit & assignment == 1]) / observed[2, 2],
      PS = sum(q[benefit & assignment == 0]) / observed[1, 1]))
}
nested_queries <- function() list(
  PNS = poc_query(c(`0` = "0", `1` = "1")),
  PN = poc_query(c(`0` = "0"), "1", "1", conditional = TRUE),
  PS = poc_query(c(`1` = "1"), "0", "0", conditional = TRUE))
nested_endpoints <- function(fit) c(fit$lower, fit$upper)

test_that("the known DGP and domain-neutral labels retain their meaning", {
  p <- nested_population(.5 * c(.5, .2, .2, .1, .3, .4, .1, .2))
  expect_equal(unname(p$truth), c(.3, 2/3, 2/7))
  expect_equal(unname(p$o), matrix(c(.35, .15, .20, .30), 2, byrow = TRUE))
  expect_equal(unname(p$a[, 2]), c(.3, .45))
  queries <- nested_queries()
  for (backend in c("reference", "zig")) {
    m <- poc_model(p$o, p$a, lp_backend = backend)
    for (fit in list(poc_bounds(m, queries$PNS), poc_shu2026(m, queries$PNS),
                     poc_exact(m, queries$PNS)))
      expect_equal(nested_endpoints(fit), c(.15, .45), tolerance = 1e-8)
    o <- p$o; a <- p$a
    dimnames(o) <- dimnames(a) <- list(c("reference", "active"), c("failure", "success"))
    q <- poc_query(c(reference = "failure", active = "success"))
    for (order in list(1:2, 2:1)) {
      labelled <- poc_model(o[order, rev(order)], a[order, rev(order)], lp_backend = backend)
      expect_equal(nested_endpoints(poc_exact(labelled, q)), c(.15, .45), tolerance = 1e-8)
    }
  }
})

test_that("nested oracle-centered risk bands retain truth and narrow all queries", {
  # Deterministic positive populations avoid modifying the caller's RNG state.
  queries <- nested_queries()
  for (index in 1:16) {
    w <- ((seq_len(8) * (2 * index + 1) + index^2) %% 29) + 1
    p <- nested_population(w / sum(w))
    previous <- NULL
    for (width in c(1, .5, .25, .1, 0)) {
      constraints <- lapply(0:1, function(x) {
        group <- p$x == 1 - x
        risk <- sum(p$q[group & p$r[, x + 1L] == 1]) / sum(p$q[group])
        poc_constraint(poc_query(setNames("1", as.character(x)),
          observed_x = as.character(1 - x), conditional = TRUE),
          max(0, risk - width), min(1, risk + width))
      })
      m <- poc_partial_model(p$o, constraints, lp_backend = "zig")
      current <- vapply(queries, function(q) nested_endpoints(poc_partial_bounds(m, q)), c(0, 0))
      expect_true(all(current[1, ] <= p$truth + 1e-8 & current[2, ] >= p$truth - 1e-8))
      for (name in names(queries))
        expect_equal(unname(current[, name]), nested_endpoints(poc_partial_bounds(m,
          queries[[name]], lp_backend = "reference")), tolerance = 1e-8)
      if (!is.null(previous)) {
        expect_true(all(current[1, ] >= previous[1, ] - 1e-8))
        expect_true(all(current[2, ] <= previous[2, ] + 1e-8))
      }
      previous <- current
    }
  }
})

test_that("joint equal-risk information is not its marginal box or its corners", {
  A <- rbind(rep(1, 4), c(0, 1, -1, 0))
  objective <- c(0, 1, 0, 0)
  solve_native <- function(A, b, maximize) {
    answer <- poc_lp_zig(nrow(A), ncol(A), as.numeric(A), b,
      rep(0L, nrow(A)), objective, maximize, 200000L)
    expect_equal(answer[1], 0)
    answer[2]
  }
  for (maximize in c(FALSE, TRUE)) {
    reference <- lpSolve::lp(if (maximize) "max" else "min", objective, A,
                             rep("=", nrow(A)), c(1, 0))
    expect_equal(reference$status, 0)
    expect_equal(solve_native(A, c(1, 0), maximize), reference$objval, tolerance = 1e-8)
    expect_equal(reference$objval, if (maximize) .5 else 0, tolerance = 1e-8)
    expect_equal(solve_native(matrix(1, 1, 4), 1, maximize),
                 if (maximize) 1 else 0, tolerance = 1e-8)
  }
  witness <- c(0, .5, .5, 0)
  expect_equal(as.vector(A %*% witness), c(1, 0))
  expect_equal(sum(witness * objective), .5)
  # Equal-risk parameter corners both force PNS=0; interior risk .5 attains .5.
  expect_equal(pmin(c(0, 1), 1 - c(0, 1)), c(0, 0))
})

test_that("post-outcome selection changes truth and preserves undefined PS", {
  full <- nested_population(.5 * c(.5, .2, .2, .1, .3, .4, .1, .2))
  retained <- full$q * (full$y == 1)
  selected <- nested_population(retained / sum(retained))
  expect_equal(unname(selected$truth["PNS"]), 4/9)
  expect_equal(unname(full$truth["PNS"]), .3)
  for (backend in c("reference", "zig")) {
    m <- poc_partial_model(selected$o, lp_backend = backend)
    result <- poc_partial_bounds(m, nested_queries()$PS)
    expect_identical(result$status, "undefined_condition")
    expect_true(is.na(result$lower) && is.na(result$upper))
  }
})

test_that("the installed tutorial executes without changing the RNG", {
  before <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  path <- system.file("doc", "nested-identification.R", package = "pnsbib")
  expect_true(nzchar(path))
  expect_output(sys.source(path, envir = new.env(parent = globalenv())), "width")
  after <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  expect_identical(after, before)
})
