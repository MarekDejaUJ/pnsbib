test_that("selected single-term fixtures agree with LP and Zig matches R", {
  expect_true("poc_single_zig" %in%
                names(getDLLRegisteredRoutines("pnsbib")[[".Call"]]))
  expect_equal(poc_single_zig(4L, 0.4, 0.2, 0.3, 0, 0.5, 0,
                              0, 0, numeric()), c(0.2, 0.3))
  f <- fixture_from_response_types()
  model <- poc_model(f$o, f$a)
  queries <- list(
    poc_query(c(x1 = "y1"), observed_y = "y1"),
    poc_query(c(x1 = "y1"), observed_y = "y2"),
    poc_query(c(x1 = "y1"), observed_x = "x2"),
    poc_query(c(x1 = "y1"), observed_x = "x2", observed_y = "y3")
  )
  for (q in queries) {
    paper_r <- poc_bounds(model, q, use_zig = FALSE)
    paper_zig <- poc_bounds(model, q, use_zig = TRUE)
    exact <- poc_exact(model, q)
    expect_equal(paper_zig$lower, paper_r$lower, tolerance = 1e-12)
    expect_equal(paper_zig$upper, paper_r$upper, tolerance = 1e-12)
    expect_equal(paper_r$lower, exact$lower, tolerance = 1e-8)
    expect_equal(paper_r$upper, exact$upper, tolerance = 1e-8)
    truth <- fixture_truth(f, q)
    expect_lte(paper_r$lower, truth + 1e-9)
    expect_gte(paper_r$upper, truth - 1e-9)
  }
})

test_that("all four multi-term theorem families contain LP bounds and truth", {
  for (seed in c(19L, 61L)) {
    f <- fixture_from_response_types(seed = seed)
    model <- poc_model(f$o, f$a)
    queries <- list(
      poc_query(c(x1 = "y1", x2 = "y3")),
      poc_query(c(x1 = "y1", x2 = "y3"), observed_x = "x3"),
      poc_query(c(x1 = "y1", x2 = "y3"), observed_y = "y2"),
      poc_query(c(x1 = "y1", x2 = "y3"), observed_x = "x3", observed_y = "y2"),
      poc_query(c(x1 = "y1", x2 = "y2", x3 = "y3"))
    )
    for (q in queries) {
      paper <- poc_bounds(model, q)
      exact <- poc_exact(model, q)
      truth <- fixture_truth(f, q)
      expect_lte(paper$lower, exact$lower + 1e-8)
      expect_gte(paper$upper, exact$upper - 1e-8)
      expect_lte(exact$lower, truth + 1e-8)
      expect_gte(exact$upper, truth - 1e-8)
    }
  }
})

test_that("Li-Pearl's three-treatment example reproduces its reported interval", {
  o <- rbind(c(238, 20, 7), c(10, 77, 259), c(147, 72, 70)) / 900
  a <- rbind(c(80, 7, 213), c(184, 29, 87), c(87, 189, 24)) / 300
  dimnames(o) <- dimnames(a) <- list(paste0("x", 1:3), paste0("y", 1:3))
  model <- poc_model(o, a)
  q <- poc_query(c(x1 = "y3", x2 = "y1", x3 = "y2"))
  paper <- poc_bounds(model, q)
  expect_equal(paper$lower, 0, tolerance = 1e-10)
  expect_equal(round(paper$upper, 3), 0.099)
})
