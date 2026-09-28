test_that("installed adjusted-risk example connects all three information models", {
  env <- new.env(parent=globalenv())
  sys.source(system.file("doc","adjusted-risks.R",package="pnsbib"),envir=env)
  z <- env$pnsbib_adjusted_examples
  expect_identical(z$annual$backend,"zig")
  expect_equal(dim(z$annual$risks),c(2L,2L))
  expect_true(all(z$curve$status=="ok"))
  expect_lte(z$partial$lower,z$complete$lower+1e-8)
  expect_gte(z$partial$upper,z$complete$upper-1e-8)
  expect_true(all(z$cumulative$risks[,2]>=z$cumulative$risks[,1]))
  exports <- getNamespaceExports("pnsbib")
  expect_true(all(c("poc_risk_fit","poc_adjusted_risks") %in% exports))
  expect_false(any(c("poc_risk_fit_zig","poc_risk_standardize_zig") %in% exports))
})
