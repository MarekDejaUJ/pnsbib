plot_rows <- function() {
  data.frame(estimand = c("PNS", "PNS", "PN"), lower = c(.1, .15, NA_real_),
    upper = c(.8, .7, NA_real_), method = c("8", "exact_lp", "exact_lp"),
    sharpness = c("valid_recursive", "sharp", "sharp"),
    status = c("ok", "ok", "undefined_condition"), denominator = c(1, 1, 0),
    input_status = "synthetic_known_margins", stringsAsFactors = FALSE)
}

test_that("plots preserve raw results, undefined rows and provenance", {
  rows <- plot_rows()
  original <- serialize(rows, NULL)
  analysis <- poc_analysis(rows, provenance = list(population = "synthetic"))
  p <- poc_plot(analysis, draw = FALSE)
  expect_s3_class(p, "pnsbib_plot")
  expect_identical(p$data, rows)
  expect_identical(p$layout$source_row, 1:3)
  expect_identical(p$layout$defined, c(TRUE, TRUE, FALSE))
  expect_identical(p$layout$style, c(1L, 2L, 2L))
  expect_identical(p$provenance, analysis$provenance)
  expect_identical(serialize(rows, NULL), original)
  expect_match(p$layout$label[1], "Li-Pearl 8")
  expect_match(p$layout$label[2], "Exact LP")
  expect_match(p$subtitle, "synthetic_known_margins")
  expect_identical(poc_plot(rows, labels = c("A", "B", "C"), draw = FALSE)$layout$label,
                   c("A", "B", "C"))
})

test_that("analytical and LP outputs feed plots without losing endpoints", {
  o <- matrix(c(.2, .1, .3, .4), 2, byrow = TRUE,
              dimnames = list(c("off", "on"), c("0", "1")))
  a <- matrix(c(.6, .4, .3, .7), 2, byrow = TRUE, dimnames = dimnames(o))
  m <- poc_model(o, a); q <- poc_query(c(off = "0", on = "1"))
  fits <- list(poc_bounds(m, q), poc_exact(m, q))
  p <- poc_plot(poc_analysis(fits, provenance = list(inputs = "synthetic")), draw = FALSE)
  expect_equal(p$data$lower, vapply(fits, function(z) z$lower, 0))
  expect_equal(p$data$upper, vapply(fits, function(z) z$upper, 0))
  expect_equal(p$data$method, vapply(fits, function(z) z$method, ""))
})

test_that("curve grouping retains order, decreasing bounds and undefined gaps", {
  rows <- plot_rows()[rep(1, 3), ]
  rownames(rows) <- NULL
  rows$horizon <- c(3, 1, 2)
  rows$lower <- rows$upper <- c(0, 1, NA_real_)
  rows$status <- c("ok", "ok", "undefined_condition")
  p <- poc_plot(rows, "curve", draw = FALSE)
  expect_identical(p$data, rows)
  expect_equal(p$layout$horizon, c(3, 1, 2))
  expect_identical(p$layout$defined, c(TRUE, TRUE, FALSE))
  expect_equal(p$layout$group, c(1L, 1L, 1L))
  duplicate <- rbind(rows, rows[1, ])
  expect_error(poc_plot(duplicate, "curve", draw = FALSE), "Duplicate horizons")
  mixed <- rows; mixed$method[1] <- "exact_lp"
  expect_error(poc_plot(mixed, "curve", series = rep("one", 3), draw = FALSE), "mixes method")
  mixed <- rows; mixed$assumptions <- c("strong", "weak", "weak")
  auto <- poc_plot(mixed, "curve", draw = FALSE)
  expect_equal(length(unique(auto$layout$series)), 2)
  expect_error(poc_plot(mixed, "curve", series = rep("one", 3), draw = FALSE), "mixes assumptions")
  all_undefined <- rows; all_undefined$lower <- all_undefined$upper <- NA_real_
  all_undefined$status <- "undefined_condition"
  expect_false(any(poc_plot(all_undefined, "curve", draw = FALSE)$layout$defined))
})

test_that("plot input validation never repairs invalid or missing bounds", {
  rows <- plot_rows()
  for (x in list(NULL, list(), data.frame(), rows[, -1:-4]))
    expect_error(poc_plot(x, draw = FALSE), "result rows")
  for (value in c(-.01, 1.1, Inf, NA_real_)) {
    bad <- rows; bad$lower[1] <- value
    expect_error(poc_plot(bad, draw = FALSE), "Bounds")
  }
  bad <- rows; bad$lower <- as.character(bad$lower)
  expect_error(poc_plot(bad, draw = FALSE), "numeric")
  bad <- rows; bad$method[1] <- NA_character_
  expect_error(poc_plot(bad, draw = FALSE), "method")
  expect_error(poc_plot(rows, labels = "one", draw = FALSE), "labels")
  expect_error(poc_plot(rows, title = NA_character_, draw = FALSE), "title")
  expect_error(poc_plot(rows, monochrome = NA, draw = FALSE), "logical")
  expect_error(poc_plot(rows, "curve", draw = FALSE), "horizon")
})

test_that("export preserves devices, data and files in PDF and PNG", {
  directory <- tempfile("pnsbib-plot-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  p <- poc_plot(plot_rows(), draw = FALSE)
  caller <- file.path(directory, "caller.pdf")
  grDevices::pdf(caller)
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device), add = TRUE)
  old <- graphics::par("mar")
  expect_identical(plot(p), p)
  expect_identical(graphics::par("mar"), old)
  formats <- if (capabilities("png")) c("png", "pdf") else "pdf"
  files <- poc_save_plot(p, file.path(directory, "bounds"), formats = formats, dpi = 100)
  expect_identical(grDevices::dev.cur(), device)
  expect_true(all(file.exists(files)))
  expect_true(all(file.info(files)$size > 0))
  expect_identical(readRDS(files[["plot.rds"]]), p)
  restored <- read.delim(files[["data.tsv"]], stringsAsFactors = FALSE)
  expect_equal(restored, p$data)
  expect_identical(rawToChar(readBin(files[["pdf"]], "raw", 5)), "%PDF-")
  if ("png" %in% formats)
    expect_identical(as.integer(readBin(files[["png"]], "raw", 8)), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))
  hashes <- tools::md5sum(files)
  expect_error(poc_save_plot(p, file.path(directory, "bounds"), formats = formats), "already exist")
  expect_identical(tools::md5sum(files), hashes)
  # Even a sidecar collision blocks the complete export before creating images.
  writeLines("sentinel", file.path(directory, "collision.data.tsv"))
  expect_error(poc_save_plot(p, file.path(directory, "collision")), "already exist")
  expect_false(file.exists(file.path(directory, "collision.pdf")))
  expect_identical(readLines(file.path(directory, "collision.data.tsv")), "sentinel")
  for (width in list(NA_real_, -1, 4, Inf))
    expect_error(poc_save_plot(p, file.path(directory, "bad"), width = width))
  expect_error(poc_save_plot(p, file.path(directory, "bad.png")), "prefix")
  expect_error(poc_save_plot(p, file.path(directory, "bad%03d")), "prefix")
  expect_error(poc_save_plot(p, file.path(directory, "bad"), formats = "svg"), "formats")
  expect_error(poc_save_plot(p, file.path(directory, "bad"), formats = c("png", "png")), "formats")
  bad <- p; bad$data$lower <- "invalid"
  expect_error(poc_save_plot(bad, file.path(directory, "failure"), formats = "pdf"))
  expect_identical(grDevices::dev.cur(), device)
})

test_that("curve graphics render all-undefined and single-horizon inputs", {
  directory <- tempfile("pnsbib-curve-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  rows <- plot_rows()[3, ]; rows$horizon <- 2
  for (mono in c(FALSE, TRUE)) {
    p <- poc_plot(rows, "curve", monochrome = mono, draw = FALSE)
    files <- poc_save_plot(p, file.path(directory, paste0("curve", mono)), formats = "pdf")
    expect_true(file.info(files[["pdf"]])$size > 100)
    expect_identical(readRDS(files[["plot.rds"]])$data, rows)
  }
})

test_that("a device startup warning cannot close or draw on the caller device", {
  directory <- tempfile("pnsbib-device-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  grDevices::pdf(file.path(directory, "caller.pdf"))
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device), add = TRUE)
  testthat::local_mocked_bindings(pdf = function(...) warning("synthetic device failure"),
                                .package = "grDevices")
  p <- poc_plot(plot_rows(), draw = FALSE)
  expect_error(poc_save_plot(p, file.path(directory, "failed"), formats = "pdf"),
               "synthetic device failure")
  expect_identical(grDevices::dev.cur(), device)
  expect_false(file.exists(file.path(directory, "failed.data.tsv")))
})
