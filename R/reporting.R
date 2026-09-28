#' Collect probability-of-causation results for plain-text reporting
#'
#' The object carries already-computed bounds and explicit provenance. It never
#' converts observational risks into intervention margins.
#'
#' @param bounds A data frame of result rows, or a list of one-row result data frames.
#' @param cohort,observed,intervention,sensitivity,diagnostics Rectangular data
#'   frames; empty tables are allowed.
#' @param provenance Named list or two-column data frame of metadata.
#' @return A pnsbib_analysis object.
#' @export
poc_analysis <- function(bounds, cohort = data.frame(), observed = data.frame(),
                         intervention = data.frame(), sensitivity = data.frame(),
                         diagnostics = data.frame(), provenance = list()) {
  if (is.list(bounds) && !is.data.frame(bounds)) {
    if (!length(bounds) || !all(vapply(bounds, is.data.frame, logical(1L))))
      stop("bounds list must contain result data frames.", call. = FALSE)
    columns <- unique(unlist(lapply(bounds, names), use.names = FALSE))
    bounds <- do.call(rbind, lapply(bounds, function(z) {
      z <- as.data.frame(z)
      for (name in setdiff(columns, names(z))) z[[name]] <- NA
      z[, columns, drop = FALSE]
    }))
    rownames(bounds) <- NULL
  }
  if (!is.data.frame(bounds) ||
      !all(c("lower", "upper", "method", "status", "denominator") %in%
           names(bounds)))
    stop("bounds needs lower, upper, method, status, denominator columns.",
         call. = FALSE)
  tables <- list(cohort = cohort, observed = observed,
                 intervention = intervention, bounds = bounds,
                 sensitivity = sensitivity, diagnostics = diagnostics)
  if (!all(vapply(tables, is.data.frame, logical(1L))))
    stop("Every analysis table must be a data frame.", call. = FALSE)
  if (is.list(provenance) && !is.data.frame(provenance)) {
    if (is.null(names(provenance)) || any(!nzchar(names(provenance))) ||
        anyDuplicated(names(provenance)))
      stop("provenance must be a named list.", call. = FALSE)
    provenance <- data.frame(key = names(provenance),
      value = vapply(provenance, function(z) paste(as.character(z), collapse = ";"), ""),
      stringsAsFactors = FALSE)
  }
  if (!is.data.frame(provenance) ||
      !all(c("key", "value") %in% names(provenance)))
    stop("provenance needs key and value columns.", call. = FALSE)
  structure(c(tables, list(provenance = provenance)), class = "pnsbib_analysis")
}

.poc_plain_table <- function(x, digits = 6L) {
  if (!nrow(x)) return("(no rows)")
  y <- as.data.frame(x, stringsAsFactors = FALSE)
  for (j in seq_along(y)) if (is.numeric(y[[j]])) {
    y[[j]] <- ifelse(is.na(y[[j]]), "NA",
                     formatC(y[[j]], digits = digits, format = "f"))
  }
  utils::capture.output(print(y, row.names = FALSE, right = FALSE, na.print = "NA"))
}

#' @export
summary.pnsbib_result <- function(object, ...) {
  summary(poc_analysis(bounds = as.data.frame(object), provenance = list(
    result_type = "static_probability_bound")))
}

#' @export
summary.pnsbib_longitudinal_result <- function(object, ...) {
  summary(poc_analysis(bounds = as.data.frame(object), provenance = list(
    result_type = "longitudinal_probability_bound")))
}

#' @export
summary.pnsbib_analysis <- function(object, ...) {
  structure(unclass(object), class = "summary.pnsbib_analysis")
}

#' @export
print.summary.pnsbib_analysis <- function(x, ..., digits = 6L) {
  for (name in names(x)) {
    cat("[", name, "]\n", sep = "")
    cat(paste(.poc_plain_table(x[[name]], digits = digits), collapse = "\n"), "\n",
        sep = "")
  }
  invisible(x)
}

#' @export
print.pnsbib_analysis <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Export all analysis summary tables and a human-readable log
#'
#' @param x A pnsbib_analysis or its summary.
#' @param directory Existing or new output directory; files must not exist.
#' @param format `tsv` or `csv`.
#' @return Invisibly, the written paths.
#' @export
poc_export_summary <- function(x, directory, format = c("tsv", "csv")) {
  format <- match.arg(format)
  if (inherits(x, "pnsbib_analysis")) x <- summary(x)
  if (!inherits(x, "summary.pnsbib_analysis") ||
      !is.character(directory) || length(directory) != 1L ||
      is.na(directory) || !nzchar(directory))
    stop("Expected an analysis summary and output directory.", call. = FALSE)
  paths <- file.path(directory, paste0(names(x), ".", format))
  log_path <- file.path(directory, "analysis.log")
  if (any(file.exists(c(paths, log_path))))
    stop("Output files already exist; use a new directory.", call. = FALSE)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(x)) {
    utils::write.table(x[[i]], paths[i], sep = if (format == "tsv") "\t" else ",",
                       row.names = FALSE, quote = TRUE, na = "NA")
  }
  utils::capture.output(print(x), file = log_path)
  invisible(c(paths, log_path))
}
