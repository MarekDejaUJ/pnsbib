.poc_plot_text <- function(x, n, name) {
  if (!is.character(x) || length(x) != n || anyNA(x) || any(!nzchar(x)))
    stop(name, " must contain one nonempty character label per row.", call. = FALSE)
  x
}

.poc_plot_method <- function(x) {
  out <- as.character(x)
  analytical <- out %in% c(as.character(4:11), "margin", "consistency")
  out[analytical] <- paste("Li-Pearl", out[analytical])
  out[out == "exact_lp"] <- "Exact LP"
  out[out == "longitudinal_exact_lp"] <- "Longitudinal LP"
  out
}

#' Visualize already-computed probability bounds
#'
#' No model is fitted and no bound is recomputed. Source rows, assumptions,
#' statuses and provenance are retained in the returned plot object. Intervals
#' are probability bounds, with the interpretation of the supplied inputs;
#' they are not automatically sampling confidence intervals. Curve endpoints
#' are connected without smoothing or imposing monotonicity. Undefined rows
#' remain visible and break curves.
#'
#' @param x A result data frame or a [poc_analysis()] object.
#' @param type Horizontal interval comparison or bounds over numeric horizons.
#' @param labels Optional character vector of row labels for interval plots.
#' @param series Optional character vector assigning curve rows to series.
#'   Within a series, methods, estimands, assumptions and declared scenarios
#'   must agree. Default grouping uses these fields when present.
#' @param title,subtitle Optional single character strings. The default subtitle
#'   explicitly describes probability bounds rather than standard errors.
#' @param monochrome Use grayscale with redundant markers and line types.
#' @param draw Draw immediately; FALSE constructs a serializable plot object.
#' @return Invisibly, a pnsbib_plot with unchanged source data, layout and
#'   provenance. Use [poc_save_plot()] to export it.
#' @examples
#' o <- matrix(c(.2, .1, .3, .4), 2, byrow = TRUE,
#'             dimnames = list(c("off", "on"), c("0", "1")))
#' a <- matrix(c(.6, .4, .3, .7), 2, byrow = TRUE, dimnames = dimnames(o))
#' m <- poc_model(o, a)
#' q <- poc_query(c(off = "0", on = "1"))
#' results <- poc_analysis(list(poc_bounds(m, q), poc_exact(m, q)),
#'                         provenance = list(inputs = "synthetic example"))
#' p <- poc_plot(results, draw = FALSE)
#' plot(p)
#' @export
poc_plot <- function(x, type = c("interval", "curve"), labels = NULL,
                     series = NULL, title = NULL, subtitle = NULL,
                     monochrome = FALSE, draw = TRUE) {
  type <- match.arg(type)
  for (flag in list(monochrome, draw))
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("monochrome and draw must be single logical values.", call. = FALSE)
  provenance <- data.frame(key = character(), value = character())
  if (inherits(x, "pnsbib_analysis")) {
    provenance <- x$provenance
    x <- x$bounds
  }
  if (!is.data.frame(x) || !nrow(x) ||
      !all(c("lower", "upper", "method", "status") %in% names(x)) ||
      anyDuplicated(names(x)) || any(vapply(x, is.list, logical(1))))
    stop("x needs nonempty result rows with lower, upper, method and status.", call. = FALSE)
  data <- as.data.frame(x)
  n <- nrow(data)
  .poc_plot_text(as.character(data$method), n, "method")
  .poc_plot_text(as.character(data$status), n, "status")
  if (!is.numeric(data$lower) || !is.numeric(data$upper))
    stop("Bounds must be numeric.", call. = FALSE)
  ok <- data$status == "ok"
  finite <- is.finite(data$lower) & is.finite(data$upper)
  if (any(ok & !finite) || any(finite & (data$lower < 0 | data$upper > 1 |
                                        data$lower > data$upper)) ||
      any(xor(is.na(data$lower), is.na(data$upper))) ||
      any(is.infinite(data$lower) | is.infinite(data$upper)))
    stop("Bounds must be ordered probabilities; ok rows require finite endpoints.", call. = FALSE)
  method <- .poc_plot_method(data$method)
  event <- if ("estimand" %in% names(data)) as.character(data$estimand) else
    if ("event" %in% names(data)) as.character(data$event) else paste("Result", seq_len(n))
  event[is.na(event) | !nzchar(event)] <- paste("Result", which(is.na(event) | !nzchar(event)))
  if (is.null(labels)) {
    labels <- paste(event, method, sep = " / ")
    if ("sharpness" %in% names(data)) labels <- paste0(labels, " [", data$sharpness, "]")
    if ("scenario" %in% names(data)) labels <- paste(labels, data$scenario, sep = " / ")
    if ("horizon" %in% names(data)) labels <- paste0(labels, " / h=", data$horizon)
  }
  labels <- .poc_plot_text(labels, n, "labels")
  semantic <- intersect(c("estimand", "method", "sharpness", "regime", "comparison",
    "condition_on", "regime_outcome", "reference_outcome", "outcome_order", "assumptions",
    "scenario", "input_status", "identification", "outcome_source", "population",
    "observed_delta", "path_delta", "time_delta", "minimum_condition_probability"), names(data))
  if (is.null(series)) {
    # Length prefixes prevent separator characters in scientific labels merging groups.
    signature <- do.call(paste0, lapply(data[semantic], function(v) {
      value <- as.character(v)
      ifelse(is.na(value), "N;", paste0(nchar(value), ":", value, ";"))
    }))
    group <- match(signature, unique(signature))
    series <- paste(event, method, sep = " / ")
    for (field in intersect(c("scenario", "input_status", "condition_on"), names(data)))
      series <- paste(series, data[[field]], sep = " / ")
    # Different assumption sets may have identical short labels: keep distinct series.
    names_by_group <- series[match(seq_len(max(group)), group)]
    names_by_group <- make.unique(names_by_group, sep = " #")
    series <- names_by_group[group]
  } else {
    series <- .poc_plot_text(series, n, "series")
    group <- match(series, unique(series))
  }
  horizon <- rep(NA_real_, n)
  if (type == "curve") {
    if (!"estimand" %in% names(data) || !"horizon" %in% names(data) ||
        !is.numeric(data$horizon) || any(!is.finite(data$horizon)))
      stop("Curve plots require estimand and finite numeric horizon columns.", call. = FALSE)
    horizon <- data$horizon
    if (anyDuplicated(data.frame(group = group, horizon = horizon)))
      stop("Duplicate horizons within a series; separate methods or scenarios.", call. = FALSE)
    for (indices in split(seq_len(n), group)) for (field in semantic)
      if (length(unique(data[[field]][indices])) > 1L)
        stop("A curve series mixes ", field, "; use separate series.", call. = FALSE)
    if (length(unique(group)) > 8L)
      stop("Use at most eight curve series per figure; split the result table.", call. = FALSE)
  }
  if (is.null(title)) title <- if (type == "curve") "Probability-of-causation curves" else
    "Probability-of-causation bounds"
  if (is.null(subtitle)) {
    subtitle <- "Probability bounds; interpretation follows the supplied inputs and assumptions."
    for (field in intersect(c("input_status", "identification"), names(data)))
      subtitle <- paste(subtitle, paste(unique(as.character(data[[field]])), collapse = "; "))
  }
  for (value in list(title, subtitle))
    if (!is.character(value) || length(value) != 1L || is.na(value))
      stop("title and subtitle must be single character strings.", call. = FALSE)
  out <- structure(list(data = data, layout = data.frame(source_row = seq_len(n),
    label = labels, series = series, group = group,
    style = if (type == "interval") match(data$method, unique(data$method)) else group,
    horizon = horizon, defined = ok),
    provenance = provenance, type = type, title = title, subtitle = subtitle,
    monochrome = monochrome), class = "pnsbib_plot")
  if (draw) plot(out)
  invisible(out)
}

#' @export
plot.pnsbib_plot <- function(x, ...) {
  if (length(list(...))) stop("Set plot options through poc_plot().", call. = FALSE)
  d <- x$data
  layout <- x$layout
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  palette <- if (x$monochrome) rep("#222222", 8) else
    c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#666666")
  col <- palette[(layout$style - 1L) %% 8L + 1L]
  pch <- c(16, 17, 15, 18, 1, 2, 0, 5)[(layout$style - 1L) %% 8L + 1L]
  note <- strwrap(x$subtitle, width = max(30, floor((grDevices::dev.size("in")[1] - .4) * 14)))
  bottom <- max(1.1, .85 + .15 * length(note))
  if (x$type == "interval") {
    labels <- vapply(layout$label, function(s) paste(strwrap(s, width = 43), collapse = "\n"), "")
    graphics::par(mai = c(bottom, 3.0, .65, .95), mgp = c(2.1, .6, 0), family = "sans")
    y <- rev(seq_len(nrow(d)))
    graphics::plot(NA_real_, NA_real_, xlim = c(0, 1), ylim = c(.35, nrow(d) + .65),
      xaxs = "i", yaxs = "i", xlab = "Probability", ylab = "", axes = FALSE, main = x$title)
    graphics::axis(1, at = seq(0, 1, .2))
    graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = .7)
    graphics::abline(v = seq(0, 1, .2), col = "#E8E8E8", lwd = .7)
    ok <- layout$defined
    graphics::segments(d$lower[ok], y[ok], d$upper[ok], y[ok], col = col[ok], lwd = 2)
    for (endpoint in c("lower", "upper"))
      graphics::points(d[[endpoint]][ok], y[ok], col = col[ok], pch = pch[ok], cex = .65, xpd = NA)
    if (any(!ok)) graphics::text(.02, y[!ok], d$status[!ok], adj = 0, cex = .72, col = "#555555")
    value <- ifelse(ok, sprintf("[%.3f, %.3f]", d$lower, d$upper), "")
    graphics::text(1.035, y, value, adj = 0, xpd = NA, cex = .7)
  } else {
    groups <- unique(layout$group)
    legend_rows <- ceiling(length(groups) / 2)
    graphics::par(mar = c(5.5, 4.5, 3.2 + 1.1 * legend_rows, 1),
                  mgp = c(2.3, .7, 0), family = "sans")
    limits <- range(layout$horizon)
    if (diff(limits) == 0) limits <- limits + c(-.5, .5)
    graphics::plot(NA_real_, NA_real_, xlim = limits, ylim = c(0, 1), yaxs = "i",
      xlab = "Horizon (model period)", ylab = "Probability", axes = FALSE)
    graphics::axis(1, at = sort(unique(layout$horizon)))
    graphics::axis(2, at = seq(0, 1, .2), las = 1)
    graphics::abline(h = seq(0, 1, .2), col = "#E8E8E8", lwd = .7)
    for (g in groups) {
      rows <- which(layout$group == g)
      rows <- rows[order(layout$horizon[rows])]
      lty <- (g - 1L) %% 6L + 1L
      for (endpoint in c("lower", "upper")) {
        values <- d[[endpoint]][rows]
        values[!layout$defined[rows]] <- NA_real_
        graphics::lines(layout$horizon[rows], values, col = col[rows[1]], lty = lty, lwd = 1.6)
        graphics::points(layout$horizon[rows], values, col = col[rows[1]], pch = pch[rows[1]], cex = .65)
      }
    }
    first <- match(groups, layout$group)
    graphics::legend("top", inset = c(0, -.075), xpd = NA, bty = "n", ncol = 2,
      legend = layout$series[first], col = col[first], pch = pch[first],
      lty = (groups - 1L) %% 6L + 1L, cex = .65)
    graphics::title(main = x$title, line = 1.3 + 1.1 * legend_rows)
    if (any(!layout$defined)) {
      note <- paste0("Gaps: ", paste(unique(d$status[!layout$defined]), collapse = ", "))
      graphics::mtext(note, side = 3, line = .1, adj = 0, cex = .65)
    }
  }
  graphics::mtext(paste(note, collapse = "\n"), side = 1, line = 3.3, cex = .65,
                  at = graphics::grconvertX(.5, from = "ndc", to = "user"))
  invisible(x)
}

#' @export
print.pnsbib_plot <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}

#' Save probability-bound graphics and their source tables
#'
#' The filename is a prefix without a PNG/PDF extension. All destinations are
#' checked before writing. Existing files are never overwritten. Each graphics
#' device opened here is closed, including when drawing fails; caller devices
#' are retained. Source data, display mapping, provenance and the complete plot
#' object are exported alongside the figures for reproducibility.
#'
#' @param x A pnsbib_plot returned by [poc_plot()].
#' @param filename Output filename prefix in an existing directory.
#' @param formats One or both of png and pdf; default writes both.
#' @param width,height Figure dimensions in inches. Default height adapts to rows.
#' @param dpi PNG resolution; default 600. PDF remains vector graphics.
#' @return Invisibly, named output paths.
#' @export
poc_save_plot <- function(x, filename, formats = c("png", "pdf"),
                          width = 8, height = NULL, dpi = 600) {
  if (!inherits(x, "pnsbib_plot")) stop("x must be a pnsbib_plot.", call. = FALSE)
  if (!is.character(filename) || length(filename) != 1L || is.na(filename) ||
      !nzchar(basename(filename)) || basename(filename) %in% c(".", "..") ||
      dir.exists(filename) || !dir.exists(dirname(filename)) || grepl("%", filename, fixed = TRUE) ||
      grepl("\\.(png|pdf)$", filename, ignore.case = TRUE))
    stop("filename must be a prefix without an extension in an existing directory.", call. = FALSE)
  if (!is.character(formats) || !length(formats) || anyNA(formats) ||
      any(!formats %in% c("png", "pdf")) || anyDuplicated(formats))
    stop("formats must contain png and/or pdf without duplicates.", call. = FALSE)
  if (is.null(height)) height <- if (x$type == "curve") 5 else max(3.5, 1.7 + .65 * nrow(x$data))
  for (size in list(width, height, dpi))
    if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size <= 0)
      stop("width, height and dpi must be positive finite numbers.", call. = FALSE)
  if (width < 5 || height < 3 || dpi < 72 || dpi > 1200 ||
      width > 30 || height > 60 || ("png" %in% formats && width * height * dpi^2 > 1e8))
    stop("Figure dimensions/resolution are outside readable, bounded export limits.", call. = FALSE)
  suffixes <- c(formats, "data.tsv", "layout.tsv", "provenance.tsv", "plot.rds")
  paths <- stats::setNames(paste0(filename, ".", suffixes), suffixes)
  links <- Sys.readlink(paths)
  if (any(file.exists(paths)) || any(!is.na(links) & nzchar(links)))
    stop("Output files already exist; use a new filename prefix.", call. = FALSE)
  if ("png" %in% formats && !capabilities("png"))
    stop("PNG graphics are unavailable on this R installation.", call. = FALSE)
  render <- function(format, path) {
    before <- grDevices::dev.list()
    caller <- grDevices::dev.cur()
    on.exit({
      # A compiled-in bitmap capability can still fail to load at runtime.
      # Never close the caller's device when no export device was opened.
      added <- setdiff(as.integer(grDevices::dev.list()), as.integer(before))
      for (device in rev(added)) grDevices::dev.off(device)
      if (caller %in% grDevices::dev.list()) invisible(grDevices::dev.set(caller))
    }, add = TRUE)
    png_type <- if (.Platform$OS.type == "windows") "windows" else
      if (identical(Sys.info()[["sysname"]], "Darwin") && capabilities("aqua")) "quartz" else
        if (capabilities("cairo")) "cairo" else getOption("bitmapType")
    withCallingHandlers({
      if (format == "pdf") grDevices::pdf(path, width = width, height = height,
        family = "Helvetica", useDingbats = FALSE, title = x$title) else
        grDevices::png(path, width = width, height = height, units = "in", res = dpi,
                      type = png_type)
    }, warning = function(w) stop(conditionMessage(w), call. = FALSE))
    if (grDevices::dev.cur() == 1L || grDevices::dev.cur() %in% before)
      stop("The graphics device could not be opened.", call. = FALSE)
    plot(x)
  }
  for (format in formats) render(format, paths[[format]])
  for (name in c("data", "layout", "provenance"))
    utils::write.table(x[[name]], paths[[paste0(name, ".tsv")]], sep = "\t",
                      row.names = FALSE, quote = TRUE, na = "NA")
  saveRDS(x, paths[["plot.rds"]], version = 2)
  invisible(paths)
}
