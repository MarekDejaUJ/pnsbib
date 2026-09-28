.poc_shu_subsets <- function(d) {
  if (length(d) < 2L) return(numeric())
  # For a fixed cardinality, the smallest sum uses the smallest distinct d's.
  cumsum(sort(d))[-1L] / seq_len(length(d) - 1L)
}

.poc_shu_interval <- function(model, x, y, ox, oy) {
  o <- model$o
  if (!is.na(ox) && ox %in% x) {
    i <- match(ox, x)
    if (!is.na(oy) && oy != y[i])
      return(list(bounds = c(0, 0), family = "consistency contradiction"))
    reduced <- .poc_shu_interval(model, x[-i], y[-i], ox, y[i])
    reduced$family <- paste0(reduced$family, " after consistency")
    return(reduced)
  }
  k <- length(x)
  if (!k) {
    mass <- if (!is.na(ox) && !is.na(oy)) o[ox, oy] else
      if (!is.na(ox)) model$px[ox] else
        if (!is.na(oy)) model$py[oy] else 1
    return(list(bounds = rep(unname(mass), 2L), family = "factual margin"))
  }
  a <- model$a[cbind(x, y)]
  own <- o[cbind(x, y)]
  d <- a - own
  b <- d + model$px[x]
  outside <- setdiff(seq_len(nrow(o)), x)
  if (!is.na(ox)) {
    mass <- if (is.na(oy)) model$px[ox] else o[ox, oy]
    lower <- c(0, sum(b) + mass - k)
    upper <- c(mass, d)
    family <- if (is.na(oy)) "theorem 3" else "theorem 5"
  } else if (is.na(oy)) {
    lower <- c(0, sum(a) - k + 1, sum(b) - b + own - k + 1)
    upper <- c(sum(own) + sum(model$px[outside]), a, .poc_shu_subsets(d))
    family <- "theorem 2"
  } else {
    matching <- y == oy
    mass <- sum(o[c(outside, x[matching]), oy])
    lower <- c(0, sum(b) + mass - k,
               sum(b) - b[matching] + own[matching] - k + 1)
    upper <- c(mass, a[matching], d[!matching],
               .poc_shu_subsets(d[matching]))
    family <- if (sum(matching) > 1L)
      "theorem 4 + project repeated-outcome extension" else "theorem 4"
  }
  list(bounds = c(max(lower), min(upper)), family = family)
}

#' Shu-Wang-Li closed-form probability-of-causation comparator
#'
#' Evaluates Theorems 2-5 of Shu, Wang and Li (IJCAI 2026) on complete
#' observational and interventional margins. Distinct-subset upper candidates
#' are computed by sorted prefix sums, without enumerating response types.
#' Repeated matches between potential and factual outcome labels use an
#' explicitly labelled project extension of Theorem 4: compatible factual
#' cells, anchor union bounds and distinct matching-subset upper bounds.
#' All intervals are labelled valid_closed_form; general sharpness is not
#' asserted. Numerical agreement with an exact LP is a case-specific check.
#'
#' @param model A complete-margin model returned by [poc_model()].
#' @param query A static query returned by [poc_query()].
#' @param use_zig Use the native closed form. FALSE selects the preserved R
#'   oracle. The native input table is limited to 1048576 cells.
#' @return A one-row pnsbib_result with the theorem/extension in method.
#' @references Shu, X., Wang, S., and Li, A. (2026). Identification of
#'   Probabilities of Causation: From Recursive to Closed-Form Bounds.
#'   IJCAI, 4063-4070. \doi{10.24963/ijcai.2026/452}.
#' @export
poc_shu2026 <- function(model, query, use_zig = TRUE) {
  idx <- .poc_indices(model, query)
  denominator <- unname(.poc_denominator(model, query, idx))
  answer <- if (isTRUE(use_zig)) {
    native <- poc_shu_zig(nrow(model$o), ncol(model$o), as.numeric(model$o),
      as.numeric(model$a), as.integer(idx$x), as.integer(idx$y),
      if (is.na(idx$ox)) 0L else as.integer(idx$ox),
      if (is.na(idx$oy)) 0L else as.integer(idx$oy), model$tolerance)
    families <- c("factual margin", "consistency contradiction", "theorem 2", "theorem 3",
      "theorem 4", "theorem 5", "theorem 4 + project repeated-outcome extension")
    list(bounds = native[1:2], family = paste0(families[native[3L] + 1L],
      if (native[4L]) " after consistency" else ""))
  } else {
    .poc_shu_interval(model, idx$x, idx$y, idx$ox, idx$oy)
  }
  method <- paste0("Shu-Wang-Li 2026: ", answer$family)
  if (query$conditional && denominator <= model$tolerance)
    return(.poc_result(NA_real_, NA_real_, method, "valid_closed_form",
                       "undefined_condition", denominator, query))
  bounds <- answer$bounds
  if (bounds[1L] > bounds[2L] + model$tolerance ||
      bounds[1L] < -model$tolerance || bounds[2L] > 1 + model$tolerance)
    stop("Closed-form bounds are inconsistent for the supplied margins.", call. = FALSE)
  bounds <- pmin(1, pmax(0, bounds))
  .poc_result(bounds[1L] / denominator, bounds[2L] / denominator,
              method, "valid_closed_form", "ok", denominator, query)
}
