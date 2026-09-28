.poc_result <- function(lower, upper, method, sharpness, status,
                        denominator, query) {
  structure(data.frame(lower = lower, upper = upper, method = method,
             sharpness = sharpness, status = status,
             denominator = denominator,
             event = paste(paste0("Y_", names(query$counterfactual), "=",
                                  unname(query$counterfactual)), collapse = " & "),
             observed_x = if (is.null(query$observed_x)) NA_character_ else query$observed_x,
             observed_y = if (is.null(query$observed_y)) NA_character_ else query$observed_y,
             conditional = query$conditional, stringsAsFactors = FALSE),
            class = c("pnsbib_result", "data.frame"))
}

.poc_single_r <- function(kind, aji, oji, vi, vk, uj, up, opk, ojk, other_ok) {
  if (kind == 4L) return(c(max(oji, aji + vi - 1), min(aji, vi)))
  if (kind == 5L) {
    partition <- sum(pmax(0, aji + other_ok - 1 + uj - oji))
    return(c(max(0, aji + vk - 1, partition),
             min(aji - oji, vk - ojk)))
  }
  if (kind == 6L) {
    return(c(max(0, aji - oji - 1 + uj + up), min(aji - oji, up)))
  }
  c(max(0, aji + opk - 1 + uj - oji), min(aji - oji, opk))
}

.poc_single <- function(model, kind, j, i, p = NA_integer_,
                        k = NA_integer_, use_zig = TRUE) {
  o <- model$o
  a <- model$a
  other_ok <- if (kind == 5L) o[-j, k] else numeric()
  args <- list(as.integer(kind), a[j, i], o[j, i], model$py[i],
               if (is.na(k)) 0 else model$py[k], model$px[j],
               if (is.na(p)) 0 else model$px[p],
               if (is.na(p) || is.na(k)) 0 else o[p, k],
               if (is.na(k)) 0 else o[j, k],
               as.numeric(other_ok))
  if (isTRUE(use_zig)) {
    return(do.call(poc_single_zig, args))
  }
  do.call(.poc_single_r, args)
}

.poc_paper_interval <- function(model, x, y, ox = NA_integer_, oy = NA_integer_,
                                cache = new.env(parent = emptyenv()),
                                use_zig = TRUE) {
  key <- paste(paste(x, y, sep = ":", collapse = ","), ox, oy, sep = "|")
  if (exists(key, envir = cache, inherits = FALSE)) return(get(key, cache))
  o <- model$o
  a <- model$a
  m <- nrow(o)
  k <- length(x)
  if (k == 0L) {
    answer <- if (!is.na(ox) && !is.na(oy)) c(o[ox, oy], o[ox, oy]) else
      if (!is.na(ox)) c(model$px[ox], model$px[ox]) else
      if (!is.na(oy)) c(model$py[oy], model$py[oy]) else c(1, 1)
  } else if (!is.na(ox) && ox %in% x) {
    t <- match(ox, x)
    if (!is.na(oy) && oy != y[t]) {
      answer <- c(0, 0)
    } else {
      answer <- .poc_paper_interval(model, x[-t], y[-t], ox, y[t], cache, use_zig)
    }
  } else if (k == 1L) {
    j <- x[1L]
    i <- y[1L]
    if (is.na(ox) && is.na(oy)) {
      answer <- c(a[j, i], a[j, i])
    } else {
      kind <- if (is.na(ox)) if (oy == i) 4L else 5L else
        if (is.na(oy)) 6L else 7L
      answer <- .poc_single(model, kind, j, i, ox, oy, use_zig)
    }
  } else {
    margins <- a[cbind(x, y)]
    reduced <- lapply(seq_len(k), function(t)
      .poc_paper_interval(model, x[-t], y[-t], cache = cache, use_zig = use_zig))
    lower_candidates <- 0
    upper_candidates <- c(margins, vapply(reduced, function(z) z[2L], 0.0))
    if (is.na(ox) && is.na(oy)) {
      lower_candidates <- c(lower_candidates, sum(margins) - k + 1,
        vapply(seq_len(k), function(t) reduced[[t]][1L] + margins[t] - 1, 0.0))
      by_x <- lapply(seq_len(m), function(p) {
        t <- match(p, x)
        if (!is.na(t)) {
          .poc_paper_interval(model, x[-t], y[-t], p, y[t], cache, use_zig)
        } else {
          .poc_paper_interval(model, x, y, p, NA_integer_, cache, use_zig)
        }
      })
      lower_candidates <- c(lower_candidates,
                            sum(vapply(by_x, function(z) z[1L], 0.0)))
      upper_candidates <- c(upper_candidates,
                            sum(vapply(by_x, function(z) z[2L], 0.0)))
    } else if (!is.na(ox) && is.na(oy)) {
      singles <- lapply(seq_len(k), function(t)
        .poc_paper_interval(model, x[t], y[t], ox, NA_integer_, cache, use_zig))
      lower_candidates <- c(lower_candidates, model$px[ox] - k + sum(margins),
        vapply(seq_len(k), function(t)
          reduced[[t]][1L] + singles[[t]][1L] - 1, 0.0))
      upper_candidates <- c(upper_candidates, model$px[ox],
        vapply(singles, function(z) z[2L], 0.0))
    } else if (is.na(ox) && !is.na(oy)) {
      singles <- lapply(seq_len(k), function(t)
        .poc_paper_interval(model, x[t], y[t], NA_integer_, oy, cache, use_zig))
      lower_candidates <- c(lower_candidates, model$py[oy] - k + sum(margins),
        vapply(seq_len(k), function(t)
          reduced[[t]][1L] + singles[[t]][1L] - 1, 0.0))
      by_x <- lapply(seq_len(m), function(p) {
        t <- match(p, x)
        if (!is.na(t)) {
          if (y[t] != oy) return(c(0, 0))
          .poc_paper_interval(model, x[-t], y[-t], p, oy, cache, use_zig)
        } else {
          .poc_paper_interval(model, x, y, p, oy, cache, use_zig)
        }
      })
      lower_candidates <- c(lower_candidates,
                            sum(vapply(by_x, function(z) z[1L], 0.0)))
      upper_candidates <- c(upper_candidates, model$py[oy],
        vapply(singles, function(z) z[2L], 0.0),
        sum(vapply(by_x, function(z) z[2L], 0.0)))
    } else {
      singles <- lapply(seq_len(k), function(t)
        .poc_paper_interval(model, x[t], y[t], ox, oy, cache, use_zig))
      lower_candidates <- c(lower_candidates, o[ox, oy] - k + sum(margins),
        vapply(seq_len(k), function(t)
          reduced[[t]][1L] + singles[[t]][1L] - 1, 0.0))
      upper_candidates <- c(upper_candidates, o[ox, oy],
        vapply(singles, function(z) z[2L], 0.0))
    }
    answer <- c(max(lower_candidates), min(upper_candidates))
  }
  if (answer[1L] > answer[2L] + model$tolerance ||
      answer[1L] < -model$tolerance ||
      answer[2L] > 1 + model$tolerance) {
    stop("Paper bounds are inconsistent for the supplied margins.", call. = FALSE)
  }
  answer <- pmin(1, pmax(0, answer))
  assign(key, answer, envir = cache)
  answer
}

#' Li-Pearl bounds for multivalued probabilities of causation
#'
#' Implements the printed Theorems 4–11 of Li and Pearl (2024). Theorem 4
#' and 5 intervals are labelled valid_analytical: they can be wider than the
#' exact LP under the supplied complete margins. Theorem 6 and 7 intervals
#' are sharp. Multi-term recursive intervals are labelled valid_recursive.
#'
#' @param model A model returned by [poc_model()].
#' @param query A query returned by [poc_query()].
#' @param method "paper" or "exact".
#' @param use_zig Use the complete Zig analytical recursion. FALSE uses the
#'   preserved R oracle. The native memo is limited to 1048576 states,
#'   `(nrow(model$o)+1)*(ncol(model$o)+1)*2^length(query$counterfactual)`.
#' @return A one-row data frame.
#' @export
poc_bounds <- function(model, query, method = c("paper", "exact"),
                       use_zig = TRUE) {
  method <- match.arg(method)
  if (method == "exact") return(poc_exact(model, query))
  idx <- .poc_indices(model, query)
  denominator <- .poc_denominator(model, query, idx)
  theorem <- .poc_theorem(idx)
  sharpness <- if (theorem %in% c("4", "5")) "valid_analytical" else
    if (theorem %in% c("6", "7", "margin", "consistency")) "sharp" else
      "valid_recursive"
  if (query$conditional && denominator <= model$tolerance) {
    return(.poc_result(NA_real_, NA_real_, theorem, sharpness,
                       "undefined_condition", denominator, query))
  }
  bounds <- if (isTRUE(use_zig)) {
    poc_paper_zig(nrow(model$o), ncol(model$o), as.numeric(model$o),
      as.numeric(model$a), as.integer(idx$x), as.integer(idx$y),
      if (is.na(idx$ox)) 0L else as.integer(idx$ox),
      if (is.na(idx$oy)) 0L else as.integer(idx$oy), model$tolerance)
  } else {
    .poc_paper_interval(model, idx$x, idx$y, idx$ox, idx$oy, use_zig = FALSE)
  }
  .poc_result(bounds[1L] / denominator, bounds[2L] / denominator,
              theorem, sharpness, "ok", denominator, query)
}
