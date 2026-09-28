fixture_from_response_types <- function(m = 3L, n = 3L, seed = 19L) {
  set.seed(seed)
  types <- as.matrix(do.call(expand.grid, rep(list(seq_len(n)), m)))
  q <- matrix(rexp(m * nrow(types)), nrow = m)
  q <- q / sum(q)
  o <- matrix(0, m, n, dimnames = list(paste0("x", seq_len(m)),
                                       paste0("y", seq_len(n))))
  a <- o
  for (x in seq_len(m)) for (y in seq_len(n)) {
    o[x, y] <- sum(q[x, types[, x] == y])
    a[x, y] <- sum(q[, types[, x] == y])
  }
  list(o = o, a = a, q = q, types = types)
}

fixture_truth <- function(fixture, query) {
  x <- match(names(query$counterfactual), rownames(fixture$o))
  y <- match(unname(query$counterfactual), colnames(fixture$o))
  factual_x <- seq_len(nrow(fixture$o))
  total <- 0
  for (z in factual_x) for (r in seq_len(nrow(fixture$types))) {
    matches <- all(fixture$types[r, x] == y)
    if (!is.null(query$observed_x)) {
      matches <- matches && rownames(fixture$o)[z] == query$observed_x
    }
    if (!is.null(query$observed_y)) {
      matches <- matches && colnames(fixture$o)[fixture$types[r, z]] == query$observed_y
    }
    if (matches) total <- total + fixture$q[z, r]
  }
  if (query$conditional) {
    den <- if (!is.null(query$observed_x) && !is.null(query$observed_y))
      fixture$o[query$observed_x, query$observed_y] else
      if (!is.null(query$observed_x)) sum(fixture$o[query$observed_x, ]) else
        sum(fixture$o[, query$observed_y])
    total <- total / den
  }
  total
}
