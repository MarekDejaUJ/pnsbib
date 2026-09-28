#' Estimate static margins by a discrete baseline g-formula
#'
#' Computes the weighted observational joint distribution and standardizes
#' each treatment-specific outcome distribution to the same empirical baseline
#' covariate distribution. Every treatment must have positive weight in every
#' target covariate stratum; this function does not extrapolate unsupported
#' treatment outcomes. `assume_identification = TRUE` records the analyst's
#' declaration of consistency, conditional exchangeability given the supplied
#' baseline covariates, positivity, and an acceptable interference model. The
#' data and support checks cannot establish those assumptions. Covariates are
#' treated as discrete categories, including numeric covariates. Factor levels
#' define the modeled treatment and outcome state spaces; encode all
#' prespecified treatment states as factor levels so an unobserved state fails
#' the positivity check instead of silently disappearing.
#'
#' @param data A data frame with one row per observational unit.
#' @param treatment,outcome Names of treatment and outcome columns.
#' @param covariates Names of baseline covariate columns; may be empty only when
#'   unconditional exchangeability is defensible.
#' @param weights `NULL`, a numeric vector, or the name of a weight column.
#'   Weights must be finite and nonnegative, with positive total weight.
#' @param assume_identification Set to `TRUE` only after defending consistency,
#'   conditional exchangeability, positivity and interference for the design.
#' @param check_lp,max_variables Passed to [poc_model()].
#' @return A list containing `model`, `support` (one row per covariate stratum
#'   and treatment), and `diagnostics`. `model$o` is the weighted observational
#'   joint table and `model$a` is the standardized interventional table under
#'   the declared identification assumptions. This function gives point
#'   estimates and support diagnostics; it does not produce confidence
#'   intervals or validate the causal assumptions.
#' @export
poc_gformula <- function(data, treatment, outcome, covariates = character(),
                         weights = NULL, assume_identification = FALSE,
                         check_lp = TRUE, max_variables = 10000L) {
  if (!is.data.frame(data) || !nrow(data) ||
      !is.character(treatment) || length(treatment) != 1L ||
      !is.character(outcome) || length(outcome) != 1L ||
      !is.character(covariates) || anyNA(covariates) ||
      anyNA(c(treatment, outcome, covariates)) ||
      anyDuplicated(c(treatment, outcome, covariates)) ||
      !all(c(treatment, outcome, covariates) %in% names(data))) {
    stop("Specify distinct, existing treatment, outcome and covariate columns in a nonempty data frame.",
         call. = FALSE)
  }
  if (!identical(assume_identification, TRUE)) {
    stop("Set assume_identification = TRUE only after defending the causal identification assumptions.",
         call. = FALSE)
  }
  columns <- c(treatment, outcome, covariates)
  if (any(vapply(data[columns], function(x) anyNA(x) ||
                 any(!nzchar(as.character(x))), logical(1)))) {
    stop("Treatment, outcome and covariate columns must have no missing or empty values.",
         call. = FALSE)
  }
  if (is.null(weights)) {
    w <- rep(1, nrow(data))
  } else if (is.character(weights) && length(weights) == 1L &&
             !is.na(weights) && weights %in% names(data)) {
    w <- data[[weights]]
  } else {
    w <- weights
  }
  if (!is.numeric(w) || length(w) != nrow(data) || anyNA(w) ||
      any(!is.finite(w)) || any(w < 0) || !is.finite(sum(w)) ||
      sum(w) <= 0) {
    stop("weights must be a nonnegative finite vector with positive total weight.",
         call. = FALSE)
  }
  x <- as.character(data[[treatment]])
  y <- as.character(data[[outcome]])
  xs <- if (is.factor(data[[treatment]])) levels(data[[treatment]]) else sort(unique(x))
  ys <- if (is.factor(data[[outcome]])) levels(data[[outcome]]) else sort(unique(y))
  if (length(xs) < 2L || length(ys) < 2L ||
      any(!nzchar(xs)) || any(!nzchar(ys))) {
    stop("Treatment and outcome must each have at least two nonempty categories.",
         call. = FALSE)
  }
  if (length(covariates)) {
    factors <- lapply(data[covariates], function(z) factor(as.character(z)))
    stratum <- as.integer(do.call(interaction,
                                 c(factors, list(drop = TRUE, lex.order = TRUE))))
  } else {
    stratum <- rep.int(1L, nrow(data))
  }
  ids <- sort(unique(stratum[w > 0]))
  total <- sum(w)
  o <- matrix(0, length(xs), length(ys), dimnames = list(xs, ys))
  a <- o
  for (i in seq_along(xs)) {
    for (j in seq_along(ys)) {
      o[i, j] <- sum(w[x == xs[i] & y == ys[j]]) / total
    }
  }
  support <- vector("list", length(ids) * length(xs))
  k <- 0L
  for (s in ids) {
    in_stratum <- stratum == s & w > 0
    stratum_weight <- sum(w[in_stratum])
    for (i in seq_along(xs)) {
      k <- k + 1L
      selected <- in_stratum & x == xs[i]
      selected_weight <- sum(w[selected])
      row <- data.frame(stratum = s, treatment = xs[i],
                        n = sum(selected), weight = selected_weight,
                        target_weight = stratum_weight,
                        propensity = selected_weight / stratum_weight,
                        stringsAsFactors = FALSE)
      if (length(covariates)) {
        row <- cbind(row, data[which(in_stratum)[1L], covariates,
                               drop = FALSE])
      }
      support[[k]] <- row
      if (selected_weight > 0) {
        for (j in seq_along(ys)) {
          a[i, j] <- a[i, j] + stratum_weight / total *
            sum(w[selected & y == ys[j]]) / selected_weight
        }
      }
    }
  }
  support <- do.call(rbind, support)
  rownames(support) <- NULL
  if (any(support$weight == 0)) {
    failing <- support[support$weight == 0, c("stratum", "treatment")]
    detail <- paste(sprintf("%s/%s", failing$stratum, failing$treatment),
                    collapse = ", ")
    stop(sprintf("Positivity fails in stratum/treatment: %s.", detail),
         call. = FALSE)
  }
  model <- poc_model(o, a, check_lp = check_lp,
                     max_variables = max_variables)
  list(model = model, support = support,
       diagnostics = list(n = nrow(data), positive_weight_n = sum(w > 0),
                          effective_n = 1 / sum((w / total)^2),
                          strata = length(ids),
                          minimum_propensity = min(support$propensity),
                          identification_assumed = TRUE,
                          formula = "sum_l P(Y=y|X=x,L=l) P(L=l)"))
}
