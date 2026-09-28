# Runnable companion to own-data.md. Synthetic examples, no API or installation.
# R prepares data, asserts results and draws; LP calls explicitly request Zig.
library(pnsbib)

run_pnsbib_examples <- function(output = tempfile("pnsbib-own-data-"), graphics = TRUE) {
  stopifnot(!file.exists(output)); dir.create(output, recursive = TRUE)
  # 1. Replace these synthetic rows with your independently prepared data.
  d <- data.frame(
    treatment = rep(c("reference", "reference", "active", "active"), c(35, 15, 20, 30)),
    outcome = rep(c("failure", "success", "failure", "success"), c(35, 15, 20, 30)))
  d$treatment <- factor(d$treatment, levels = c("reference", "active"))
  d$outcome <- factor(d$outcome, levels = c("failure", "success"))
  stopifnot(!anyNA(d), nrow(d) > 0L)
  o <- unclass(table(d$treatment, d$outcome)) / nrow(d)
  dimnames(o) <- list(levels(d$treatment), levels(d$outcome))
  stopifnot(abs(sum(o) - 1) < 1e-12)
  a <- rbind(reference = c(failure = .70, success = .30),
             active = c(failure = .55, success = .45))
  # These are supplied hypothetical intervention risks, not observed row risks.
  m <- poc_model(o, a, lp_backend = "zig")
  validated <- validate_poc_inputs(o, a, lp_backend = "zig")
  stopifnot(validated$checked_lp)
  q <- list(PNS = poc_query(c(reference = "failure", active = "success")),
    PN = poc_query(c(reference = "failure"), "active", "success", conditional = TRUE),
    PS = poc_query(c(active = "success"), "reference", "failure", conditional = TRUE))
  fits <- lapply(q, function(query) poc_exact(m, query))
  stopifnot(abs(fits$PNS$lower - .15) < 1e-8, abs(fits$PNS$upper - .45) < 1e-8)
  comparison <- list(poc_bounds(m, q$PNS), poc_shu2026(m, q$PNS), fits$PNS)
  stopifnot(all(vapply(comparison, function(z) z$status == "ok", logical(1))))
  ref <- poc_exact(m, q$PNS, lp_backend = "reference")
  stopifnot(max(abs(c(ref$lower, ref$upper) - c(fits$PNS$lower, fits$PNS$upper))) < 1e-8)
  # Equivalent complete-margin long-table input.
  long <- expand.grid(treatment = rownames(o), outcome = colnames(o), stringsAsFactors = FALSE)
  long$observed_joint <- as.vector(o); long$interventional <- as.vector(a)
  stopifnot(identical(poc_model(long, lp_backend = "zig")$o, m$o))

  # 2. Retain missing intervention information instead of filling it with data.
  free <- poc_partial_model(o, lp_backend = "zig")
  free_bounds <- poc_partial_bounds(free, q$PNS, return_witness = TRUE)
  cs <- list(
    poc_constraint(poc_query(c(reference = "success")), .25, .35, label = "reference_risk_band"),
    poc_constraint(poc_query(c(active = "success")), .40, .50, label = "active_risk_band"))
  partial <- poc_partial_model(o, cs, lp_backend = "zig")
  partial_bounds <- poc_partial_bounds(partial, q$PNS, return_witness = TRUE)
  stopifnot(partial_bounds$lower >= free_bounds$lower - 1e-8,
            partial_bounds$upper <= free_bounds$upper + 1e-8)
  subgroup <- poc_constraint(poc_query(c(reference = "success"),
    observed_x = "active", conditional = TRUE), .2, .4, label = "active_subgroup_missing_risk")
  subgroup_fit <- poc_partial_bounds(poc_partial_model(o, list(subgroup), lp_backend = "zig"), q$PN)
  strata <- list(baseline_A = .4 * o, baseline_B = .6 * o)
  sc <- poc_constraint(poc_query(c(active = "success")), .40, .50, stratum = "baseline_A")
  sm <- poc_partial_model(strata, list(sc), lp_backend = "zig")
  stratum_fit <- poc_partial_bounds(sm, q$PNS, stratum = "baseline_A")

  # 3. Explicit categorical states and several counterfactual requirements.
  mo <- matrix(1/9, 3, 3, dimnames = list(c("off", "on", "alternative"), c("low", "middle", "high")))
  mm <- poc_model(mo, mo * 3, lp_backend = "zig")
  mq <- poc_query(c(off = "low", on = "high", alternative = "middle"))
  multivalued <- list(poc_bounds(mm, mq), poc_shu2026(mm, mq), poc_exact(mm, mq))

  # 4. Fixed input bands are sensitivity assumptions, not confidence intervals.
  sensitivity <- poc_sensitivity(m, q$PNS, delta = c(0, .025, .05, .1))
  input_band <- poc_input_sensitivity(m, q$PN, observed_delta = .02,
    interventional_delta = .05, minimum_condition_probability = .1)

  # 5. Complete factual path cells; all model periods have the declared meaning.
  regimes <- rbind(off = c(0, 0), on = c(1, 1))
  paths <- expand.grid(history = c("off", "on"), path = c("00", "01", "11"), stringsAsFactors = FALSE)
  paths$prob <- c(.2, .1, .1, .1, .2, .3)
  margins <- expand.grid(regime = c("off", "on"), horizon = 1:2, stringsAsFactors = FALSE)
  margins$prob1 <- c(.4, .6, .6, .8)
  lm <- poc_longitudinal_model(paths, regimes, time_margins = margins, absorbing = TRUE, lp_backend = "zig")
  lq <- poc_longitudinal_query("pns", "on", "off", 2)
  horizon <- poc_longitudinal_bounds(lm, lq)
  curve <- poc_curve(lm, "on", "off")
  band_curve <- poc_curve(lm, "on", "off", bands = list(time_delta = .05))
  lband <- poc_longitudinal_sensitivity(lm, lq, time_delta = .05)
  other_queries <- list(
    poc_longitudinal_query("persistent", "on", "off", 2, start = 1),
    poc_longitudinal_query("event_time", "on", "off", 2),
    poc_longitudinal_query("lagged", "on", "off", 2, lag = 1),
    poc_longitudinal_query("trajectory", "on", "off", 2, regime_path = "11", reference_path = "00"))
  other_paths <- lapply(other_queries, function(query) poc_longitudinal_bounds(lm, query))
  ip <- data.frame(regime = paths$history, path = paths$path, prob = paths$prob * 2)
  path_model <- poc_longitudinal_model(paths, regimes, interventional_paths = ip, absorbing = TRUE, lp_backend = "zig")
  path_band <- poc_longitudinal_sensitivity(path_model, lq, path_delta = .05, observed_delta = .01)

  # Prefix conditioning pools factual continuations sharing the treatment prefix.
  pr <- rbind(off = c(0, 0), later = c(0, 1), early = c(1, 0), always = c(1, 1))
  po <- expand.grid(history = rownames(pr), path = c("00", "11"), stringsAsFactors = FALSE)
  po$prob <- c(.1, .2, 0, 0, 0, 0, .3, .4)
  pm <- poc_longitudinal_model(po, pr, allowed_paths = c("00", "11"), lp_backend = "zig")
  prefix <- poc_longitudinal_bounds(pm, poc_longitudinal_query("pn", "early", "off", 1, condition_on = "prefix"))

  # Named categorical potential paths and explicitly ordered absorption.
  op <- rbind(LL = c("low", "low"), LM = c("low", "middle"),
              LH = c("low", "high"), MM = c("middle", "middle"),
              MH = c("middle", "high"), HH = c("high", "high"))
  co <- expand.grid(history = c("off", "on"), path = rownames(op), stringsAsFactors = FALSE)
  co$prob <- 1/nrow(co)
  state_margins <- data.frame(regime = c("off", "on"), horizon = 2,
                             outcome = "high", prob = .5)
  cm <- poc_longitudinal_model(co, regimes, time_margins = state_margins, outcome_paths = op, absorbing = TRUE,
    outcome_order = c("low", "middle", "high"), lp_backend = "zig")
  cq <- poc_longitudinal_query("pns", "on", "off", 2,
    regime_outcome = c("middle", "high"), reference_outcome = "low")
  categorical <- poc_longitudinal_bounds(cm, cq)

  # 6. Sampling vectors are one contribution per independent unit, not counts.
  map <- poc_margin_map(m)
  samples <- list()
  for(i in which(map$type == "observed")) {
    # Use explicit map labels, never rely on table storage order.
    samples[[map$cell_id[i]]] <- as.numeric(d$treatment == map$treatment[i] & d$outcome == map$outcome[i])
  }
  stopifnot(all(lengths(samples) == nrow(d)), all(rowSums(do.call(cbind, samples)) == 1))
  region <- poc_confidence_region(m, samples, target = "Illustrative iid target",
    sampling_unit = "individual", method = "clopper_pearson", assume_sampling = TRUE)
  confidence <- poc_confidence_bounds(region, q$PNS, witnesses = TRUE, lp_backend = "zig")
  # Hypothetical independent, equal-weight clusters of ten units each.
  clusters <- lapply(samples, function(v) colMeans(matrix(v, nrow = 10)))
  cluster_region <- poc_confidence_region(m, clusters, target = "Illustrative equal-cluster target",
    sampling_unit = "cluster", method = "hoeffding", assume_sampling = TRUE)
  cluster_confidence <- poc_confidence_bounds(cluster_region, q$PNS, lp_backend = "zig")

  # 7. Optional baseline discrete standardization, with external identification.
  gd <- expand.grid(treatment = c("reference", "active"), outcome = c("failure", "success"),
                    baseline = c("A", "B"), replicate = 1:10, stringsAsFactors = FALSE)
  gd$treatment <- factor(gd$treatment, levels = c("reference", "active"))
  gd$outcome <- factor(gd$outcome, levels = c("failure", "success"))
  gf <- poc_gformula(gd, "treatment", "outcome", covariates = "baseline", assume_identification = TRUE)
  adjusted <- poc_exact(gf$model, q$PNS, lp_backend = "zig")

  # 8. Reports and graphics retain source probabilities and statuses.
  report <- poc_analysis(comparison, observed = long,
    sensitivity = sensitivity, provenance = list(target = "Synthetic tutorial", margins = "Declared hypothetical risks"))
  summary(report)
  exports <- poc_export_summary(report, file.path(output, "tables"))
  plot_files <- list()
  if(graphics) {
    fig <- poc_plot(report, draw = FALSE, monochrome = TRUE)
    plot_files$interval <- poc_save_plot(fig, file.path(output, "interval"), width = 6, dpi = 150)
    cf <- poc_plot(curve, type = "curve", draw = FALSE, monochrome = TRUE)
    plot_files$curve <- poc_save_plot(cf, file.path(output, "curve"), width = 6, dpi = 150)
  }
  all_fits <- c(fits, comparison, list(free_bounds, partial_bounds, subgroup_fit, stratum_fit),
    multivalued, list(input_band, horizon, lband), other_paths,
    list(path_band, prefix, categorical, confidence, cluster_confidence, adjusted))
  stopifnot(all(vapply(all_fits, function(z) z$status == "ok", logical(1))),
    nrow(curve) == 6L, nrow(band_curve) == 6L)
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
  list(output = output, observed = o, complete = fits, partial = partial_bounds,
       curve = curve, confidence = confidence, gformula = gf,
       interval_count = length(all_fits), exports = exports, plot_files = plot_files)
}

pnsbib_examples <- run_pnsbib_examples()
cat("Own-data examples passed; outputs:", pnsbib_examples$output, "\n")
