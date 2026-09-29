# Supplementary documentation retained when native bindings are regenerated.

#' @name poc_partial_model
#' @rdname poc_partial_model
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
NULL

#' @name poc_constraint
#' @rdname poc_constraint
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' q <- poc_query(c(reference = "success"))
#' poc_constraint(q, lower = .25, upper = .35, label = "reference_risk")
NULL

#' @name poc_longitudinal_bounds
#' @rdname poc_longitudinal_bounds
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' regimes <- rbind(off = c(0, 0), on = c(1, 1))
#' observed <- expand.grid(history = c("off", "on"),
#'                         path = c("00", "01", "11"), stringsAsFactors = FALSE)
#' observed$prob <- c(.2, .1, .1, .1, .2, .3)
#' m <- poc_longitudinal_model(observed, regimes, absorbing = TRUE)
#' q <- poc_longitudinal_query("pns", "on", "off", 2)
#' poc_longitudinal_bounds(m, q, lp_backend = "zig")
NULL

#' @name poc_shu_zig
#' @rdname poc_shu_zig
#' @references
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' poc_shu_zig(2L, 2L, as.numeric(o), as.numeric(a),
#'             c(1L, 2L), c(1L, 2L), 0L, 0L, 1e-9)
NULL

#' @name validate_poc_inputs
#' @rdname validate_poc_inputs
#' @references
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' validate_poc_inputs(o, a)
#' poc_model(o, a, lp_backend = "zig")
NULL

#' @name poc_paper_zig
#' @rdname poc_paper_zig
#' @references
#' Li, A. and Pearl, J. (2024). Probabilities of Causation with Nonbinary
#' Treatment and Effect. \emph{Proceedings of the AAAI Conference on Artificial
#' Intelligence}, 38(18), 20465--20472. \doi{10.1609/aaai.v38i18.30030}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' poc_paper_zig(2L, 2L, as.numeric(o), as.numeric(a),
#'               c(1L, 2L), c(1L, 2L), 0L, 0L, 1e-9)
NULL

#' @name poc_export_summary
#' @rdname poc_export_summary
#' @references
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' report <- poc_analysis(poc_exact(m, q), provenance = list(target = "Synthetic population"))
#' out <- tempfile("pnsbib-report-")
#' poc_export_summary(report, out)
#' unlink(out, recursive = TRUE)
NULL

#' @name poc_input_sensitivity
#' @rdname poc_input_sensitivity
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' poc_input_sensitivity(m, q, observed_delta = .01, interventional_delta = .05)
NULL

#' @name poc_static_event_zig
#' @rdname poc_static_event_zig
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' poc_static_event_zig(2L, 2L, c(1L, 2L, 2L, 1L),
#'   c(1L, 2L), c(1L, 1L), 1L, 2L, 0L, 0L)
NULL

#' @name poc_longitudinal_model
#' @rdname poc_longitudinal_model
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' regimes <- rbind(off = c(0, 0), on = c(1, 1))
#' observed <- expand.grid(history = c("off", "on"),
#'                         path = c("00", "01", "11"), stringsAsFactors = FALSE)
#' observed$prob <- c(.2, .1, .1, .1, .2, .3)
#' m <- poc_longitudinal_model(observed, regimes, absorbing = TRUE)
#' q <- poc_longitudinal_query("pns", "on", "off", 2)
#' poc_longitudinal_bounds(m, q)
NULL

#' @name poc_longitudinal_sensitivity
#' @rdname poc_longitudinal_sensitivity
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
NULL

#' @name poc_bounds
#' @rdname poc_bounds
#' @references
#' Li, A. and Pearl, J. (2024). Probabilities of Causation with Nonbinary
#' Treatment and Effect. \emph{Proceedings of the AAAI Conference on Artificial
#' Intelligence}, 38(18), 20465--20472. \doi{10.1609/aaai.v38i18.30030}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' poc_bounds(m, q)
#' poc_bounds(m, q, method = "exact")
NULL

#' @name poc_gformula
#' @rdname poc_gformula
#' @references
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' d <- expand.grid(treatment = c("reference", "active"),
#'                  outcome = c("failure", "success"), baseline = c("A", "B"))
#' # Identification is an external design assumption, not a result of this check.
#' g <- poc_gformula(d, "treatment", "outcome", covariates = "baseline",
#'                  assume_identification = TRUE)
#' g$model$a
#' g$support
NULL

#' @name poc_query
#' @rdname poc_query
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' poc_query(c(reference = "failure", active = "success"))
#' # PN conditions on factual active-policy success.
#' poc_query(c(reference = "failure"), observed_x = "active",
#'           observed_y = "success", conditional = TRUE)
NULL

#' @name poc_partial_bounds
#' @rdname poc_partial_bounds
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' p <- poc_partial_model(o, list(
#'   poc_constraint(poc_query(c(active = "success")), .40, .50)))
#' poc_partial_bounds(p, q, return_witness = TRUE)
NULL

#' @name poc_longitudinal_query
#' @rdname poc_longitudinal_query
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
NULL

#' @name poc_lp_zig
#' @rdname poc_lp_zig
#' @references
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' # Maximize x1 + 2*x2 subject to x1 + x2 <= 1 and x >= 0.
#' poc_lp_zig(1L, 2L, c(1, 1), 1, -1L, c(1, 2), TRUE, 1000L)[1:6]
NULL

#' @name poc_single_zig
#' @rdname poc_single_zig
#' @references
#' Li, A. and Pearl, J. (2024). Probabilities of Causation with Nonbinary
#' Treatment and Effect. \emph{Proceedings of the AAAI Conference on Artificial
#' Intelligence}, 38(18), 20465--20472. \doi{10.1609/aaai.v38i18.30030}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' # Theorem 4 for a binary synthetic table.
#' poc_single_zig(4L, .3, .15, .45, .45, .5, .5, .15, .15, .3)
NULL

#' @name poc_analysis
#' @rdname poc_analysis
#' @references
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' report <- poc_analysis(list(poc_bounds(m, q), poc_exact(m, q)),
#'                        provenance = list(target = "Synthetic population"))
#' summary(report)
NULL

#' @name poc_curve
#' @rdname poc_curve
#' @references
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
NULL

#' @name poc_exact
#' @rdname poc_exact
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' poc_exact(m, q)
#' poc_exact(m, q, lp_backend = "zig")
NULL

#' @name poc_sensitivity
#' @rdname poc_sensitivity
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' poc_sensitivity(m, q, delta = c(0, .05))
NULL

#' @name poc_margin_map
#' @rdname poc_margin_map
#' @references
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' poc_margin_map(m)
NULL

#' @name poc_static_system_zig
#' @rdname poc_static_system_zig
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' head(poc_static_system_zig(2L, 2L, as.numeric(o), as.numeric(a), 10000L))
NULL

#' @name poc_confidence_bounds
#' @rdname poc_confidence_bounds
#' @references
#' Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and
#' identification. \emph{Annals of Mathematics and Artificial Intelligence},
#' 28, 287--313. \doi{10.1023/A:1018912507879}.
#' 
#' See \link{pnsbib-methods} for the method-to-function map, assumptions and
#' references. Graphics and reporting functions preserve supplied results;
#' they do not perform causal identification.
#' @examples
#' o <- rbind(reference = c(failure = .35, success = .15),
#'            active = c(failure = .20, success = .30))
#' a <- rbind(reference = c(failure = .70, success = .30),
#'            active = c(failure = .55, success = .45))
#' m <- poc_model(o, a)
#' q <- poc_query(c(reference = "failure", active = "success"))
#' map <- poc_margin_map(m)
#' samples <- setNames(lapply(1:4, function(j) as.numeric(rep(1:4, 25) == j)),
#'                     map$cell_id[map$type == "observed"])
#' region <- poc_confidence_region(m, samples, target = "Illustrative iid population",
#'   sampling_unit = "individual", assume_sampling = TRUE)
#' poc_confidence_bounds(region, q)
NULL
