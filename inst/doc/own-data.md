# Using pnsbib with your own data

Version 0.0.21 also provides `poc_risk_fit()` for a native penalized binary-risk fit and `poc_adjusted_risks()` for group-specific annual or cumulative risks standardized to a common target. The [adjusted-risk guide](adjusted-risks.md) documents these additional functions and their connections to all three information models. Every workflow below remains available unchanged.

This guide covers the supported user workflows in pnsbib 0.0.21. The [executable companion](own-data.R) exercises complete, partial, multivalued, trajectory, sensitivity, sampling-region, standardization, reporting and graphics interfaces on small synthetic inputs. It installs nothing and requires no API key. Start with a defined population, treatment contrast, outcome and causal-input justification. The package calculates probability bounds after those choices.

## Start and choose a model

```r
library(pnsbib)
packageVersion("pnsbib")
source(system.file("doc", "own-data.R", package = "pnsbib"))
pnsbib_examples$complete$PNS
# lower = 0.15, upper = 0.45, status = "ok"
```

The companion writes reports and figures to a fresh R temporary directory and prints its location. Temporary output may disappear when R exits. For retained output, call `run_pnsbib_examples(output = "/absolute/local/new-output")` after sourcing, using an unused directory outside cloud storage. The known-population tutorial `nested-identification.md` explains why a true PNS of .30 can coexist with [.15,.45] bounds.

| Your justified information | Model | Calculation |
| --- | --- | --- |
| Factual table and every intervention distribution | `poc_model()` | Li–Pearl, Shu–Wang–Li or exact LP |
| Factual table, optional selected risks or bands | `poc_partial_model()` | Partial-information exact LP |
| Factual histories and paths, optional regime margins | `poc_longitudinal_model()` | Horizon/path exact LP and curves |

PNS, PN and PS name events, not competing model classes. Complete margins, partial information and trajectories describe information structures. Analytical formulas and LP describe calculation methods. Compare methods on the same population, event, information and assumptions.

## Prepare factual data and intervention information

For a static analysis, supply one row per target unit before aggregation. An identifier helps detect duplicates outside the package. Fix eligibility and the outcome threshold before examining attribution results. Keep genuine zero outcomes. Resolve missingness explicitly; `table()` otherwise silently excludes missing values. Missing treatment or outcome measurements are not zero probabilities.

```r
# Replace this synthetic data frame with read.csv("my-study.csv").
d <- data.frame(
  treatment = rep(c("reference", "reference", "active", "active"), c(35, 15, 20, 30)),
  outcome = rep(c("failure", "success", "failure", "success"), c(35, 15, 20, 30)))
d$treatment <- factor(d$treatment, levels = c("reference", "active"))
d$outcome <- factor(d$outcome, levels = c("failure", "success"))
stopifnot(!anyNA(d), nrow(d) > 0L)
o <- unclass(table(d$treatment, d$outcome)) / nrow(d)
dimnames(o) <- list(levels(d$treatment), levels(d$outcome))

# Declared hypothetical intervention probabilities for this example.
a <- rbind(reference = c(failure = .70, success = .30),
           active = c(failure = .55, success = .45))
m <- poc_model(o, a, lp_backend = "zig")
```

The whole factual matrix `o` must sum to one. Each row of intervention matrix `a` must sum to one. Both matrices need identical row and column names in identical order, at least two states per dimension, finite entries and every cell including zeros. A within-treatment observed proportion estimates `P(Y | X)`, which becomes `P(Y_x)` only under a justified identification argument. A treated-group effect alone does not fill a complete population table. With unknown intervention risks, omit them through the partial model instead of inserting guesses or zero rows.

`validate_poc_inputs(o, a, lp_backend = "zig")` performs the same validation as `poc_model()`. The returned `checked_lp` records whether full feasibility was checked within `max_variables`; bypassing or exceeding that check is not a feasibility certificate. A complete-margin long-table alternative has exactly one row per treatment/outcome pair and columns `treatment`, `outcome`, `observed_joint`, `interventional`; call `poc_model(long_table, lp_backend = "zig")` without a second matrix.

Continuous measurements require externally defined finite categories. The package neither ranks outcomes nor chooses cut points. Labels such as `low`, `middle`, `high` are nominal unless an explicit temporal ordering is supplied. For overlapping database indicators, encode mutually exclusive joint states such as `00`, `01`, `10`, `11`, not overlapping database names as if they were disjoint treatments. The package accepts general categorical data beyond bibliometrics.

## Static queries, partial information and sensitivity

```r
q <- list(
  PNS = poc_query(c(reference = "failure", active = "success")),
  PN = poc_query(c(reference = "failure"),
                 observed_x = "active", observed_y = "success", conditional = TRUE),
  PS = poc_query(c(active = "success"),
                 observed_x = "reference", observed_y = "failure", conditional = TRUE))
poc_bounds(m, q$PNS)       # Printed Li–Pearl bound, native by default
poc_shu2026(m, q$PNS)     # Separate closed-form comparator
poc_exact(m, q$PNS)       # Exact LP, inherits the explicit Zig backend
poc_exact(m, q$PN)
poc_exact(m, q$PS)
```

The named counterfactual vector can contain several treatment requirements. Adding factual `observed_x` and/or `observed_y` defines a joint event; `conditional = TRUE` divides by its factual condition. The representation covers the implemented preservation, replacement, substitute and multi-hypothetical Li–Pearl query family. A three-treatment example is `poc_query(c(off = "low", on = "high", alternative = "middle"))` on conformable three-treatment matrices. A static query names one outcome per queried treatment; the longitudinal interface additionally supports outcome sets.

```r
free <- poc_partial_model(o, lp_backend = "zig")
poc_partial_bounds(free, q$PNS)
constraints <- list(
  poc_constraint(poc_query(c(reference = "success")), .25, .35,
                 label = "reference_population_risk"),
  poc_constraint(poc_query(c(active = "success")), .40, .50,
                 label = "active_population_risk"))
partial <- poc_partial_model(o, constraints, lp_backend = "zig")
fit <- poc_partial_bounds(partial, q$PNS, return_witness = TRUE)
fit
attr(fit, "witnesses")
```

`upper` defaults to `lower` in `poc_constraint()`, giving an equality. A conditional constraint such as `poc_query(c(reference = "success"), observed_x = "active", conditional = TRUE)` concerns the active factual subgroup. It must not be relabeled as a whole-population intervention risk. Every supplied constraint is tested jointly for compatibility.

For disjoint baseline strata, supply a named list of matrices containing `P(S=s,X=x,Y=y)`. Their total across all strata is one. For example, `list(A = .4 * o, B = .6 * o)` fixes stratum weights .4 and .6. A constraint with `stratum = "A"` is conditional on that stratum; `poc_partial_bounds(model, query, stratum = "A")` also targets that stratum. Independently row-normalizing or stratum-normalizing the factual matrices changes the input contract.

```r
poc_sensitivity(m, q$PNS, delta = c(0, .025, .05, .1))
poc_input_sensitivity(m, q$PN, observed_delta = .02,
  interventional_delta = .05, minimum_condition_probability = .1)
```

`poc_sensitivity()` varies intervention cells and holds factual cells fixed. `poc_input_sensitivity()` permits both factual and intervention cell variation, preserving normalization and consistency. Conditional queries optimize numerator and denominator together. A conditioning floor is an additional assumption and must be disclosed. These absolute-probability bands are neither confidence intervals, exposure misclassification rates nor automatically calibrated confounding parameters. Bounds-in-bounds describes joint optimization over compatible inputs and response distributions; the package has no universal arbitrary joint-region or DiD estimator.

## Longitudinal histories, horizons and categorical paths

Prepare a factual history/path table with `history`, `path`, `prob`, including all declared cells and zeros. Its probabilities sum to one over the full table. The named `regimes` matrix has one row per fixed treatment history and one column per period. Factual history labels must match regime names. Period numbers are positions in that matrix; retain your calendar mapping separately. The potential path domain includes possible unobserved paths, not just paths seen in the data.

```r
regimes <- rbind(off = c(0, 0), on = c(1, 1))
paths <- expand.grid(history = c("off", "on"), path = c("00", "01", "11"),
                     stringsAsFactors = FALSE)
paths$prob <- c(.2, .1, .1, .1, .2, .3)
tm <- expand.grid(regime = c("off", "on"), horizon = 1:2,
                  stringsAsFactors = FALSE)
tm$prob1 <- c(.4, .6, .6, .8)  # Explicit hypothetical intervention risks
lm <- poc_longitudinal_model(paths, regimes, time_margins = tm,
                             absorbing = TRUE, lp_backend = "zig")
lq <- poc_longitudinal_query("pns", "on", "off", horizon = 2)
poc_longitudinal_bounds(lm, lq)
curve <- poc_curve(lm, "on", "off", horizons = 1:2)
curve
poc_curve(lm, "on", "off", bands = list(time_delta = .05))
```

Omit `time_margins` for unrestricted intervention risks. Binary horizon margins use `regime`, `horizon`, `prob1`; categorical horizon margins use `regime`, `horizon`, `outcome`, `prob`. Do not mix those schemas. Intervention-path constraints use `regime`, `path`, `prob` and can specify selected cells; a complete supplied path distribution sums to one within each regime. Unspecified path, state and horizon information remains free subject to joint consistency. The model must be feasible before sensitivity bands are applied.

| Query kind | Event at horizon `t` |
| --- | --- |
| `pns` | Active-set response under the active regime and reference-set response under the reference regime |
| `pn` | Reference-set potential response among factually active-set outcomes under the active history |
| `ps` | Active-set potential response among factually reference-set outcomes under the reference history |
| `trajectory` | The specified pair of complete path identifiers |
| `persistent` | Both required responses at every period from `start` through `t` |
| `event_time` | First active-set entry at `t`, and reference-set response at `t` |
| `lagged` | Active-set response at `t`, reference-set response at `t - lag` |

Exact trajectory queries require `regime_path` and `reference_path`. A lag must be smaller than the horizon. `poc_curve()` handles all listed kinds except exact trajectory pairs, preserving the requested horizon/kind order. For reverting outcomes, an event-time reference-set condition at `t` does not imply absence of earlier reference success.

PN/PS default to `condition_on = "history"`. Use `"prefix"` only when the question conditions on all factual histories sharing the treatment prefix through the horizon. Later continuations are pooled. With `no_anticipation = TRUE`, potential outcomes agree for regimes sharing full treatment prefixes; a thresholded outcome set is not the basis of that restriction. The companion supplies a four-regime prefix example.

For categorical paths, supply a named matrix `outcome_paths` with a full state in each period, and use its row identifiers in the factual `path` column. Set `outcome_order` explicitly when `absorbing = TRUE`. The companion includes six ordered two-period paths over `low`, `middle`, `high` and queries `regime_outcome = c("middle", "high")`, `reference_outcome = "low"`. Absorption constrains outcomes within a regime over time; it does not imply treatment monotonicity or monotone PNS curves.

Use `poc_longitudinal_sensitivity()` with `observed_delta`, `path_delta` and `time_delta` for the corresponding supplied cells. Missing intervention margins remain free. The optional `minimum_condition_probability` applies to PN/PS and changes their admissible region. `poc_curve(..., bands = list(...))` applies these same controls to the requested query grid.

## Sampling uncertainty and optional input estimation

Confidence projection has a different input contract. `poc_margin_map(template)` gives cell identifiers and event labels for a complete static or longitudinal structural template. Partial-model templates are unsupported. Supply a named list of vectors keyed by those identifiers, with one bounded contribution per independent sampling unit. The full margin map must be retained; identifiers belong to that template.

```r
map <- poc_margin_map(m)
samples <- list()
for (i in which(map$type == "observed")) {
  samples[[map$cell_id[i]]] <- as.numeric(
    d$treatment == map$treatment[i] & d$outcome == map$outcome[i])
}
region <- poc_confidence_region(m, samples,
  target = "Illustrative iid target", sampling_unit = "individual",
  method = "clopper_pearson", conf_level = .95, assume_sampling = TRUE)
poc_confidence_bounds(region, q$PNS, witnesses = TRUE, lp_backend = "zig")
```

This example supplies only factual-cell observations. The template's intervention values are ignored, so they cannot tighten the confidence calculation. All numerical template probabilities are ignored, including cells lacking samples; unsupplied cells and `numeric(0)` remain unrestricted. Separate randomized or otherwise justified risk samples can supply intervention cells for the same target.

`method = "clopper_pearson"` requires unweighted iid binary individual contributions. `method = "hoeffding"` accepts independent bounded contributions, optionally with externally fixed nonnegative weights. For clusters, use `sampling_unit = "cluster"` and one bounded contribution per independent cluster. Copying a cluster mean onto every member creates no extra independent information. Equal-cluster and equal-individual weighting can define different targets. The simultaneous guarantee depends on the declared sampling and structural assumptions; empty regions must remain in coverage accounting.

The optional `poc_gformula()` helper standardizes discrete baseline-covariate cells to a common empirical target. It returns `model`, `support` and `diagnostics`. Every treatment must have positive weighted support in every target stratum. Even numeric covariates are treated as discrete categories. Set factor levels in advance to preserve unobserved treatment states. `assume_identification = TRUE` declares externally justified consistency, conditional exchangeability, positivity and interference assumptions; the function cannot establish them. It supplies point estimates, not generic g-formula standard errors or LongBet estimation. The helper is R-based; downstream exact bounds can explicitly select Zig with `poc_exact(fit$model, query, lp_backend = "zig")`.

## Inspect, plot and save results

Inspect `status`, `method`, `sharpness`, event labels, assumptions and conditioning information alongside `lower` and `upper`. A zero lower endpoint is a valid result permitting no attribution within the model. `undefined_condition` means the factual conditioning event has zero probability; keep the undefined result instead of replacing it with zero. Constructors reject incompatible inputs. LP failures or resource-limit errors require diagnosis, never deletion of inconvenient rows.

The printed Li–Pearl Theorems 4/5 retain `valid_analytical`; their expressions need not be sharp with complete inputs. Recursive bounds retain `valid_recursive`. Shu–Wang–Li results retain `valid_closed_form`, including its labeled repeated-factual-outcome extension. Exact bounds use `sharp_given_inputs` where applicable: sharpness concerns the specified feasible continuous response model. It is distinct from valid causal input identification, integer fixed-cohort schedules or an optimal confidence procedure.

```r
report <- poc_analysis(list(poc_bounds(m, q$PNS), poc_exact(m, q$PNS)),
  provenance = list(target = "Synthetic example", risks = "Declared assumptions"))
summary(report)
out <- tempfile("pnsbib-report-")
dir.create(out)
poc_export_summary(report, file.path(out, "tables"), format = "tsv")
fig <- poc_plot(report, draw = FALSE, monochrome = TRUE)
poc_save_plot(fig, file.path(out, "static"), formats = c("png", "pdf"), dpi = 600)
curves <- poc_plot(curve, type = "curve", draw = FALSE, monochrome = TRUE)
poc_save_plot(curves, file.path(out, "horizons"), formats = c("png", "pdf"))
```

`plot(fig)` draws an existing plot object. Plotting performs no estimation or optimization. Curve displays preserve undefined gaps and decreases, require unique horizons within each series and allow at most eight series. `poc_save_plot()` also writes source-data, display-mapping and provenance tables and an RDS object. Its filename is a prefix without an image extension; the parent directory must exist. Export refuses collisions. Width is 5–30 inches, height 3–60, PNG resolution 72–1200 dpi, with a 100-million-pixel cap. PDF remains vector. Summary export accepts TSV or CSV and retains all supplied report tables.

## Function reference and operational limits

The following index covers every exported function. The six generated native wrappers are developer interfaces; use high-level model/query functions for normal analysis.

| Export | Purpose |
| --- | --- |
| `poc_model`, `validate_poc_inputs` | Complete static input construction and validation |
| `poc_query` | Static counterfactual and factual conditions |
| `poc_bounds` | Printed Li–Pearl Theorems 4–11; `method="exact"` selects LP |
| `poc_shu2026` | Named closed-form comparator |
| `poc_exact` | Complete static exact LP |
| `poc_constraint` | Event equality or probability band |
| `poc_partial_model`, `poc_partial_bounds` | Partial, subgroup and disjoint-stratum analysis |
| `poc_sensitivity`, `poc_input_sensitivity` | Intervention-only or joint factual/intervention bands |
| `poc_gformula` | Discrete baseline standardization and support diagnostics |
| `poc_longitudinal_model` | Regimes, outcome-path domain and supplied margins |
| `poc_longitudinal_query`, `poc_longitudinal_bounds` | Seven longitudinal query kinds and exact LP |
| `poc_longitudinal_sensitivity` | Factual-path, intervention-path and horizon bands |
| `poc_curve` | Ordered horizon/query grids |
| `poc_margin_map`, `poc_confidence_region`, `poc_confidence_bounds` | Simultaneous sampled-margin region and projection |
| `poc_analysis`, `poc_export_summary` | Provenance-aware reports and table export |
| `poc_plot`, `poc_save_plot` | Interval/curve graphics and PNG/PDF export |
| `poc_single_zig`, `poc_paper_zig`, `poc_shu_zig` | Low-level native analytical kernels |
| `poc_lp_zig` | Low-level continuous LP and numerical certificates |
| `poc_static_system_zig`, `poc_static_event_zig` | Low-level complete static system and event construction |

Use `help("function_name", package = "pnsbib")` for exact arguments and return fields. `help("pnsbib-workflows", package = "pnsbib")` provides the workflow entry point. Standard `summary()` methods handle results and analysis collections; `print()` and `plot()` methods handle presentation objects.

Analytical calculations default to Zig; `use_zig = FALSE` selects R oracles. Continuous LP defaults remain `lp_backend = "reference"`; explicitly choose `"zig"` in constructors or override per calculation. Saved models lacking a backend field retain reference behavior. There is no silent fallback. Complete static response construction and LP can be native, while partial/trajectory construction, additional band/fractional helpers and g-formula computation still include R. Full native porting is not claimed.

For `m` treatments and `k` outcomes, complete static LP size is `m * k^m`; disjoint strata multiply the corresponding partial-model size. A longitudinal response assignment chooses one potential path per regime, so path-domain and regime counts can cause exponential growth. Default `max_variables` is 10,000; increasing the cap does not reduce memory needs. The analytical recursion has a 1,048,576-state memo limit, the closed form a 1,048,576-input-cell limit, and the native LP a 16,777,216-cell tableau/scratch budget and normally 200,000 pivots. Use supported scientifically justified path restrictions, or smaller declared questions, before simply raising limits. Do not remove possible paths merely because none were observed.

Common errors have specific remedies: a normalization error requires checking joint versus conditional probabilities; a consistency error requires reconciling information about the same target; a positivity error requires an honest redesign or supported partial information; a path error requires matching regime/path labels and full domain support; an export collision requires a new prefix. Preserve the failed inputs and status. A solver's successful completion establishes numerical compatibility, not the study's causal assumptions.
