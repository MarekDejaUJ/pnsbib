# Adjusted binary risks and probability bounds

Version 0.0.21 adds `poc_risk_fit()` and `poc_adjusted_risks()` alongside every existing probability-table function. Fitting, prediction and common-target averaging use Zig through rzig. R validates inputs and assembles models and plots. Existing bounds functions never fit a risk model implicitly.

## Run the complete example

```r
library(pnsbib)
source(system.file("doc", "adjusted-risks.R", package="pnsbib"))
pnsbib_adjusted_examples$annual$risks
pnsbib_adjusted_examples$complete
pnsbib_adjusted_examples$partial
pnsbib_adjusted_examples$curve
```

The [companion script](adjusted-risks.R) is a small synthetic example, with both groups sharing a response distribution. It fits annual and cumulative risks, sends a complete table to `poc_model()` and `poc_exact()`, sends risk bands to `poc_partial_model()`, and sends horizon margins to `poc_longitudinal_model()` and `poc_curve()`. It checks the analytical interval and fixed-centre band containment. No files are written. The [own-data guide](own-data.md) covers the unchanged categorical, subgroup, path, sensitivity, confidence and plotting interfaces.

## Prepare your data

Supply a finite numeric design matrix `x` with a first column of ones. Use only baseline covariates justified by the intervention design. Compute centring/scaling in training data and apply the same constants to `target`. The function does not build formulas, select variables, standardize columns, impute missingness or tune penalties. Training groups labelled `0` and `1` must both be represented. Labels identify reference and active regimes, not ranks.

Supply `y` as a binary numeric vector or a matrix with one annual outcome per column. For citation success, an externally chosen threshold such as at least one citation defines 0/1; ranks are unnecessary. Keep zero-outcome units unless a baseline eligibility rule independently excludes them. Selecting on later success changes the target and can introduce selection bias. Handle missingness explicitly before fitting, never by silently substituting zero.

`target` defaults to the full training design. Both groups are predicted on these same target rows. `target_weights` default to equal weights and are normalized during averaging. The factual probability table supplied to the LP must describe that same target and weighting convention. A different target population requires a corresponding observed table and transport assumptions. Group-specific observed row risks are not substitutes for common-target risks.

## Annual and cumulative modes

```r
x <- cbind(intercept=1, baseline=rep(c(-1,0,1,0,-1,0,1,0),2))
group <- rep(0:1,each=8)
y <- cbind(t1=rep(c(0,0,1,1,0,1,0,1),2),
           t2=rep(c(0,1,0,1,1,0,1,0),2))
annual <- poc_adjusted_risks(x,y,group,mode="annual")
cumulative <- poc_adjusted_risks(x,y,group,mode="cumulative")
annual$risks
cumulative$risks
single <- poc_risk_fit(x,y[,1])
single[c("coefficients","objective","iterations","gradient","backend")]
```

Annual mode fits separate penalized logistic models within each group and horizon, then averages predictions on the common target. Cumulative mode takes the same annual `y`, fits first-success hazards among units with no earlier success, converts individual hazards to `1 - product(1 - hazard)` through each horizon, then averages. It does not multiply averaged hazards. Annual success need not be absorbing. An empty positive-weight training risk set fails explicitly.

The `risks` matrix has rows `0`, `1` and one column per horizon. `predictions` has dimensions target rows by horizons by groups. In cumulative mode these are cumulative probabilities. `fits` contains coefficients, training counts and optimization diagnostics. Individual predictions do not identify a pairing of potential outcomes.

## Penalties, weights and computation

The criterion is summed weighted binary negative log likelihood plus `lambda/2 * sum(penalty * coefficients^2)`. Defaults are `lambda=1`, intercept multiplier `0.1` and other multipliers `1`. Every multiplier and `lambda` must be strictly positive. Intercept regularization keeps separated estimates finite; it is a modelling choice. Training `weights` are nonnegative frequency weights. Multiplying all training weights changes the likelihood-to-penalty balance unless `lambda` is scaled too. Changing target weights changes standardization, not coefficients.

The lower-level fitter exposes `max_iterations=100L` and `tolerance=1e-10`. Native optimization uses stable log likelihoods, Newton steps, Cholesky solves and line search. Errors and resource limits are explicit, with no R fallback. At most 256 columns and 8,000,000 cells in each training/target matrix are accepted. `poc_risk_fit()` accepts a general finite numeric design; if its first column is not an intercept, specify appropriate penalty multipliers instead of relying on the first-column default.

## Interpretation and uncertainty

Causal interpretation requires a defined intervention, consistency, conditional exchangeability given baseline covariates, overlap, a common target and an interference argument. Retrospective adoption labels require a justified connection to the intervention of interest. Out-of-fold Brier scores, log loss, calibration and support diagnostics assess prediction, not these causal assumptions. This is outcome-regression standardization, not a doubly robust estimator or automatic adjustment for treatment-confounder feedback.

Exact fitted margins yield plug-in identification bounds conditional on those values. Fitted risk plus/minus `0.05`, clipped to `[0,1]`, is an assumption band unless separately calibrated. Changing fitted centres can move or widen bounds. Nesting holds for nested information sets with a fixed centre and unchanged remaining constraints. Incompatible margins stay incompatible; inspect statuses instead of forcing an interval.

For input-stability analysis, resample whole independent units, refit the estimation procedure and preserve all horizon/regime risks together in each replicate. A percentile box around fixed factual cells is not automatically a confidence region for causal bounds. Existing sampling-region functions retain their separate contracts. The new functions provide no standard errors or confidence intervals. Cross-world dependence remains a question for the existing analytical or exact LP models.

## Binding maintenance

Builds use the bundled rzig manifest. The two native risk routines are internal; only the validated R wrappers are public exports. Package developers who regenerate bindings should run `Rscript tools/sync-risk-bindings.R .` from the package source directory afterward. The helper restores the two internal registrations. Preserve R-header discovery and platform linker configuration and rerun the tests. Ordinary package installation requires neither regeneration nor this helper.
