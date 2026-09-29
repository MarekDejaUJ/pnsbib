# pnsbib: probabilities of causation

`pnsbib` bounds probabilities of necessity and sufficiency (PNS), necessity (PN), and sufficiency (PS) for **finite categorical treatments and outcomes**. Bibliometric indexing is one application. The inputs are named probability tables and declared constraints, not DOIs, citation counts or provider records.

The package includes Zig-backed binary-risk estimation and common-target standardization alongside the complete, partial and longitudinal bounds interfaces. Probability tables remain the inputs to the bounds engines; risk-estimation functions help estimate those tables from externally prepared binary outcomes and baseline covariates. See [NEWS.md](NEWS.md) for version history.

## Install from GitHub

Install R and its platform-specific build tools, then run:

```r
install.packages("remotes")
remotes::install_github("MarekDejaUJ/pnsbib", upgrade = "never")
```

The repository root is the R package; no `subdir` argument is needed. Source installation automatically obtains Zig 0.16.0 when that version is unavailable. The R-only installer downloads the official archive, checks its pinned size and SHA-256 before extraction, and keeps the compiler in the build R session's temporary directory, outside the package source. Temporary compiler files are removed after compilation. A first download is approximately 49–93 MiB depending on platform; extraction needs additional disk space. An existing compatible compiler is reused; otherwise each source installation obtains a temporary copy. There is no Python requirement, manual Zig download or special project folder to configure. An installed binary requires no compiler.

The generated rzig framework is bundled; `rzig` is a suggested development dependency for regenerating bindings, not a runtime requirement. Automatic compiler setup is provided by pnsbib, independently of rzig 0.2.3. Imported R dependencies are installed automatically. Building requires a C toolchain suitable for R: Xcode Command Line Tools on macOS, R development headers and a C compiler on Linux, or the matching Rtools distribution on Windows. Automatic archive selection covers aarch64 and x86_64 on macOS, Linux and Windows. The current release is locally validated on macOS; other platforms require their own checks.

For managed or offline installations, `ZIG` can select an existing Zig 0.16.0 executable, `PNSBIB_ZIG_CACHE` can explicitly opt into a persistent cache at a user-chosen absolute path, and `PNSBIB_ZIG_ARCHIVE` can supply the matching official archive locally. The local archive undergoes the same pinned size and SHA-256 checks. `PNSBIB_ZIG_DOWNLOAD=false` disables network downloads. A wrong explicit compiler, corrupt archive or invalid cache entry causes an error; the installer does not silently substitute an unverified compiler. These settings are optional and are unnecessary for ordinary online installation.

## Installed documentation

After installation, `vignette("using-pnsbib", package = "pnsbib")` provides a runnable, installation-free workflow. `help("pnsbib-methods")` maps every public function to its method and references; `citation("pnsbib")` gives the software citation. Function help pages include small executable examples. The installed own-data and adjusted-risk guides provide further cases.

## Estimate risks, then bound attribution

Use `poc_adjusted_risks()` for group-specific annual risks or cumulative first-success risks standardized to one declared target. Use `poc_risk_fit()` for a single penalized logistic fit. Both compute in Zig through rzig. The [adjusted-risk guide](inst/doc/adjusted-risks.md) and [runnable example](inst/doc/adjusted-risks.R) connect estimated risks to complete, partial and longitudinal models, with explicit assumption bands. These functions do not replace the existing LP or infer causal assumptions from predictive performance.

```r
library(pnsbib)
source(system.file("doc", "adjusted-risks.R", package = "pnsbib"))
pnsbib_adjusted_examples$complete
pnsbib_adjusted_examples$partial
pnsbib_adjusted_examples$curve
```

## Start with your own data

The [complete own-data guide](inst/doc/own-data.md) covers all supported functions, data schemas, causal-input requirements, examples and troubleshooting. The [runnable companion](inst/doc/own-data.R) uses small synthetic data and exercises complete, partial, multivalued and longitudinal models, subgroup/stratum constraints, all seven path query kinds, prefix conditioning, sensitivity, sampling regions, discrete standardization, reports and PNG/PDF exports.

```r
library(pnsbib)
help("pnsbib-workflows", package = "pnsbib")
source(system.file("doc", "own-data.R", package = "pnsbib"))
pnsbib_examples$complete$PNS
# [0.15, 0.45] under the supplied hypothetical complete risks
```

The script writes only to a new temporary directory and prints its location. To retain its outputs, use `run_pnsbib_examples(output = "/absolute/local/new-output")` after sourcing. No acquisition, credentials or installation is performed by the examples. Replace synthetic rows with your data only after declaring the population, treatment, outcome, time origin, missingness rule and source of intervention information. The guide explains why observed row risks are not automatically intervention risks.

## Choose an information model, then a method

| Information model | Constructor | Main calculation |
|---|---|---|
| Complete-margin static | `poc_model()`: joint factual distribution and every intervention row | `poc_bounds()` (Li-Pearl), `poc_shu2026()`, `poc_exact()` |
| Partial-information static | `poc_partial_model()`: fixed factual cells and only justified additional equalities/bands | `poc_partial_bounds()` |
| Longitudinal trajectories | `poc_longitudinal_model()`: factual history/path cells, categorical regimes, optional justified path/time margins | `poc_longitudinal_bounds()`, `poc_curve()` |

PNS, PN and PS are estimands. Li-Pearl, Shu-Wang-Li and exact LP are computational methods. Complete, partial and longitudinal are information models. These distinctions prevent an analytical/LP comparison from silently changing its scientific question.

The observational table sums to one over all cells. Each complete intervention row sums to one. Unspecified intervention margins in partial or longitudinal models remain unknown. Observational group risks do not automatically identify intervention risks. Factual-subgroup constraints refer to that subgroup, not the whole population. Stratum tables contain joint masses across an exhaustive disjoint population partition.

## A general-purpose example with known truth

Consider a synthetic population, half receiving each policy. Conditional response-type probabilities for `(Y_reference,Y_active)=(0,0),(0,1),(1,0),(1,1)` are `(.5,.2,.2,.1)` in the reference group and `(.3,.4,.1,.2)` in the active group. Assignment depends on response type. The resulting population has PNS `.30`, PN `2/3` and PS `2/7`.

```r
library(pnsbib)
observed <- rbind(reference = c(failure = .35, success = .15),
                 active = c(failure = .20, success = .30))
intervention <- rbind(reference = c(failure = .70, success = .30),
                     active = c(failure = .55, success = .45))
model <- poc_model(observed, intervention, lp_backend = "zig")
query <- poc_query(c(reference = "failure", active = "success"))
poc_bounds(model, query)
poc_shu2026(model, query)
poc_exact(model, query)
# All three give PNS [.15, .45] for these identical inputs.
```

The generating distribution supplies both intervention rows. They are not estimated by normalizing the observed rows. Matching bounds in this example do not establish general sharpness of either analytical family. The [installed tutorial](inst/doc/nested-identification.md) reconstructs the population and explains missing-risk bands; its [executable companion](inst/doc/nested-identification.R) is checked in the test suite.

```r
source(system.file("doc", "nested-identification.R", package = "pnsbib"))
```

## Outcomes, zeros and population selection

Construct outcomes outside the package, using categories and follow-up windows fixed by the scientific question. A binary success indicator is sufficient for binary PNS; ordinal categories allow multivalued events. Labels are categorical and are never automatically ranked, normalized to percentiles or interpreted as citation counts. Continuous-outcome identification is outside this finite API.

Retain zero outcomes when they belong to the target population. Filtering on a post-treatment outcome, such as retaining only cited papers, changes the population and can create selection bias. In the example, retaining factual successes changes true PNS from `.30` to `4/9` and removes the conditioning group needed for PS. A failed acquisition is unknown, not zero. Prespecified baseline eligibility can define a different target, but its input probabilities must all refer to that same target.

## Bounds-in-bounds: nested information, not a new estimator

An inner interval optimizes a causal query for fixed input probabilities and structural assumptions. An outer interval optimizes over a declared region of admissible inputs. Shrinking that region cannot widen the outer interval for a fixed query and model. This descriptive "bounds-in-bounds" view unifies existing sensitivity and confidence-projection interfaces; it is not a new identification theorem or a replacement for a study design.

Use `poc_constraint()` with `poc_partial_model()` for justified event equalities/bands, `poc_input_sensitivity()` or `poc_longitudinal_sensitivity()` for their documented cell-band models, and `poc_confidence_region()` with `poc_confidence_bounds()` for supported sampling regions. The package does not expose a universal arbitrary joint-region adapter or a DiD estimator.

Input dependence matters. Equal binary intervention risks permit PNS at most `.5`. Replacing their joint equality by two unrestricted marginal intervals allows PNS up to `1`; checking only the two equal-risk corners incorrectly gives zero. Conditional ratios require one compatible numerator and denominator, with undefined zero conditions retained. A rectangular relaxation can supply a conservative outer interval, not the sharp interval for a discarded joint restriction.

Assumption bands, sampling confidence regions, measurement/source envelopes and posterior regions have different interpretations. A selection-sensitivity parameter names one missing-risk restriction; "selection" is not a name for all probability bounds.

## Analytical validity and exactness

The full printed Li-Pearl Theorems 4-11 query family is implemented. Theorems 4/5 retain `valid_analytical` labels because complete-information counterexamples refute general sharpness of those printed expressions. Recursive bounds retain `valid_recursive`. The separate Shu-Wang-Li comparator retains `valid_closed_form`; its repeated factual-outcome extension is explicitly labelled. Numerical agreement does not prove general sharpness.

Exact LP endpoints are sharp for the declared feasible continuous finite-response model, within numerical tolerance. They do not certify historical labels, causal identification or deterministic finite-cohort integer support. Incompatible and undefined cases remain explicit. Fixed-cohort integer schedules require a separate integer analysis when additional restrictions make discreteness consequential.

## Longitudinal outcomes and uncertainty

Declare the full potential path domain, not just paths observed in the sample. Named `outcome_paths` support categorical paths. Absorption requires a scientifically justified ordering; otherwise paths may revert. No anticipation compares full states along shared treatment prefixes, not merely thresholded query sets. Horizon positions are model periods, not automatically time since indexing. Absorbing outcomes need not produce monotone PNS curves.

`poc_curve()` preserves ordered horizon/query combinations and undefined statuses. Event-time queries describe first entry into the active outcome set with reference-set membership at that horizon. For reverting paths this does not assert that the reference never entered the active set earlier.

Confidence projection uses its model as a structural template only: the template's numerical probabilities are ignored and unsupplied sampled cells remain free. Hoeffding-Bonferroni regions require independent bounded contributions and externally fixed weights. Clopper-Pearson-Bonferroni regions require unweighted iid Bernoulli contributions. Within-unit dependence across cells is allowed. Empty regions remain in coverage accounting.

Supply one contribution per genuinely independent sampling unit. Repeating a cluster mean over its papers creates no additional independent information. Estimated weights, generic g-formula standard errors and dependent-cluster inference need a separate justification. The finite discrete g-formula supplies regime margins only under its stated identification and support assumptions; it does not supply a generic sampling inference procedure.

## Native computation and reference oracles

| Route | Current native scope |
|---|---|
| Li-Pearl analytical family | Zig through rzig; `use_zig=FALSE` selects the R oracle |
| Shu-Wang-Li comparator | Zig through rzig; `use_zig=FALSE` selects the R oracle |
| Complete static exact | Native response construction, event masks and LP with `lp_backend="zig"` |
| Partial/trajectory LP and curves | Native optimization available; response construction remains R |
| Sensitivity/confidence routes | Native optimization and static base construction available; additional band/fractional construction remains R |
| Other numerical helpers | See individual help; the full native port remains incomplete |
| Public validation, reporting, graphics | R |

The explicit model `lp_backend` is inherited by LP calls and can be overridden per call. The default remains `"reference"`; saved models lacking that field retain reference behavior. There is no silent fallback. Resource limits fail explicitly; language choice supplies no automatic speed claim. Production numerical development targets Zig through rzig, with R mathematical oracles retained for validation.

## Plot and export

```r
comparison <- poc_analysis(list(poc_bounds(model, query), poc_exact(model, query)),
  provenance = list(inputs = "known synthetic population"))
figure <- poc_plot(comparison, draw = FALSE, monochrome = TRUE)
plot(figure)
# Existing local directory, unused prefix:
poc_save_plot(figure, tempfile("pnsbib-comparison-"),
              formats = c("png", "pdf"), dpi = 600)
```

Graphics consume computed results and do not invoke solvers. Export preserves endpoints, status labels, source tables, display mapping and provenance. Longitudinal plots preserve decreases and gaps. PNG supports raster delivery and PDF preserves vectors. Quartz-backed PNG on supported Macs does not require XQuartz.

## Installation and checks

Build from source in a local working directory. Install the R dependencies and suggested test/development packages first; compiler setup follows the automatic process above. No bibliographic API keys are required.

```r
install.packages(c("lpSolve", "digest", "testthat", "rzig", "knitr", "rmarkdown"))
```

```sh
R CMD build pnsbib
R CMD check --no-manual pnsbib_0.2.0.tar.gz
R CMD INSTALL pnsbib_0.2.0.tar.gz
```

To obtain the source, run `git clone https://github.com/MarekDejaUJ/pnsbib.git` from the parent directory before these commands. `R CMD check` runs the R regression suite and installed examples. To run the standalone kernel tests from the cloned `pnsbib` directory, use the same managed build helper so that the temporary compiler remains available until testing finishes:

```sh
(cd src/rzig && Rscript --vanilla ../../inst/toolchain/zig-toolchain.R --build test --release=safe -j2 --summary all --cache-dir .zig-cache --global-cache-dir .zig-global-cache)
```

## Reproduce analyses

The installed examples above reproduce small known-population calculations and all three information models. The [reproduction scripts](inst/reproduction/README.md) also reproduce 570 adjusted institutional bounds from supplied aggregate risks and path counts, verify their native LP certificates, and run 600 seeded numerical experiments. They require no API key or private panel. Re-estimating the institutional risks themselves requires the underlying researcher panel, which is not distributed.

## License

GPL-3. See [COPYING](COPYING) for the full license and [copyright notices](inst/COPYRIGHTS) for the bundled rzig framework.
