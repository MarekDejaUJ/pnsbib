# Reproduce probability bounds and numerical experiments

These R scripts reproduce 570 adjusted institutional bounds from aggregate risks and path counts, verify native LP certificates and run 600 known-population repetitions. They need pnsbib 0.0.22 and the suggested `digest` package. No Python, API credentials or individual researcher records are used.

## Run from R

```r
library(pnsbib)
reproduction <- system.file("reproduction", package="pnsbib")
output <- tempfile("pnsbib-reproduction-")
dir.create(output)
previous <- setwd(reproduction)
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, c("bounds.R", shQuote(file.path(output,"bounds")), "data")) == 0)
stopifnot(system2(rscript, c("audit-lp.R", shQuote(file.path(output,"bounds")),
                           shQuote(file.path(output,"audit")))) == 0)
stopifnot(system2(rscript, c("plot.R", shQuote(file.path(output,"bounds","bounds.tsv")),
                           shQuote(file.path(output,"figure")))) == 0)
stopifnot(system2(rscript, c("simulation.R", shQuote(file.path(output,"simulation")))) == 0)
setwd(previous)
output
```

Alternatively, from the source `inst/reproduction` directory, run `Rscript bounds.R /absolute/new-output data`, then `Rscript audit-lp.R /absolute/new-output /absolute/new-audit`. The simulation command is `Rscript simulation.R /absolute/new-simulation`. Every output directory must be new, with an existing writable parent. The bounds calculation may take several minutes and writes intermediate optimization certificates to the chosen output directory.

## Data and interpretation

The institutional target contains 873 researchers, including 303 switchers and 570 controls. The opportunity origin is 2021 and annual/cumulative outcome horizons are 2021–2025. The aggregate data contain twenty common-target estimated risks, the fitted-input stability radius and joint group/path counts. The original researcher-level panel is not included. Re-estimating its risks requires that panel; these scripts reproduce attribution conditional on the supplied fitted inputs.

The synthetic design uses seed 2026092803 and 200 repetitions for each of three mechanisms: correctly specified outcome regression, omitted quadratic outcome structure and poor treatment overlap. Exact and five-point risk-band PNS, PN and PS intervals are checked against known generating values. Misspecification failures are retained. Truth containment in an identification interval is distinct from sampling-confidence coverage.

The main outputs are `bounds.tsv`, `checks.tsv`, `endpoint-witnesses.tsv`, `lp-certificates.tsv`, `simulation-risks.tsv`, `simulation-bounds.tsv` and `simulation-checks.tsv`. `audit-lp.R` checks feasibility, dual feasibility and optimality gaps from the saved native programs without calling a solver. Sampling stability and assumed risk bands retain separate interpretations. No estimator or solver silently substitutes for a failed calculation.

`plot.R` produces a six-panel PNG/PDF figure: PNS, PN and PS in rows, annual and cumulative outcomes in columns. It plots partial static LP intervals for exact adjusted risks and five-point risk bands, checks their agreement with the matched trajectory LP, and checks complete analytical/exact agreement at exact risks. The plotted intervals are identification bounds, not confidence intervals. `plot-data.tsv` preserves every plotted endpoint and `plot-checks.tsv` records the comparisons. A grayscale preview is included. Run it separately with `Rscript plot.R /absolute/bounds/bounds.tsv /absolute/new-figure`.
