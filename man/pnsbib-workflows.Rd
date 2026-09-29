\name{pnsbib-workflows}
\alias{pnsbib-workflows}
\title{Analyze Your Own Finite Categorical Data}
\description{
A workflow guide to the complete-margin, partial-information and longitudinal
interfaces, with an installed executable tutorial and a full export index.
Outcome construction, eligibility and causal-input identification precede the
probability calculation. Bibliometric records are one possible application.
}
\details{
For binary annual outcomes and prepared baseline design matrices,
\code{\link{poc_adjusted_risks}} estimates annual or cumulative first-success
risks standardized to one common target. \code{\link{poc_risk_fit}} fits a
single penalized logistic model. Both use Zig through rzig and require
external causal justification. The installed \code{doc/adjusted-risks.md}
and \code{doc/adjusted-risks.R} connect estimated margins to all three models.

Construct a named treatment-by-outcome factual matrix summing to one.
Complete intervention matrices have matching labels and each row sums to one.
Use \code{poc_model} for complete margins, \code{poc_partial_model} for
selected justified constraints, and \code{poc_longitudinal_model} for named
regimes and full categorical potential-path domains. Preserve observed zeros
and handle unknown measurements explicitly before aggregation.

\code{poc_query} specifies PNS, PN, PS or other finite counterfactual events.
\code{poc_bounds}, \code{poc_shu2026} and \code{poc_exact} compare methods
under complete information. Partial models use \code{poc_partial_bounds}.
Longitudinal queries include horizon PNS/PN/PS, exact trajectory pairs,
persistent, first-entry and lagged events. \code{poc_curve} preserves the
declared horizon grid; model periods require an external calendar mapping.

Sensitivity bands describe assumed input departures. Simultaneous sampling
regions instead require independent bounded unit contributions and a common
target. The confidence interface ignores all numerical template probabilities.
The discrete baseline g-formula supplies point estimates and support checks
under declared identification assumptions, without generic standard errors.

Explicit \code{lp_backend="zig"} chooses native optimization; the default
remains \code{"reference"}. There is no silent fallback. Partial/trajectory
construction and several helpers still contain R computation. Inspect
statuses, conditioning populations and sharpness labels with every interval.
Graphics use computed results and support PNG/PDF and source-table export.

The installed \file{doc/own-data.md} explains data schemas, every exported
function, subgroup/stratum inputs, categorical paths, prefix conditions,
uncertainty, backends, resource limits and error handling. Its executable
companion uses synthetic data, installs nothing and writes to a fresh temporary
directory. Use \code{run_pnsbib_examples(output="/absolute/local/new-output")}
after sourcing to retain results in an unused local directory.
}
\seealso{
\code{\link{pnsbib-package}}, \code{\link{poc_model}},
\code{\link{poc_partial_model}}, \code{\link{poc_longitudinal_model}},
\code{\link{poc_confidence_region}}, \code{\link{poc_gformula}},
\code{\link{poc_plot}}, \code{\link{poc_save_plot}}
}
\examples{
system.file("doc", "own-data.md", package = "pnsbib")
source(system.file("doc", "own-data.R", package = "pnsbib"))
pnsbib_examples$complete$PNS
unlink(pnsbib_examples$output, recursive = TRUE)
}
\references{
See \link{pnsbib-methods} for the method-to-function map, assumptions and
references. Graphics and reporting functions preserve supplied results;
they do not perform causal identification.
}
