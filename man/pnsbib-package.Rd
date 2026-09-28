\name{pnsbib-package}
\alias{pnsbib-package}
\alias{pnsbib}
\docType{package}
\title{Finite Probabilities of Causation with Partial Information}
\description{
Bounds on probabilities of necessity and sufficiency, necessity and sufficiency
for finite categorical treatments and outcomes. Bibliometric coverage is one
application; the probability-table interface is domain-independent.
}
\details{
\code{\link{poc_adjusted_risks}} and \code{\link{poc_risk_fit}} provide
opt-in Zig-backed binary-risk estimation from prepared baseline designs.
They supplement the probability-table interfaces and do not automatically
establish causal identification or sampling confidence intervals.

Choose an information model: \code{poc_model} for complete intervention
margins, \code{poc_partial_model} for selected event constraints and
\code{poc_longitudinal_model} for categorical regimes and outcome paths.
PNS, PN and PS are estimands, while Li-Pearl, Shu-Wang-Li and exact LP are
calculation methods. Exactness is relative to the specified feasible model.

Outcomes and eligibility are defined outside the package. Labels are not
automatically ranked. Retain genuine zero outcomes and distinguish unknown
measurements from zeros. Post-outcome selection changes the population.

Bounds-in-bounds describes optimization over a declared input region followed
by compatible response distributions, or equivalent joint optimization.
Assumption bands and confidence regions have different interpretations.
The installed tutorial explains nesting, coupled input restrictions,
conditional denominators and a known generating population. No general
DiD estimator or arbitrary joint-region adapter is provided.

The explicit \code{lp_backend="zig"} selects native optimization through rzig.
Complete static construction is also native; partial and trajectory
construction and additional numerical helpers still contain R computation.
The default remains \code{"reference"}. There is no silent fallback.
}
\seealso{
\code{\link{pnsbib-workflows}},
\code{\link{poc_model}}, \code{\link{poc_partial_model}},
\code{\link{poc_longitudinal_model}}, \code{\link{poc_confidence_bounds}},
\code{\link{poc_input_sensitivity}}, \code{\link{poc_plot}}
}
\examples{
system.file("doc", "own-data.md", package = "pnsbib")
system.file("doc", "nested-identification.md", package = "pnsbib")
source(system.file("doc", "nested-identification.R", package = "pnsbib"))
}
