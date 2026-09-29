\name{pnsbib-methods}
\alias{pnsbib-methods}
\title{Methods, Assumptions and References}
\description{Method-to-function guidance for finite probabilities of causation. Each calculation is conditional on the supplied information and structural restrictions.}
\details{
\strong{Events and complete margins.} Tian and Pearl define binary PNS, PN and PS. Li and Pearl extend finite categorical queries. Use \code{\link{poc_model}}, \code{\link{validate_poc_inputs}} and \code{\link{poc_query}} to declare the target's factual and intervention probabilities. \code{\link{poc_bounds}} implements Li--Pearl Theorems 4--11. Theorems 4 and 5 are valid analytical bounds, not generally sharp; recursive results retain their own validity label. \code{\link{poc_shu2026}} implements Shu--Wang--Li with a separately labelled repeated-factual-outcome extension and no general sharpness claim.

\strong{Exact and partial models.} \code{\link{poc_exact}} optimizes a finite nonnegative response-type distribution subject to normalization, factual consistency and the supplied margins. \code{\link{poc_constraint}}, \code{\link{poc_partial_model}} and \code{\link{poc_partial_bounds}} add selected event equalities or bands, including factual-subgroup and fixed-stratum constraints. Missing margins remain free. Charnes--Cooper transformation handles varying positive denominators in fractional queries. Exactness refers to the declared continuous finite response model, not integer finite-cohort schedules.

\strong{Longitudinal models.} \code{\link{poc_longitudinal_model}}, \code{\link{poc_longitudinal_query}}, \code{\link{poc_longitudinal_bounds}} and \code{\link{poc_curve}} apply the same finite-response construction to categorical outcome paths under named treatment regimes. Potential-path support, no anticipation and ordered absorption are explicit restrictions. These package constructions are not additional Li--Pearl theorems. More horizons alone need not tighten a bound. Identification of regime margins remains an external causal-design question, as discussed by Hernan and Robins.

\strong{Sensitivity and sampling.} \code{\link{poc_sensitivity}}, \code{\link{poc_input_sensitivity}} and \code{\link{poc_longitudinal_sensitivity}} optimize across documented hypothetical input bands. \code{\link{poc_margin_map}}, \code{\link{poc_confidence_region}} and \code{\link{poc_confidence_bounds}} instead use simultaneous regions built from independent bounded unit contributions: weighted Hoeffding or unweighted Bernoulli Clopper--Pearson intervals, with a union-bound allocation. Numerical template probabilities are ignored. Empty projected regions remain failures in coverage assessment. Bounds-in-bounds is a descriptive organization of nested feasible sets, not a new estimator or identification theorem.

\strong{Risk estimation.} \code{\link{poc_gformula}} standardizes discrete conditional risks to a common baseline distribution under declared identification assumptions. \code{\link{poc_risk_fit}} and \code{\link{poc_adjusted_risks}} implement penalized logistic prediction and target-weighted standardization, or products of discrete first-success survival probabilities. Le Cessie and van Houwelingen provide the ridge-logistic background; the package's strictly positive penalty on every coefficient, including the intercept, and its supplied penalty weights are specific implementation choices. These functions do not implement that paper's penalty-selection procedures. Hernan and Robins explain standardization and the consistency, exchangeability and positivity requirements. Predictions, identification bounds and sampling confidence intervals are distinct objects.

\strong{Numerical and reporting interfaces.} \code{\link{poc_single_zig}}, \code{\link{poc_paper_zig}} and \code{\link{poc_shu_zig}} are low-level analytical kernels. \code{\link{poc_lp_zig}} solves nonnegative continuous LPs; \code{\link{poc_static_system_zig}} and \code{\link{poc_static_event_zig}} construct static constraints and event masks. Use validated model interfaces for ordinary analysis. \code{\link{poc_analysis}}, \code{\link{poc_export_summary}}, \code{\link{poc_plot}} and \code{\link{poc_save_plot}} summarize or display existing calculations. They introduce no new statistical method and do not convert sensitivity intervals into confidence intervals. Export functions require explicit destinations.
}
\references{
Tian, J. and Pearl, J. (2000). Probabilities of causation: Bounds and identification. \emph{Annals of Mathematics and Artificial Intelligence}, 28, 287--313. \doi{10.1023/A:1018912507879}.

Li, A. and Pearl, J. (2024). Probabilities of Causation with Nonbinary Treatment and Effect. \emph{Proceedings of the AAAI Conference on Artificial Intelligence}, 38(18), 20465--20472. \doi{10.1609/aaai.v38i18.30030}.

Shu, X., Wang, S. and Li, A. (2026). Identification of Probabilities of Causation: From Recursive to Closed-Form Bounds. \emph{Proceedings of IJCAI}, 4063--4070. \doi{10.24963/ijcai.2026/452}.

Charnes, A. and Cooper, W. W. (1962). Programming with linear fractional functionals. \emph{Naval Research Logistics Quarterly}, 9(3--4), 181--186. \doi{10.1002/nav.3800090303}.

Hoeffding, W. (1963). Probability inequalities for sums of bounded random variables. \emph{Journal of the American Statistical Association}, 58(301), 13--30. \doi{10.1080/01621459.1963.10500830}.

Clopper, C. J. and Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the case of the binomial. \emph{Biometrika}, 26(4), 404--413. \doi{10.1093/biomet/26.4.404}.

le Cessie, S. and van Houwelingen, J. C. (1992). Ridge estimators in logistic regression. \emph{Applied Statistics}, 41(1), 191--201. \doi{10.2307/2347628}.

Hernan, M. A. and Robins, J. M. (2020). \emph{Causal Inference: What If}. Chapman & Hall/CRC. \url{https://miguelhernan.org/whatifbook}.
}
\seealso{\link{pnsbib-workflows}, \code{vignette("using-pnsbib", package="pnsbib")}.}
