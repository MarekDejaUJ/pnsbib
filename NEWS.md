# pnsbib 0.0.22

Package-only distribution under GPL-3, with retained notices for the bundled MIT-licensed rzig framework. Added automatic R-only, checksum-verified Zig 0.16.0 setup for source installation, GitHub installation instructions, portable reproduction scripts, six-panel adjusted attribution graphics and aggregate inputs. Numerical functions, signatures, tolerances and defaults are unchanged from 0.0.21.

# pnsbib 0.0.21

Additive native risk-estimation release. New `poc_risk_fit()` and `poc_adjusted_risks()` implement strictly penalized binary logistic fitting and common-target annual or cumulative first-success risks in Zig through rzig. Existing exports, LP code, analytical formulas, tolerances and defaults remain. Added installed guide, runnable complete/partial/longitudinal example and regression tests. Estimated margins require external causal justification; no automatic confidence interval or cross-world identification is supplied.

# pnsbib 0.0.20

Own-data documentation release. Added an installed workflow guide, runnable synthetic examples, workflow help and a complete exported-function index. Examples cover all public information models, categorical/path/prefix queries, subgroup/stratum constraints, sensitivity, confidence projection, baseline standardization and plot/table exports. Numerical source, signatures, tolerances and backend defaults are unchanged from 0.0.19. Empirical supplements retain their original package pins.

# pnsbib 0.0.19

Native LP stability fix for a five-period annual sufficiency sensitivity program. Phase I stops when both artificial mass and its objective reach zero within the existing feasibility tolerance. A safeguarded two-pass ratio test avoids a tiny pivot when a stable alternative fits within the unchanged, stricter pivot tolerance. Original-matrix primal/dual certificates, resource limits and explicit failure behavior remain mandatory. Added primitive and five-period annual/cumulative horizon regressions. No public signature, backend default or R numerical implementation changes.

# pnsbib 0.0.18

Documentation and regression-test release. General finite-data scope, model/method/estimand distinctions, a known generating population and bounds-in-bounds exposition are now explicit. A dependency-free installed tutorial includes executable assertions. Regression tests cover nested information, coupled-input dependence, label invariance and outcome-dependent selection. No numerical algorithms, exported signatures or backend defaults change.

# pnsbib 0.0.17

Explicit native static response-type/base-constraint/event construction. Static exact calls with `lp_backend="zig"` use native construction and optimization. Partial/longitudinal construction and additional band/fractional helpers remain R-based. Default remains reference. Static output cap: 16,777,216 doubles.

# pnsbib 0.0.16

Explicit native LP selection across all public continuous LP routes, inherited from the model and overridable per call. No silent fallback; old saved objects retain reference behavior.

# pnsbib 0.0.15

Native LP pivot policy corrected using the most improving coefficient and lexicographic minimum-ratio selection. Two retained fractional-program pivot-limit failures from 0.0.14 are regression fixtures. Tolerances and 200,000-pivot limit are unchanged; tableau and ordering scratch share the 16,777,216-cell budget.

# pnsbib 0.0.14

Experimental native continuous LP with primal/dual/Farkas/ray diagnostics. Initial cross-route validation found two pivot-limit failures; this version did not establish adoption readiness.

# pnsbib 0.0.13

Full native Shu-Wang-Li comparator, preserved R oracle and labelled repeated-outcome extension. Results remain `valid_closed_form`. Native input table cap: 1,048,576 cells.

# pnsbib 0.0.12

Full printed Li-Pearl Theorems 4-11 recursion in Zig through rzig, with an explicit R oracle. Native memo cap: 1,048,576 states with size `(m+1)*(n+1)*2^k`.

# pnsbib 0.0.11

R graphics and PNG/PDF export with unchanged numerical endpoints and source-table provenance.

# pnsbib 0.0.10

Portable R-header discovery through `R.home("include")`; numerical source unchanged from 0.0.9.

# pnsbib 0.0.9

General Shu-Wang-Li comparator on complete margins, with explicitly labelled repeated factual-outcome extension. Numerical agreement with LP is not a general sharpness proof.

# pnsbib 0.0.8

Li-Pearl Theorems 4/5 labels corrected to `valid_analytical` after complete-information counterexamples. Arithmetic unchanged. Historical saved labels require that qualification.

# pnsbib 0.0.7

Conservative simultaneous sampling-region projection with declared sampling units, fixed weights and structural templates whose numerical margins are ignored.

# pnsbib 0.0.6

Categorical potential paths, ordered absorption, state-specific margins and horizon/query curves. No anticipation compares full states.

# pnsbib 0.0.5

Partial-information model with selected probability equalities/bands, factual-subgroup constraints and jointly weighted disjoint strata.

Earlier releases established the finite analytical/exact interface, longitudinal binary paths and sensitivity calculations.
