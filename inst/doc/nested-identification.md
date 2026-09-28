# A known population and nested information

This tutorial uses finite binary potential outcomes for two generic policies. It requires no bibliographic records and performs no rank transformation. Run its executable companion with `source(system.file("doc", "nested-identification.R", package="pnsbib"))`.

## Population, observed cells and intervention risks

Let half the population receive each policy. The table gives response-type probabilities conditional on factual assignment, not intervention risks.

| `(Y_reference,Y_active)` | Given reference | Given active |
|---|---:|---:|
| `(0,0)` | .5 | .3 |
| `(0,1)` | .2 | .4 |
| `(1,0)` | .2 | .1 |
| `(1,1)` | .1 | .2 |

The joint response masses are half those entries. PNS is `.5*.2 + .5*.4 = .30`; PN is `.4/(.4+.2) = 2/3`; PS is `.2/(.5+.2) = 2/7`. The factual table is `rbind(c(.35,.15),c(.20,.30))`, and the success risks under reference and active intervention are `.30` and `.45`. The latter use both assignment groups. Factual group success risks are instead `.30` and `.60`.

With the factual table and complete intervention risks, Li-Pearl, Shu-Wang-Li and exact LP all give PNS `[.15,.45]` in this example. The known truth `.30` lies inside. The result illustrates partial identification even with known population inputs. It establishes agreement here, not universal analytical sharpness.

## Bounds within admissible input bounds

For a fixed query `theta(q)`, let `Q(eta,A)` contain the response distributions compatible with input probabilities `eta` and structural assumptions `A`. Inner bounds are `L(eta)=inf theta(q)` and `U(eta)=sup theta(q)` over `Q(eta,A)`. For a joint admissible input region `H`, the outer endpoints are `inf_{eta in H} L(eta)` and `sup_{eta in H} U(eta)`, retaining only compatible pairs `(eta,q)`. The result is an interval hull; a nonconvex region can have gaps that the hull fills.

If `H_small` is contained in `H_large`, their feasible response sets have the same containment. Therefore `L_large <= L_small <= U_small <= U_large` whenever both sets are nonempty and the query is defined. This is the sense of "bounds-in-bounds" used here. It describes existing projection and sensitivity operations. It is not a new causal-identification assumption or a DiD estimator.

The executable example centers bands on known missing-subgroup risks and narrows their widths through `1,.5,.25,.1,0`. Those oracle-centered bands demonstrate geometry, not an empirical way to estimate missing risks. For empirical use, the band center, scale, target and justification must come from the design or an explicitly hypothetical sensitivity scenario.

## Why joint input dependence matters

Let `q00,q01,q10,q11` denote binary response masses and suppose the only intervention restriction is equality of the two success risks. Then `q01=q10`, so `PNS=q01<=.5`. Both zero and `.5` are attained. Each risk separately ranges over `[0,1]`; replacing their equality by two separate intervals discards the dependence and permits PNS equal to one. Searching only the equal-risk endpoints `(0,0)` and `(1,1)` incorrectly gives maximum zero. The maximum occurs at the interior pair `(.5,.5)`.

The regression suite checks this tiny linear program against direct feasible witnesses with the native and reference solvers. The public partial-model constraint API supports its documented probability events and bands, not every arbitrary joint linear input restriction. The low-level native LP interface is not a general scientific region adapter.

## Conditional events and uncertainty labels

PN and PS are ratios for factual conditioning populations. Fixed positive denominators permit scaling optimized numerators. Variable denominators require joint fractional optimization; independently chosen numerator and denominator extrema can describe different populations. Zero denominators stay undefined. Incompatible input regions stay empty and must not disappear from coverage accounting.

An assumption band states what missing risks are allowed. A confidence region requires a sampling argument. A source envelope describes admissible measurements. A posterior region depends on its probability model. Their labels are not interchangeable. Projecting a valid simultaneous confidence region through a correctly specified causal model inherits its coverage; projection alone supplies neither premise.

## Preserve the intended population

Retaining only factual successes in this population gives PNS `4/9` and makes PS undefined because the reference/failure group disappears. Post-outcome filtering changes the scientific target. Preserve zero outcomes when they belong to the target, and keep missing acquisition outcomes distinct from observed zeros. Prespecified baseline eligibility is possible, with all inputs then referring to the same selected population.

Complete-margin static, partial-information static and longitudinal trajectory models all support probability-of-causation questions. Only the supplied information differs. A trajectory model additionally declares regimes, complete potential path support and any justified temporal restrictions. The package treats category labels as names; substantively meaningful thresholds and rankings are constructed outside it.
