\name{poc_shu2026}
\alias{poc_shu2026}
\title{Shu-Wang-Li closed-form probability-of-causation comparator}
\usage{poc_shu2026(model, query, use_zig = TRUE)}
\arguments{
\item{model}{A complete-margin model returned by \code{poc_model()}.}
\item{query}{A static query returned by \code{poc_query()}.}
\item{use_zig}{Use the native closed form. FALSE selects the preserved R oracle. The native input table is limited to 1048576 cells.}
}
\value{A one-row \code{pnsbib_result} with the theorem/extension in \code{method}.}
\description{
Evaluates Theorems 2--5 of Shu, Wang and Li (IJCAI 2026) on complete
observational and interventional margins. Distinct-subset upper candidates
are computed by sorted prefix sums, without enumerating response types.
}
\details{
Repeated matches between potential and factual outcome labels use an
explicitly labelled project extension of Theorem 4. Let \eqn{S} be the queried
treatments whose requested outcome equals the factual outcome, \eqn{C} the
observed mass of compatible factual cells, and for each queried treatment
let \eqn{a_j=P(Y_{x_j}=y_j)}, \eqn{o_j=P(X=x_j,Y=y_j)},
\eqn{u_j=P(X=x_j)}, \eqn{d_j=a_j-o_j}, \eqn{b_j=d_j+u_j}.
For \eqn{k} potential requirements, lower candidates are zero,
\eqn{C+\sum_j b_j-k}, and
\eqn{o_i+\sum_{j\ne i}b_j-k+1} for \eqn{i\in S}.
Upper candidates are \eqn{C}, \eqn{a_i} for \eqn{i\in S},
\eqn{d_j} for \eqn{j\notin S}, and
\eqn{\sum_{i\in T}d_i/(|T|-1)} for distinct subsets \eqn{T\subseteq S}
of size at least two. The last family follows by summing
\eqn{P(E)-P(E,X=x_i)\le d_i} over disjoint factual treatments.
Union bounds over off-own-treatment failures establish the lower candidates.

All intervals are labelled \code{valid_closed_form}; general sharpness is not
asserted. Numerical agreement with an exact LP is a case-specific check.
Consistency reductions retain the original query's conditioning denominator.
Rectangular state spaces and repeated potential labels are supported.
No causal identification or sampling-uncertainty claim follows from the inputs.
}
\references{
Shu, X., Wang, S., and Li, A. (2026). Identification of Probabilities of
Causation: From Recursive to Closed-Form Bounds. IJCAI, 4063--4070.
\doi{10.24963/ijcai.2026/452}.
}
\seealso{\code{\link{poc_bounds}}, \code{\link{poc_exact}}}
