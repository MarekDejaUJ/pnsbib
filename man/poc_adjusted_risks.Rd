\name{poc_adjusted_risks}
\alias{poc_adjusted_risks}
\alias{poc_risk_fit}
\title{Native Binary-Risk Estimation and Common-Target Standardization}
\usage{
poc_risk_fit(x, y, target = x, weights = rep(1, nrow(x)),
  lambda = 1, penalty = c(.1, rep(1, ncol(x)-1)),
  max_iterations = 100L, tolerance = 1e-10)
poc_adjusted_risks(x, y, group, target = x,
  weights = rep(1, nrow(x)), target_weights = rep(1, nrow(target)),
  mode = c("annual", "cumulative"), lambda = 1,
  penalty = c(.1, rep(1, ncol(x)-1)))
}
\arguments{
\item{x,target}{Finite numeric training and common-target design matrices. For standardization, the first column must be an intercept of ones. Scale baseline predictors externally using training data only.}
\item{y}{Binary numeric outcomes. Standardization accepts a matrix of annual outcomes, one horizon per column, including for cumulative mode.}
\item{group}{Training labels 0 and 1, with both groups represented.}
\item{weights,target_weights}{Nonnegative finite training frequency weights and target weights. Target weights are normalized in native computation.}
\item{lambda}{Strictly positive ridge coefficient on summed negative log likelihood.}
\item{penalty}{Strictly positive coefficient penalty multipliers.}
\item{max_iterations,tolerance}{Iteration limit and normalized maximum-gradient convergence tolerance.}
\item{mode}{Separate annual risks or cumulative first-success risks estimated from discrete hazards.}
}
\description{An opt-in estimation layer for complete, partial and longitudinal probability-of-causation models. Zig through rzig performs optimization, prediction and weighted common-target standardization. Existing probability-table interfaces remain unchanged. No R numerical fallback is used.}
\details{The criterion is summed weighted binomial negative log likelihood plus one-half lambda times the penalty-weighted squared coefficients. Positive intercept regularization keeps separated fits finite. Cumulative mode fits each group and horizon among researchers with no earlier success and transforms hazards by survival products. Empty training risk sets fail explicitly. Cross-world dependence remains unrestricted. These inputs require consistency, an appropriate intervention, conditional exchangeability, overlap and an interference argument for causal interpretation. No sampling confidence interval is supplied, and compatibility with factual cells must be checked separately.}
\value{For standardization, a list with \code{risks} (two rows labelled 0 and 1, one column per horizon), \code{predictions} (target rows by horizons by groups), nested \code{fits}, \code{mode}, \code{backend}, \code{target_n}, \code{target_weight_sum} and \code{interpretation}. Cumulative predictions are first-success probabilities, not raw hazards. The lower-level fitter returns \code{coefficients}, \code{predictions}, \code{objective}, \code{iterations}, \code{gradient} (maximum absolute gradient divided by the larger of one and total training weight), \code{backend}, \code{lambda} and \code{penalty}.}
\seealso{\code{\link{poc_model}}, \code{\link{poc_partial_model}}, \code{\link{poc_longitudinal_model}}. The installed \code{doc/adjusted-risks.md} and \code{doc/adjusted-risks.R} show all three connections.}
\examples{
x <- cbind(intercept=1, baseline=c(-1,0,1,-1,0,1))
y <- cbind(year1=c(0,0,1,0,1,1), year2=c(0,1,1,1,1,1))
a <- poc_adjusted_risks(x,y,group=c(0,0,0,1,1,1))
a$risks
poc_risk_fit(x, y[, 1])
}

\references{
le Cessie, S. and van Houwelingen, J. C. (1992). Ridge estimators in
logistic regression. \emph{Applied Statistics}, 41(1), 191--201.
\doi{10.2307/2347628}.

Hernan, M. A. and Robins, J. M. (2020). \emph{Causal Inference: What If}.
Chapman & Hall/CRC. \url{https://miguelhernan.org/whatifbook}.

See \link{pnsbib-methods} for the method-to-function map, assumptions and
references. Graphics and reporting functions preserve supplied results;
they do not perform causal identification.
}
