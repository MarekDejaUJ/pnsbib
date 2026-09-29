\name{poc_plot}
\alias{poc_plot}
\alias{plot.pnsbib_plot}
\alias{print.pnsbib_plot}
\title{Visualize already-computed probability bounds}
\usage{
poc_plot(x, type = c("interval", "curve"), labels = NULL,
  series = NULL, title = NULL, subtitle = NULL,
  monochrome = FALSE, draw = TRUE)
\method{plot}{pnsbib_plot}(x, ...)
\method{print}{pnsbib_plot}(x, ...)
}
\arguments{
\item{x}{A result data frame or \code{poc_analysis()} object for \code{poc_plot}; a plot object for the methods.}
\item{type}{Horizontal interval comparison or numeric-horizon curves.}
\item{labels}{Optional character vector of interval row labels.}
\item{series}{Optional character vector assigning curve rows to series. Methods, estimands, assumptions and declared scenarios must agree within series. Default grouping uses available semantic fields.}
\item{title, subtitle}{Single character strings; defaults describe probability bounds and supplied input interpretation.}
\item{monochrome}{Use grayscale, retaining redundant marker and line types.}
\item{draw}{Draw immediately; FALSE constructs the plot object only.}
\item{...}{Reserved; additional arguments are rejected. Set options through \code{poc_plot}.}
}
\value{Invisibly, a \code{pnsbib_plot} containing unchanged source rows, display mapping, provenance and presentation settings.}
\description{Displays previously computed bounds on a full probability axis. No solver or margin estimator is called. Undefined results remain labelled and interrupt curves; endpoints are never smoothed or forced to be monotone.}
\details{Curves require an estimand and finite numeric horizon for every row, unique within each series. At most eight series can be drawn together. Input row order is retained; drawing uses increasing horizons. The plot does not establish that different methods received identical models or information: the caller must ensure comparability. Bounds are not automatically sampling confidence intervals. All source metadata survive export even when too detailed for a figure label.}
\seealso{\code{\link{poc_save_plot}}, \code{\link{poc_curve}}, \code{\link{poc_analysis}}}
\examples{
o <- matrix(c(.2, .1, .3, .4), 2, byrow = TRUE,
            dimnames = list(c("off", "on"), c("0", "1")))
a <- matrix(c(.6, .4, .3, .7), 2, byrow = TRUE, dimnames = dimnames(o))
m <- poc_model(o, a)
q <- poc_query(c(off = "0", on = "1"))
p <- poc_plot(poc_analysis(list(poc_bounds(m, q), poc_exact(m, q)),
  provenance = list(inputs = "synthetic example")), draw = FALSE)
plot(p)
}
\references{
See \link{pnsbib-methods} for the method-to-function map, assumptions and
references. Graphics and reporting functions preserve supplied results;
they do not perform causal identification.
}
