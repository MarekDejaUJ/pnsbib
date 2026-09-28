\name{poc_save_plot}
\alias{poc_save_plot}
\title{Save probability-bound graphics and source tables}
\usage{
poc_save_plot(x, filename, formats = c("png", "pdf"),
  width = 8, height = NULL, dpi = 600)
}
\arguments{
\item{x}{A \code{pnsbib_plot} from \code{poc_plot()}.}
\item{filename}{Filename prefix without a PNG/PDF extension in an existing directory. Percent format sequences are not accepted.}
\item{formats}{One or both of \code{png} and \code{pdf}.}
\item{width, height}{Dimensions in inches. Height defaults to 5 for curves or adapts to interval rows. Width must be 5--30; height 3--60.}
\item{dpi}{PNG resolution, 72--1200, default 600. PNG exports are capped at 100 million pixels. PDF graphics remain vector.}
}
\value{Invisibly, named paths to images, source-data/layout/provenance TSVs and the serialized plot RDS.}
\description{Exports figures with their complete source tables and presentation object. All destination collisions are checked before writing and existing files are never overwritten. Each opened graphics device is closed, including on drawing errors; the caller's device is preserved.}
\details{No numerical engine is invoked. Partial new files may remain after an I/O or graphics error; use a fresh prefix after diagnosing the error. Raster and PDF formats can differ in font rendering, but consume the identical source rows.}
\seealso{\code{\link{poc_plot}}}
