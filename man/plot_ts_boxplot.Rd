% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_ts_boxplot}
\alias{plot_ts_boxplot}
\title{Plot timeseries values as boxplots.}
\usage{
plot_ts_boxplot(
  data,
  y,
  horiz = NULL,
  horiz2 = NULL,
  vert = NULL,
  vert2 = NULL,
  relative.error = FALSE,
  axes.free = TRUE,
  print = TRUE
)
}
\arguments{
\item{data}{A valid data frame containing scalar or timeseries values
from a \pkg{ss3sim} simulation. That data are generated from
\code{\link{get_results_all}}.}

\item{y}{A character string denoting which column to use as the y
variable. Must be a numeric column.}

\item{horiz, horiz2}{A character string denoting which column to use as
the first (\code{horiz}) and second (\code{horiz2}) level of faceting in
the horizontal direction. E.g. "M" or "species". A value of NULL (default)
indicates no faceting.}

\item{vert, vert2}{A character string denoting which column to use as
the first (\code{vert}) and second (\code{vert2}) level of faceting in
the vertical direction. E.g. "M" or "species". A value of NULL (default)
indicates no faceting.}

\item{relative.error}{Boolean for whether the y-axis should be
interpreted as relative error. If \code{TRUE}, \code{ylim} is set to
\code{c(-1,1)}, the y axis label is changed automatically, and a
red line at y=0 is added.}

\item{axes.free}{Boolean for whether the y-axis scales should be free
in \code{facet_grid}.}

\item{print}{A logical for whether the plot is printed or not.}
}
\description{
Plot timeseries values as boxplots.
}
\details{
The \pkg{ss3sim} plotting functions are simply
wrappers for \pkg{ggplot2} code, specific to the output from
\pkg{ss3sim} simulation scalar and timeseries (ts) objects. They are
designed to quickly explore simulation output, rather than
publication-level figures. The functions use the \code{aes_string}
function within \pkg{ggplot2} such that arguments are passed as
characters that refer to columns of \code{data}.

Note that there are some subtle differences between the
functions. Scalar plots require a value for \code{x}, while for ts plots
\code{x} is invalid because it is fixed internally as 'year', since it
makes no sense to use another column. Boxplots cannot have a color
mapped to them like points or lines, and thus \code{color} is not a
valid argument. The ts point and line plots are grouped internally by
'ID', which is a combination of scenario and replicate.
}
\section{Output}{
 These functions print the \code{ggplot} object, but
also return it invisibly for saving or printing again later.
}

\examples{
\dontrun{
ts_dat$SpawnBio <- with(ts_dat, (SpawnBio_om-SpawnBio_em)/SpawnBio_om)
# Merge in max_grad, a performance metric, to use for color
ts_dat <- merge(scalar_dat[, c("ID", "max_grad")], ts_dat)
plot_ts_boxplot(ts_dat, y = "SpawnBio", horiz = "D", vert = "E",
  relative.error = TRUE)
}
}
\author{
Cole Monnahan
}
