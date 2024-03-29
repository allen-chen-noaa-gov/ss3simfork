% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_data.r
\name{calculate_data_units}
\alias{calculate_data_units}
\title{Given sampling arguments, calculate super set of fleets, years, and data
types.}
\usage{
calculate_data_units(
  index_params = NULL,
  lcomp_params = NULL,
  agecomp_params = NULL,
  calcomp_params = NULL,
  mlacomp_params = NULL,
  wtatage_params = NULL
)
}
\arguments{
\item{index_params}{Named lists containing the arguments for
\code{sample_index}.}

\item{lcomp_params}{Named lists containing the arguments for
\code{\link{sample_lcomp}}.}

\item{agecomp_params}{Named lists containing the arguments for
\code{\link{sample_agecomp}}.}

\item{calcomp_params}{Named lists containing the arguments for
\code{\link{sample_calcomp}}.}

\item{mlacomp_params}{Named lists containing the arguments for
\code{\link{sample_mlacomp}}.}

\item{wtatage_params}{Named lists containing the arguments for
\code{\link{sample_wtatage}}.}
}
\value{
An invisible list of fleets, years, and types.
}
\description{
Given sampling arguments, calculate super set of fleets, years, and data
types.
}
\note{
A superset by nature is larger than the individual sets used to
create it (unless all sampling arguments are identical), so that the
returned list will created some unnecessary combinations. This was done
intentionally for simplicity but may be changed later. See the vignette
for further information.

See further examples in \code{\link{change_data}}.
}
\examples{
## Should throw error since nothing passed
# calculate_data_units()
## Only one fleet
calculate_data_units(lcomp_params=list(fleets=1, years=c(3,4,6)))
## Add new fleet
calculate_data_units(lcomp_params=list(fleets=1, years=c(3,4,6)),
                     agecomp_params=list(fleets=2, years=5))
## If CAL data called, need other types even if not specified
calculate_data_units(calcomp_params=list(fleets=1, years=c(3,4,6)))
}
\seealso{
clean_data, change_data
}
\author{
Cole Monnahan
}
