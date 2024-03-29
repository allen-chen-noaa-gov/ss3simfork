% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_data.r
\name{clean_data}
\alias{clean_data}
\title{Given sampling arguments remove ("clean") all data in a .dat file that
is not specified}
\usage{
clean_data(
  dat_list,
  index_params = NULL,
  lcomp_params = NULL,
  agecomp_params = NULL,
  calcomp_params = NULL,
  mlacomp_params = NULL,
  verbose = FALSE
)
}
\arguments{
\item{dat_list}{An SS data list object as read in from
\code{\link[r4ss]{SS_readdat}} in the \pkg{r4ss} package. Make sure
you select option \code{section=2}.}

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

\item{verbose}{When \code{TRUE} it will print a message when rows are
deleted.}
}
\value{
An invisible cleaned data list as an object.
}
\description{
This prepares a \code{.dat} file to be used by an EM, whereas before it may
have had leftover data from sampling purposes. See examples in
\code{\link{change_data}}.
}
\note{
This function does not write the result to file.
}
\seealso{
calculate_data_units, change_data

Other sampling functions: 
\code{\link{sample_agecomp}()},
\code{\link{sample_calcomp}()},
\code{\link{sample_index}()},
\code{\link{sample_lcomp}()},
\code{\link{sample_mlacomp}()},
\code{\link{sample_wtatage}()}
}
\author{
Cole Monnahan
}
\concept{sampling functions}
