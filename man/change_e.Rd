% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_e.r
\name{change_e}
\alias{change_e}
\title{Methods to alter which parameters are estimated in a SS3 \code{.ctl} file.}
\usage{
change_e(
  ctl_file_in = "em.ctl",
  ctl_file_out = "em.ctl",
  dat_list = NULL,
  for_file_in = "forecasts.ss",
  natM_type = "1Parm",
  natM_n_breakpoints = NULL,
  natM_lorenzen = NULL,
  natM_val = c(NA, NA),
  par_name = NULL,
  par_int = "NA",
  par_phase = "NA",
  forecast_num = 0,
  run_change_e_full = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{ctl_file_in}{Input SS3 control file}

\item{ctl_file_out}{Output SS3 control file}

\item{dat_list}{An SS data list object as read in from
\code{\link[r4ss]{SS_readdat}} in the \pkg{r4ss} package. Make sure
you select option \code{section=2}.}

\item{for_file_in}{Input SS3 forecast file}

\item{natM_type}{*A character string corresponding to option 0:4 in SS3 (i.e.
"1Parm", "n_breakpoints", "Lorenzen", "agespecific",
"agespec_withseasinterpolate"). A value of \code{NA} will leave the
configuration of natural mortality as specified in \code{ctl_file_in}.}

\item{natM_n_breakpoints}{*A vector of ages at which you want breakpoints.
Only used if you specify \code{natM_type = "n_breakpoints"}.}

\item{natM_lorenzen}{*The reference age for the Lorenzen function.  Only used
if you specify \code{natM_type = "Lorenzen"}. Length should be one even
if the \code{.ctl} has two genders.}

\item{natM_val}{*A vector of numeric values. Interpretation of the values are
dependent upon \code{natM_type}. If \code{natM_type = "agespecific"} or
\code{natM_type = "agespec_withseasinterpolate"} the vector specifies the
fixed natural mortality parameters for each integer age.  Specify values
in the following order: female area 1, female area 2, male area 1, male
area, etc. Ensure that there is one entry per integer age x area x gender.
If \code{natM_type = "1Param"}, \code{natM_type = "n_breakpoints"}, or
\code{natM_type = "Lorenzen"} the vector specifies the initial and phase
values for each natM parameter (i.e. \code{c(int, phase, int, phase,
etc.)}, where the first two values could correspond to ages 0-2 natural
mortality and the third and fourth value could correspond to ages 3-8
natural mortality).  For any specified initial value, the parameter bounds
will be altered to 50 percent above and below the specified initial value,
if the initial value lies above or below the original bounds.}

\item{par_name}{*A vector of values, separated by commas.  Each value
corresponds to a parameter that you wish to turn on or off in the
\code{ctl_file_in}. The values will later be turned into character values
and used to search for specific lines for each parameter in the
\code{ctl_file_in}, therefore it is best to use full parameter names as
they are specified in \code{ctl_file_in}.}

\item{par_int}{*A vector of initial values, one for each parameter in
\code{par_name}.  Values can be \code{NA} if you do not wish to change the
initial value for a given parameter.}

\item{par_phase}{*A vector of phase values, one for each parameter in
\code{par_name}.  Values can be \code{NA} if you do not wish to change
the phase for a given parameter.}

\item{forecast_num}{*Number of years to perform forecasts. For those years,
the data will be removed from the \code{dat_list}, enabling SS3 to
generate forecasts rather than use the data to fit the model.}

\item{run_change_e_full}{*If \code{FALSE} \code{change_e} will only
manipulate for forecasting, if \code{TRUE} (default) the full function
capability will be ran.}

\item{verbose}{When \code{TRUE} messages will be returned
from the function. Often useful for debugging. The default
is \code{FALSE}.}
}
\value{
Altered versions of SS3 \code{.ctl} and \code{forecast.ss} files are written
to the disk and the altered \code{dat_list} is returned invisibly.
}
\description{
Takes SS3 \code{.ctl} and \code{forecast.ss} files, along with
  a list structure which houses the data file as read in by
  \code{\link[r4ss]{SS_readdat}}
  and changes which parameters are estimated, how natural mortality is
  estimated, and if forecasts are performed. The function can be called by
  itself or within \code{\link{run_ss3sim}} to alter an estimation model
  \code{.ctl} file.
  If used with \code{\link{run_ss3sim}} the case file should be named
  \code{E}. A suggested (default) case letter is \code{E} for estimation.
}
\details{
Turning parameters on and off is the main function of
  \code{change_e}.  \code{change_e} was not created with the capability of
  adding parameters to a \code{.ctl} file.  The function can only add
  parameters for age specific natural mortality, and only for models with
  one growth morph.  Furthermore, the function is designed to add complexity
  to the natural mortality type and not remove complexity.  Therefore, the
  function will fail if natural mortality in the \code{ctl_file_in} is not
  specified as \code{"1Param"} and \code{natM_type} is anything other than
  \code{NULL} or \code{"1Param"}.
}
\section{Which arguments to specifiy in case files}{

All function argument descriptions that start with an asterisk (*) will be passed
through the case files to \code{\link{run_ss3sim}}. If one of these arguments
is not specified in a case file, then a value of \code{NULL} will be passed,
which may or may not be an appropriate value. Other arguments will be ignored
if specified.
}

\examples{
\dontrun{
library(r4ss)
# Create a temporary folder for the output and set the working directory:
temp_path <- file.path(tempdir(), "ss3sim-tv-example")
dir.create(temp_path, showWarnings = FALSE)
wd <- getwd()
setwd(temp_path)

d <- system.file("extdata", package = "ss3sim")
ctl_file <- paste0(d, "/models/cod-om/codOM.ctl")
data.old <- r4ss::SS_readdat(file.path(d, "models", "cod-om", "codOM.dat"))
change_e(ctl_file_in = ctl_file, ctl_file_out = "change_e.ctl",
         dat_list = data.old, for_file_in = "forecast.ss",
         natM_type = "n_breakpoints", natM_n_breakpoints = c(1, 4),
         natM_lorenzen = NULL, natM_val = c(.2, 3, 0.4, 5),
         par_name = c("_steep", "SizeSel_1P_1_Fishery"),
         par_int = c(0.3, 40), par_phase = c(3, 2),
         forecast_num = 0, run_change_e_full = TRUE )
# clean up
file.remove("change_e.ctl")
setwd(wd)
}
}
\seealso{
Other change functions: 
\code{\link{change_data}()},
\code{\link{change_em_binning}()},
\code{\link{change_f}()},
\code{\link{change_maturity}()},
\code{\link{change_retro}()},
\code{\link{change_tv}()}
}
\author{
Kelli Johnson
}
\concept{change functions}
