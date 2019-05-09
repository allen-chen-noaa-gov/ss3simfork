#put in case-write.R write E case
case_E <- function(natM_type = NULL, natM_n_breakpoints = NULL, natM_lorenzen = NULL, natM_val = NULL,
  par_name = NULL, par_int = NULL, par_phase = NULL, forecast_num = NULL, case, spp) {

  old <- options()$"deparse.cutoff"
  options(deparse.cutoff = 500L)
  on.exit(options(deparse.cutoff = old))

  for (ind in seq_along(spp)){
    filename <- paste0("E", case, "-", spp[ind], ".txt")
    mapply(write, file = filename, MoreArgs = list(ncolumns = 2,
      x = c("natM_type;", natM_type)))
    mapply(write, file = filename, MoreArgs = list(ncolumns = 2, append = TRUE,
      x = c("natM_n_breakpoints;", case_deparse(natM_n_breakpoints))))
    mapply(write, file = filename, MoreArgs = list(ncolumns = 2, append = TRUE,
      x = c("natM_lorenzen;", case_deparse(natM_lorenzen))))
    mapply(write, file = filename, MoreArgs = list(ncolumns = 2, append = TRUE,
      x = c("natM_val;", case_deparse(natM_val))))
	mapply(write, file = filename, MoreArgs = list(ncolumns = 2, append = TRUE,
      x = c("par_name;", par_name)))
	mapply(write, file = filename, MoreArgs = list(ncolumns = 2, append = TRUE,
      x = c("par_int;", case_deparse(par_int))))
	mapply(write, file = filename, MoreArgs = list(ncolumns = 2, append = TRUE,
      x = c("par_phase;", case_deparse(par_phase))))
	mapply(write, file = filename, MoreArgs = list(ncolumns = 2, append = TRUE,
      x = c("forecast_num;", case_deparse(forecast_num))))
  }
}
