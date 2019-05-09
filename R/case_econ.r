#put in case-write.R write econ case
case_econ <- function(functionname = NULL, case, spp) {

  old <- options()$"deparse.cutoff"
  options(deparse.cutoff = 500L)
  on.exit(options(deparse.cutoff = old))

  for (ind in seq_along(spp)){
    filename <- paste0("econ", case, "-", spp[ind], ".txt")
    mapply(write, file = filename, MoreArgs = list(ncolumns = 2,
      x = c("functionname;", functionname)))
  }
}