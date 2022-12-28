#' readabund
#'
#'
#' @param filename description
#' @param type description
#'
#' @export

readabund <- function(filename, type, istruecpue = FALSE) {

abunddat <- read.csv(filename)

tempabund <- abunddat[order(abunddat$year), ]
tempabund$se_log.y <- NULL
tempabund$diffperc <- NULL

tempabund$Index <- type
if (istruecpue == TRUE) {
names(tempabund)[names(tempabund) == "TrueCPUE"] <- "CPUE"
} else {
names(tempabund)[names(tempabund) == "EstCPUE"] <- "CPUE"
}

tempabund$RelDiff <- tempabund$CPUE/(tempabund$CPUE[1])

tempabund <- tempabund[, c("year", "seas", "Index", "CPUE", "RelDiff")]

return(tempabund)

}
