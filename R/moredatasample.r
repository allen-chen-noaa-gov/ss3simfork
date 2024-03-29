#' moredatasample
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatasample <- function(datfile.origsave, dat_list, locnum, obsnum, betavar,
  uparams, abundse = NULL, catchscale = NULL, filename = NULL, ...) {
dattemp <- sample_index(dat_list        = datfile.origsave,
  outfile         = NULL,
  fleets          = 2,
  years           = list((100 - list(...)$obsyears):100),
  sds_obs         = list(0.2),
  write_file      = FALSE)

realcpue <- (datfile.origsave$CPUE[(100 - list(...)$obsyears):100, ])

if (is.null(catchscale) == TRUE) {
  catchscale <- min(
    datfile.origsave$catch[(100 - list(...)$obsyears):100, ])/1000
}

scaleabund <- realcpue$obs/mean(realcpue$obs)

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)

if (is.null(filename) == TRUE) {
  abundfile <- paste0(getwd(), "/", abundtitle, "/bias-abund-",
      gsub("/", "-", abundtitle), ".csv")
} else {
  abundfile <- paste0(filename,  gsub(".*/", "", abundtitle), ".csv")
}

if (file.exists(abundfile) == FALSE) {

newabund <- list()
newse <- list()
for (i in seq_along(scaleabund)) {

# normalize abundance over locations so only the OM trend affects relative
# abundance, not spatial size of the fishery, to compare scenarios. Could
# generalize to let spatial size matter in the future.
if (length(betavar) == 1) {
tempb <- runif(betavar, 0.75, 1.50)
betavarscaled <- (tempb/sum(tempb)) * 10.125
} else {
betavarscaled <- (betavar/sum(betavar)) * 10.125
}

if (list(...)$trend == TRUE) {
slope <- (max(betavar) - min(betavar))/list(...)$obsyears
slopesign <- sign(betavar[1] - betavar[length(betavar)])
betaslope <- seq((-slopesign) * slope, slopesign * slope, by = (
  slopesign * (2 * slope)/(length(betavar) - 1)))
betavarscaled <- ((betavar/sum(betavar)) * 10.125) + (betaslope * i)
}

betavarin <- as.matrix(betavarscaled) * scaleabund[i]
kk <- dim(locnum)[1]

otherdatfin <- spatial_fishery(locnum, obsnum, betavarin, uparams,
  datfile.origsave$catch[(100 - list(...)$obsyears):100, ], year = i,
  random = TRUE, list(...)$avghauls, catchscale, list(...)$catchvarV,
  list(...)$catchvarN)

choicefin <- otherdatfin$choicefin
sifin <- do.call(cbind, otherdatfin$griddat)
catchfin <- otherdatfin$catchfin

XX <- model.matrix(~as.factor(V1) - 1, choicefin) * sifin
YY <- catchfin$V1

results_savev <- lm(YY ~ XX - 1)

fin <- c(catchscale)
bstring <- paste0("~(", paste(paste("x", 1:(kk), sep = ""), collapse = "+"),
  ")*(%f)")
form <- do.call(sprintf, c(fmt = bstring, as.list(fin)))
seout <- msm::deltamethod(as.formula(form), coef(results_savev),
    vcov(results_savev))

newabund[i] <- mean(results_savev$coef) * catchscale
newse[i] <- sqrt(log(1 + (((seout)/(newabund[[i]]))^2)))

}

dattemp$CPUE$obs <- unlist(newabund)
dattemp$CPUE$se_log <- unlist(newse)
dattemp$CPUE <- dattemp$CPUE[is.na(dattemp$CPUE$obs) == FALSE, ]
dattemp$CPUE$index <- 3

abundout <- merge(realcpue, dattemp$CPUE,  by = c("year", "seas"))
abundout$index.x <- NULL
abundout$se_log.x <- NULL
abundout$index.y <- NULL
names(abundout)[names(abundout) == "obs.x"] <- "TrueCPUE"
names(abundout)[names(abundout) == "obs.y"] <- "EstCPUE"

abundout$diffperc <- (abundout$TrueCPUE - abundout$EstCPUE)/abundout$TrueCPUE

abundout <- abundout[order(abundout$year), ]

write.table(abundout, file = abundfile, sep = ",", row.names = FALSE,
  quote = FALSE)

} else {

abundout <- read.table(abundfile, sep = ",", header = TRUE)

abundout <- abundout[order(abundout$year), ]

}

dat_list$CPUE <- rbind(dat_list$CPUE,
    dattemp$CPUE[is.na(dattemp$CPUE$obs) == FALSE, ])

rownames(dat_list$CPUE) <- seq_len(nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index == 3, ])[1]

return(dat_list)

}
