#' moredatabias
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatabias <- function(datfile.origsave, dat_list, locnum, obsnum, betavar,
  uparams, abundse = NULL, catchscale = NULL, filename = NULL, ...) {

sizetrue <- 1

dattemp <- sample_index(dat_list        = datfile.origsave,
  outfile         = NULL,
  fleets          = 2,
  years           = list((100 - list(...)$obsyears):100),
  sds_obs         = list(0.2),
  write_file      = FALSE)

realcpue <- (datfile.origsave$CPUE[(100 - list(...)$obsyears):100, ])

matching_years <- realcpue$year
reallength <- datfile.origsave$lencomp[datfile.origsave$lencomp$Yr %in%
  matching_years & datfile.origsave$lencomp$FltSvy == 1, ]
num_colnames <- grep("^l\\d+$", colnames(reallength), value = TRUE)
reallength <- reallength[, num_colnames, drop = FALSE]
#num_colnames_numeric <- as.numeric(sub("^l", "", num_colnames))

# start with static fish site preference
# start with static abundance scenario one
# size preference in fisher utility

# do not collect length data
# use case length data 
# obviously don't have feedback bc om is fixed

conc <- 5
# four areas, each element (column) corresponds to an area
# each prop corresponds to one length bin
# prop1 <- c(5,1,1,1)

proplist <- list()
dirchlist <- list()
for (i in 1:dim(reallength)[2]) {
  v <- rep(1, dim(locnum)[1])
  pos <- ceiling(i * dim(locnum)[1] / dim(reallength)[2])
  v[pos] <- 5
  proplist[[i]] <- v
  # each x corresponds to the proportion of a length in an area
  # each prop corresponds to one length bin
  dirchlist[[i]] <- MCMCpack::rdirichlet(1, alpha = v * conc)
}

# columns sum to 1 so no biomass lost
mat <- do.call(rbind, dirchlist) |> t()
# after transpose each row is an area
# and each column is a length bin
# colSums(mat)

if (is.null(catchscale) == TRUE) {
  catchscale <- min(datfile.origsave$catch[
    (100 - list(...)$obsyears):100, ])/1000
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

nal <- as.numeric(reallength[i, ])
ones <- rep(1, nrow(mat))
# distribution by length and area
# X is then the numbers for a length class (column), where each row is an area
Xlengths <- mat * ones %*% t(nal)
# no fish lost
# colSums(X)-nal
# X[,2] - mat[,2]*nal[2]

Xlengths_norm <- sweep(Xlengths, 1, rowSums(Xlengths), `/`)

linear_vec <- seq(0, 2, length.out = 45)
linear_vec <- linear_vec / linear_vec[23]

weighted_Xlengths <- sweep(Xlengths_norm, 2, linear_vec, `*`)

# For each row of Xlengths_norm, calculate the expected (average) column number
# this is the average size
avg_col_num <- apply(Xlengths_norm, 1, function(prob_row) {
  sum(prob_row * seq_len(ncol(Xlengths_norm)))
})

avg_price <- rowSums(weighted_Xlengths)

if (sizetrue == 1) {
  betavar <- rowSums(Xlengths)
}

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
    random = FALSE,
  list(...)$avghauls, catchscale, list(...)$catchvarV, list(...)$catchvarN,
  avg_price)

profitrowset <- c(otherdatfin$profitfin > 10)
otherdatfin$startloc <- as.matrix(otherdatfin$startloc[profitrowset, ])
otherdatfin$choicefin <- data.frame(V1 = otherdatfin$choicefin[profitrowset, ])
otherdatfin$catchfin <- data.frame(V1 = otherdatfin$catchfin[profitrowset, ])
otherdatfin$distance <- otherdatfin$distance[profitrowset, ]
otherdatfin$intdat[[1]] <- as.matrix(otherdatfin$intdat[[1]][profitrowset, ])
otherdatfin$griddat[[1]] <- otherdatfin$griddat[[1]][profitrowset, ]

# if there are too few observations we cannot estimate locations
set <- as.numeric(unlist(dimnames(table(otherdatfin$choicefin)[
  table(otherdatfin$choicefin) >= list(...)$minobs])))
while(all(table(otherdatfin$choicefin) >= list(...)$minobs) == FALSE ||
  all(table(otherdatfin$startloc) >= list(...)$minobs) == FALSE ||
  length(set) < dim(otherdatfin$distance)[2]) {
choiceset <- as.numeric(unlist(dimnames(table(otherdatfin$choicefin)[
  table(otherdatfin$choicefin) >= list(...)$minobs])))
startset <- as.numeric(unlist(dimnames(table(otherdatfin$startloc)[
  table(otherdatfin$startloc) >= list(...)$minobs])))
set <- intersect(startset, choiceset)
rowset <- !((!(otherdatfin$choicefin$V1 %in% set)) |
  (!(otherdatfin$startloc %in% set)))

otherdatfin$startloc <- as.matrix(otherdatfin$startloc[rowset, ])
otherdatfin$choicefin <- data.frame(V1 = otherdatfin$choicefin[rowset, ])
otherdatfin$catchfin <- data.frame(V1 = otherdatfin$catchfin[rowset, ])
otherdatfin$distance <- otherdatfin$distance[rowset, ]
otherdatfin$distance <- otherdatfin$distance[, set]
otherdatfin$intdat[[1]] <- as.matrix(otherdatfin$intdat[[1]][rowset, ])
otherdatfin$griddat[[1]] <- otherdatfin$griddat[[1]][rowset, ]
otherdatfin$griddat[[1]] <- otherdatfin$griddat[[1]][, set]
}
# end of areas far away

kk <- length(betavarin[set])
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
