#' moredatadahl
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @importFrom barebones.FishSET logit_correction_polyint_estscale
#' @importFrom barebones.FishSET discretefish_subroutine
#' @importFrom barebones.FishSET explore_startparams
#'
#' @export

moredatadahl <- function(datfile.origsave, dat_list, locnum, obsnum, betavar,
  uparams, abundse = NULL, catchscale = NULL, filename = NULL, ...) {

sizetrue <- 0

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

# conc <- 5
# # four areas, each element (column) corresponds to an area
# # each prop corresponds to one length bin
# # prop1 <- c(5,1,1,1)

# proplist <- list()
# dirchlist <- list()
# for (i in 1:dim(reallength)[2]) {
#   v <- rep(1, dim(locnum)[1])
#   pos <- ceiling(i * dim(locnum)[1] / dim(reallength)[2])
#   v[pos] <- 5
#   proplist[[i]] <- v
#   # each x corresponds to the proportion of a length in an area
#   # each prop corresponds to one length bin
#   dirchlist[[i]] <- MCMCpack::rdirichlet(1, alpha = v * conc)
# }

# # columns sum to 1 so no biomass lost
# mat <- do.call(rbind, dirchlist) |> t()
# # after transpose each row is an area
# # and each column is a length bin
# # colSums(mat)

if (is.null(catchscale) == TRUE) {
  catchscale <- min(datfile.origsave$catch[
    (100 - list(...)$obsyears):100, ])/1000
}

scaleabund <- realcpue$obs/mean(realcpue$obs)

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)

if (is.null(filename) == TRUE) {
  abundfile <- paste0(getwd(), "/", abundtitle, "/cor-abund-",
      gsub("/", "-", abundtitle), ".csv")
} else {
  abundfile <- paste0(filename,  gsub(".*/", "", abundtitle), ".csv")
}

if (file.exists(abundfile) == FALSE) {

newabund <- list()
newse <- list()
choiceout <- list()
startout <- list()
paramsout <- list()
trueout <- list()
for (i in seq_along(scaleabund)) {

# nal <- as.numeric(reallength[i, ])
# ones <- rep(1, nrow(mat))
# # distribution by length and area
# # X is then the numbers for a length class (column), where each row is an area
# Xlengths <- mat * ones %*% t(nal)
# # no fish lost
# # colSums(X)-nal
# # X[,2] - mat[,2]*nal[2]

# Xlengths_norm <- sweep(Xlengths, 1, rowSums(Xlengths), `/`)

# # Create a linear vector of length 45, first element 0, median element 1
# linear_vec <- seq(0, 2, length.out = 45)
# linear_vec <- linear_vec / linear_vec[23]

# weighted_Xlengths <- sweep(Xlengths_norm, 2, linear_vec, `*`)

# # For each row of Xlengths_norm, calculate the expected (average) column number
# # this is the average size
# avg_col_num <- apply(Xlengths_norm, 1, function(prob_row) {
#   sum(prob_row * seq_len(ncol(Xlengths_norm)))
# })

# avg_price <- rowSums(weighted_Xlengths)

# if (sizetrue == 1) {
#   betavar <- rowSums(Xlengths)
# }

# normalize abundance over locations so only the OM trend affects relative
# abundance, not spatial size of the fishery, to compare scenarios. Could
# generalize to let spatial size matter in the future.
if (length(betavar) == 1) {
tempb <- runif(betavar, 0.75, 1.50)
betavarscaled <- (tempb/sum(tempb)) * 10.125
tempb <- runif(betavar, 0.75, 1.50)
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

otherdatfin <- spatial_fishery(locnum, obsnum, betavarin, uparams,
  datfile.origsave$catch[(100 - list(...)$obsyears):100, ], year = i,
  random = FALSE,
  list(...)$avghauls, catchscale, list(...)$catchvarV, list(...)$catchvarN)
polyn <- list(...)$polyn
polyintnum <- list(...)$polyintnum
regconstant <- 0
polyconstant <- 1
singlecor <- list(...)$singlecor

otherdatfin$polyn <- polyn
otherdatfin$polyintnum <- polyintnum
otherdatfin$regconstant <- regconstant
otherdatfin$polyconstant <- polyconstant
otherdatfin$singlecor <- singlecor
zifin <- do.call(cbind, otherdatfin$intdat)

#6 is about 5 percent
profitrowset <- c(otherdatfin$profitfin > 6)
otherdatfin$startloc <- as.matrix(otherdatfin$startloc[profitrowset, ])
otherdatfin$choicefin <- data.frame(V1 = otherdatfin$choicefin[profitrowset, ])
otherdatfin$catchfin <- data.frame(V1 = otherdatfin$catchfin[profitrowset, ])
otherdatfin$distance <- otherdatfin$distance[profitrowset, ]
otherdatfin$intdat[[1]] <- as.matrix(otherdatfin$intdat[[1]][profitrowset, ])
otherdatfin$griddat[[1]] <- otherdatfin$griddat[[1]][profitrowset, ]

# for not estimating areas far away with few observations
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
otherdatfin$pricefin <- data.frame(V1 = otherdatfin$pricefin[rowset, ])
otherdatfin$distance <- otherdatfin$distance[rowset, ]
otherdatfin$distance <- otherdatfin$distance[, set]
otherdatfin$intdat[[1]] <- as.matrix(otherdatfin$intdat[[1]][rowset, ])
otherdatfin$griddat[[1]] <- otherdatfin$griddat[[1]][rowset, ]
otherdatfin$griddat[[1]] <- otherdatfin$griddat[[1]][, set]
}

otherdatfin$choicefin$V1 <- as.factor(otherdatfin$choicefin$V1)
levels(otherdatfin$choicefin$V1) <-
  seq_along(unique(levels(as.factor(otherdatfin$choicefin$V1))))
otherdatfin$choicefin$V1 <- as.integer(otherdatfin$choicefin$V1)
# end of areas far away

kk <- length(betavarin[set])
if (regconstant == 1) {
  initparams <- unname(c(1, betavarin[set], rep(0,
    (((polyn + polyconstant) * (1 + (1 - singlecor))) + polyintnum) * kk),
    rep(-1, dim(zifin)[2]), 1))
  #Initial paramters for revenue then cost.
} else {
  initparams <- unname(c(1, betavarin[set], rep(0,
    (((polyn + polyconstant) * (1 + (1 - singlecor))) + polyintnum) * kk),
    rep(-1, dim(zifin)[2]), 1))
}

optimOpt <- c(100000, 1.00000000000000e-08, 1, 0)
#Optimization options for the maximum number of function evaluations,
#maximum iterations, and the relative tolerance of x. Then, how often to report
#output, and whether to report output.

methodname <- "BFGS"

bw <- -1

otherdatfin$bw <- bw
func <- barebones.FishSET::logit_correction_polyint_estscale

results_savev <- barebones.FishSET::discretefish_subroutine(
otherdatfin$catchfin,
otherdatfin$choicefin, otherdatfin$distance, otherdatfin, initparams,
optimOpt, func, methodname)

initcount <- 0
searchspace <- 1000

if (regconstant == 1) {
changevec <- unname(c(1, rep(0, length(betavarin[set])),
  rep(1, (((polyn + polyconstant) * (1 + (1 - singlecor))) + polyintnum) * kk),
  rep(1, dim(zifin)[2]), 1))
    # Initial paramters for revenue then cost.
} else {
changevec <- unname(c(1, rep(0, length(betavarin[set])),
  rep(1, 
    (((polyn + polyconstant) * (1 + (1 - singlecor))) + polyintnum) * kk),
  rep(1, dim(zifin)[2]), 1))
}

results <- barebones.FishSET::explore_startparams(searchspace, initparams,
  dev = 2,
  logit_correction_polyint_estscale, otherdatfin$catchfin,
  otherdatfin$choicefin, otherdatfin$distance, otherdatfin,
  changevec)

LLmat <- data.frame(cbind(1:searchspace, unlist(results$saveLLstarts)))
LLmatorder <- LLmat[order(LLmat$X2), ]

initparamssave <- results$savestarts[LLmatorder$X1[1:100]]

while (any(is.na(as.numeric(results_savev$OutLogit[,2]))) == TRUE &&
  initcount < 10) {

initcount <- initcount + 1

initparams <- initparamssave[[initcount]]

results_savev <- barebones.FishSET::discretefish_subroutine(
  otherdatfin$catchfin,
  otherdatfin$choicefin, otherdatfin$distance, otherdatfin, initparams,
  optimOpt, func, methodname)

}

fin <- c(catchscale)
bstring <- paste0("~(", paste(paste("x", 2:(kk + 1), sep = ""), collapse = "+"),
  ")*(%f)")
form <- do.call(sprintf, c(fmt = bstring, as.list(fin)))
seout <- msm::deltamethod(as.formula(form), results_savev$OutLogit[, 1],
  results_savev$H1)

if (is.numeric(results_savev$OutLogit[2:(kk + 1), 1]) == FALSE ||
  any(is.na(as.numeric(results_savev$OutLogit[, 2]))) == TRUE) {
  newabund[i] <- NA
  newse[i] <- NA
} else {
  newabund[i] <- mean(results_savev$OutLogit[2:(kk + 1), 1]) * catchscale
  newse[i] <- sqrt(log(1 + (((seout)/(newabund[[i]]))^2)))
}
paramsout[[i]] <-  results_savev$OutLogit[2:(kk + 1), 1]
trueout[[i]] <- t(betavarin[set])
choiceout[[i]] <- table(otherdatfin$choicefin)
startout[[i]] <- table(otherdatfin$startloc)
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
write.table(cbind(do.call(rbind, lapply(startout, `[`, 1:length(betavar))), 
  do.call(rbind, lapply(choiceout, `[`, 1:length(betavar)))),
  file = paste0(getwd(), "/", abundtitle, "/sampout-",
  gsub("/", "-", abundtitle), ".csv"),
  sep = ",", row.names = FALSE, quote = FALSE)
write.table(cbind(do.call(rbind, lapply(trueout, `[`, 1:length(betavar))), 
  do.call(rbind, lapply(paramsout, `[`, 1:length(betavar)))),
  file = paste0(getwd(), "/", abundtitle, "/paramsout-",
  gsub("/", "-", abundtitle), ".csv"),
  sep = ",", row.names = FALSE, quote = FALSE)

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
