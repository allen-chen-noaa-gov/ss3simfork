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

moredatadahl <- function(datfile.origsave,dat_list,locnum,obsnum,betavar,
    uparams,abundse=NULL,abundscale=NULL,filename=NULL,...){
dattemp <- sample_index(dat_list        = datfile.origsave,
						outfile         = NULL,
						fleets          = 2,
						years           = list((100-list(...)$obsyears):100),
						sds_obs         = list(0.2),
						write_file      = FALSE)

realcpue <- (datfile.origsave$CPUE[(100-list(...)$obsyears):100,])

if (is.null(abundscale) == TRUE) {
    scaleq <- mean(realcpue$obs)/sum(betavar)
} else {
    scaleq <- abundscale
}

scalecatch <- realcpue$obs/mean(realcpue$obs)

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)

if (is.null(filename) == TRUE) {
    abundfile <- paste0(getwd(), "/", abundtitle, "/cor-abund-", 
        gsub("/", "-", abundtitle), ".csv")
} else {
    abundfile <- paste0(filename,  gsub(".*/","",abundtitle), ".csv")
}

if (file.exists(abundfile) == FALSE) {
    
newabund <- list()
choiceout <- list()
startout <- list()
for (i in 1:length(scalecatch)) { 

set.seed(i)
betavarin <- as.matrix(betavar)*scalecatch[i]
kk <- dim(locnum)[1]

otherdatfin <- spatial_fishery(locnum,obsnum,betavarin,uparams,
    datfile.origsave$catch[(100-list(...)$obsyears):100,],year=i,random=FALSE,
    list(...)$avghauls)

polyn <- 3
polyintnum <- 1
regconstant <- 0
polyconstant <- 1
singlecor <- 0

otherdatfin$polyn <- polyn
otherdatfin$polyintnum <- polyintnum
otherdatfin$regconstant <- regconstant
otherdatfin$polyconstant <- polyconstant
otherdatfin$singlecor <- singlecor
zifin <- do.call(cbind,otherdatfin$intdat)

if (regconstant == 1) {
    initparams <- unname(c(1, betavarin, rep(0, 
        (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
        rep(-1,dim(zifin)[2]), 1))
    #Initial paramters for revenue then cost.
} else {
    initparams <- unname(c(1, betavarin, rep(0, 
        (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
        rep(-1,dim(zifin)[2]), 1))
}

optimOpt <- c(100000,1.00000000000000e-08,1,0) 
#Optimization options for the maximum number of function evaluations, 
#maximum iterations, and the relative tolerance of x. Then, how often to report 
#output, and whether to report output.

methodname = "BFGS"

bw <- -1

otherdatfin$bw <- bw
func <- barebones.FishSET::logit_correction_polyint_estscale

results_savev <- barebones.FishSET::discretefish_subroutine(otherdatfin$catchfin,
    otherdatfin$choicefin,otherdatfin$distance,otherdatfin,initparams,
    optimOpt,func,methodname)
    
initcount <- 0
searchspace <- 1000

if (regconstant == 1) {
changevec <- unname(c(1, rep(0, length(betavarin)),
    rep(1, (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
    rep(1,dim(zifin)[2]), 1)) 
    # Initial paramters for revenue then cost.
} else {
changevec <- unname(c(1, rep(0, length(betavarin)),
    rep(1, (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
    rep(1,dim(zifin)[2]), 1)) 
}

results <- barebones.FishSET::explore_startparams(searchspace, initparams, dev = 2, 
    logit_correction_polyint_estscale, otherdatfin$catchfin, 
    otherdatfin$choicefin, otherdatfin$distance, otherdatfin,
    changevec)

LLmat <- data.frame(cbind(1:searchspace,unlist(results$saveLLstarts)))
LLmatorder <- LLmat[order(LLmat$X2),]

initparamssave <- results$savestarts[LLmatorder$X1[1:100]]

while ((any(is.na(as.numeric(results_savev$OutLogit[,2]))) == TRUE ||
    results_savev$OutLogit[1,1] < 0) & initcount < 20) {

initcount <- initcount + 1
    
initparams <- initparamssave[[initcount]]

results_savev <- barebones.FishSET::discretefish_subroutine(otherdatfin$catchfin,
    otherdatfin$choicefin,otherdatfin$distance,otherdatfin,initparams,
    optimOpt,func,methodname)

}

if (is.numeric(results_savev$OutLogit[2:(kk+1),1]) == FALSE ||
    any(is.na(as.numeric(results_savev$OutLogit[,2]))) == TRUE) {
    newabund[i] <- NA
} else {
    newabund[i] <- sum(results_savev$OutLogit[2:(kk+1),1])*scaleq
}
choiceout[[i]] <- table(otherdatfin$choicefin)
startout[[i]] <- table(otherdatfin$startloc)
}

dattemp$CPUE$obs <- unlist(newabund)
dattemp$CPUE[is.na(dattemp$CPUE$obs)==FALSE,]
dattemp$CPUE$index <- 3

abundout <- merge(realcpue, dattemp$CPUE,  by = c("year", "seas"))
abundout$index.x <- NULL
abundout$se_log.x <- NULL
abundout$index.y <- NULL
names(abundout)[names(abundout) == "obs.x"] <- "TrueCPUE"
names(abundout)[names(abundout) == "obs.y"] <- "EstCPUE"

abundout$diffperc = (abundout$TrueCPUE - abundout$EstCPUE)/abundout$TrueCPUE

abundout <- abundout[order(abundout$year),]

write.table(abundout, file=abundfile, sep=",", row.names=FALSE, quote = FALSE)
write.table(cbind(do.call(rbind, startout), do.call(rbind, choiceout)), 
	file=paste0(getwd(), "/", abundtitle, "/sampout-", 
	gsub("/", "-", abundtitle), ".csv"), 
	sep=",", row.names=FALSE, quote = FALSE)

} else {

abundout <- read.table(abundfile, sep=",", header=TRUE)

abundout <- abundout[order(abundout$year),]

dattemp$CPUE$obs <- abundout$EstCPUE
dattemp$CPUE$index <- 3

}

if (is.null(abundse) == TRUE) {
    dattemp$CPUE$se_log <- sqrt(log(1+((sd(dattemp$CPUE$obs)/
        mean(dattemp$CPUE$obs))^2)))
} else {
    dattemp$CPUE$se_log <- abundse
}

dat_list$CPUE <- rbind(dat_list$CPUE, 
    dattemp$CPUE[is.na(dattemp$CPUE$obs)==FALSE,])

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)

}
