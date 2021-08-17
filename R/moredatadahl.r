#' moredatadahl
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatadahl <- function(datfile.origsave,dat_list,locnum,obsnum,betavar,
    uparams,...){
dattemp <- sample_index(dat_list        = datfile.origsave,
						outfile         = NULL,
						fleets          = 2,
						years           = list(74:100),
						sds_obs         = list(0.2),
						write_file      = FALSE)
						
realcpue <- (datfile.origsave$CPUE[74:100,])

scaleq <- mean(realcpue$obs)/sum(betavar)
		
scalecatch <- realcpue$obs/mean(realcpue$obs)

################################################################################
################################################################################
################################################################################

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)

abundfile <- paste0(getwd(), "/", abundtitle, "/rel-se-abund-", 
    gsub("/", "-", abundtitle), ".csv")

if (file.exists(abundfile) == FALSE) {
    
newabund <- list()
for (i in 1:length(scalecatch)) { 

library(barebones.FishSET)

betavar <- as.matrix(betavar)*scalecatch[i]
kk <- dim(locnum)[1]

otherdatfin <- spatial_fishery(locnum,obsnum,betavar,uparams)

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
    initparams <- unname(c(1, betavar, rep(0, 
        (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
        rep(-1,dim(zifin)[2]), 1))
    #Initial paramters for revenue then cost.
} else {
    initparams <- unname(c(1, betavar, rep(0, 
        (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
        rep(-1,dim(zifin)[2]), 1))
}

optimOpt <- c(100000,1.00000000000000e-08,1,0) 
#Optimization options for the maximum number of
#function evaluations, maximum iterations, and the relative tolerance of x.
#Then, how often to report output, and whether to report output.

methodname = "BFGS"

bw <- -1

otherdatfin$bw <- bw

func <- logit_correction_polyint_estscale

results_savev <- discretefish_subroutine(otherdatfin$catchfin,
    otherdatfin$choicefin,otherdatfin$distance,otherdatfin,initparams,
    optimOpt,func,methodname)
    
initcount <- 0
searchspace <- 1000

if (regconstant == 1) {
changevec <- unname(c(1, rep(0, length(betavar)),
    rep(1, (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
    rep(1,dim(zifin)[2]), 1)) 
    # Initial paramters for revenue then cost.
} else {
changevec <- unname(c(1, rep(0, length(betavar)),
    rep(1, (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
    rep(1,dim(zifin)[2]), 1)) 
}

results <- explore_startparams(searchspace, initparams, dev = 2, 
    logit_correction_polyint, otherdatfin$catchfin, otherdatfin$choicefin, 
    otherdatfin$distance, otherdatfin,
    changevec)

LLmat <- data.frame(cbind(1:searchspace,unlist(results$saveLLstarts)))
LLmatorder <- LLmat[order(LLmat$X2),]

initparamssave <- results$savestarts[LLmatorder$X1[1:100]]
browser()

while ((any(is.na(as.numeric(results_savev$OutLogit[,2]))) == TRUE ||
    results_savev$OutLogit[1,1] < 0) & initcount < 20) {

initcount <- initcount + 1
    
initparams <- initparamssave[[initcount]]

results_savev <- discretefish_subroutine(otherdatfin$catchfin,
    otherdatfin$choicefin,otherdatfin$distance,otherdatfin,initparams,
    optimOpt,func,methodname)

}

newabund[i] <- sum(results_savev$OutLogit[2:(kk+1),1])*scaleq

browser()

}

dattemp$CPUE$obs <- unlist(newabund)
dattemp$CPUE$index <- 3

abundout <- merge(realcpue, dattemp$CPUE,  by = c("year", "seas"))
abundout$index.x <- NULL
abundout$se_log.x <- NULL
abundout$index.y <- NULL
names(abundout)[names(abundout) == "obs.x"] <- "TrueCPUE"
names(abundout)[names(abundout) == "obs.y"] <- "DahlCPUE"

abundout$diffperc = (abundout$TrueCPUE - abundout$DahlCPUE)/abundout$TrueCPUE

abundout <- abundout[order(abundout$year),]

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)
write.table(abundout, 
file=paste0("C:\\Users\\Allen.Chen\\SS3SIM_SCRATCH\\080721_mortup4\\",
    "abund_indices\\rel-se-abund-",gsub("/", "-", abundtitle),".csv"), 
sep=",", row.names=FALSE, quote = FALSE)

} else {

abundout <- read.table(paste0("C:\\Users\\Allen.Chen\\SS3SIM_SCRATCH\\", 
    "080721_mortup4\\", "abund_indices\\rel-se-abund-",
    gsub("/", "-", abundtitle),".csv"), sep=",", header=TRUE)
# dattemp$CPUE$se_log <- mean(abs(abundout$diffperc))
# dattemp$CPUE$se_log <- sqrt(log(1+((sd(dattemp$CPUE$obs)/
    #mean(dattemp$CPUE$obs))^2)))

abundout <- abundout[order(abundout$year),]

dattemp$CPUE$obs <- abundout$DahlCPUE
dattemp$CPUE$index <- 3

}

dattemp$CPUE$se_log <- 0.2

dat_list$CPUE <- rbind(dat_list$CPUE, dattemp$CPUE)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)

}
