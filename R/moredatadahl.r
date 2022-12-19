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
    uparams,abundse=NULL,catchscale=NULL,filename=NULL,...){
dattemp <- sample_index(dat_list        = datfile.origsave,
						outfile         = NULL,
						fleets          = 2,
						years           = list((100-list(...)$obsyears):100),
						sds_obs         = list(0.2),
						write_file      = FALSE)

realcpue <- (datfile.origsave$CPUE[(100-list(...)$obsyears):100,])

if (is.null(catchscale) == TRUE) {
    catchscale <- min(datfile.origsave$catch[(100-list(...)$obsyears):100,])/1000
}

scaleabund <- realcpue$obs/mean(realcpue$obs)

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)

if (is.null(filename) == TRUE) {
    abundfile <- paste0(getwd(), "/", abundtitle, "/cor-abund-", 
        gsub("/", "-", abundtitle), ".csv")
} else {
    abundfile <- paste0(filename,  gsub(".*/","",abundtitle), ".csv")
}

if (file.exists(abundfile) == FALSE) {
    
newabund <- list()
newse <- list()
choiceout <- list()
startout <- list()
paramsout <- list()
trueout <- list()
for (i in 1:length(scaleabund)) { 

betavarscaled <- (betavar/sum(betavar))*10.125

if (list(...)$trend==TRUE) {
slope <- (max(betavar)-min(betavar))/list(...)$obsyears
slopesign <- sign(betavar[1]-betavar[length(betavar)])
betaslope <- seq((-slopesign)*slope, slopesign*slope,by=(
    slopesign*(2*slope)/(length(betavar)-1)))
betavarscaled <- ((betavar/sum(betavar))*10.125)+(betaslope*i)  
}

betavarin <- as.matrix(betavarscaled)*scaleabund[i]

otherdatfin <- spatial_fishery(locnum,obsnum,betavarin,uparams,
    datfile.origsave$catch[(100-list(...)$obsyears):100,],year=i,random=FALSE,
    list(...)$avghauls,catchscale,list(...)$catchvarV,list(...)$catchvarN)
otherdatfinsave <- otherdatfin

polyn <- 2
polyintnum <- 1
regconstant <- 0
polyconstant <- 1
singlecor <- 1

otherdatfin$polyn <- polyn
otherdatfin$polyintnum <- polyintnum
otherdatfin$regconstant <- regconstant
otherdatfin$polyconstant <- polyconstant
otherdatfin$singlecor <- singlecor
zifin <- do.call(cbind,otherdatfin$intdat)

# for not estimating areas far away with few observations
set <- as.numeric(unlist(dimnames(table(otherdatfin$choicefin)[
    table(otherdatfin$choicefin) >= list(...)$minobs])))
while(all(table(otherdatfin$choicefin)>=list(...)$minobs) == FALSE) {
set <- as.numeric(unlist(dimnames(table(otherdatfin$choicefin)[
    table(otherdatfin$choicefin) >= list(...)$minobs])))
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

otherdatfin$choicefin$V1 <- as.factor(otherdatfin$choicefin$V1)
levels(otherdatfin$choicefin$V1) <- 
    1:length(unique(levels(as.factor(otherdatfin$choicefin$V1))))
otherdatfin$choicefin$V1 <- as.integer(otherdatfin$choicefin$V1)
# end of areas far away

kk <- length(betavarin[set])
if (regconstant == 1) {
    initparams <- unname(c(1, betavarin[set], rep(0, 
        (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
        rep(-1,dim(zifin)[2]), 1))
    #Initial paramters for revenue then cost.
} else {
    initparams <- unname(c(1, betavarin[set], rep(0, 
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
changevec <- unname(c(1, rep(0, length(betavarin[set])),
    rep(1, (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
    rep(1,dim(zifin)[2]), 1)) 
    # Initial paramters for revenue then cost.
} else {
changevec <- unname(c(1, rep(0, length(betavarin[set])),
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
    results_savev$OutLogit[1,1] < 0) & initcount < 10) {

initcount <- initcount + 1
    
initparams <- initparamssave[[initcount]]

results_savev <- barebones.FishSET::discretefish_subroutine(otherdatfin$catchfin,
    otherdatfin$choicefin,otherdatfin$distance,otherdatfin,initparams,
    optimOpt,func,methodname)

}

fin <- c(catchscale)
bstring <- paste0("~(", paste(paste("x", 2:(kk+1), sep=""), collapse="+"), 
	")*(%f)")
form <- do.call(sprintf, c(fmt = bstring, as.list(fin)))
seout <- msm::deltamethod(as.formula(form), results_savev$OutLogit[,1], 
	results_savev$H1)

if (is.numeric(results_savev$OutLogit[2:(kk+1),1]) == FALSE ||
    any(is.na(as.numeric(results_savev$OutLogit[,2]))) == TRUE) {
    newabund[i] <- NA
	newse[i] <- NA
} else {
    newabund[i] <- sum(results_savev$OutLogit[2:(kk+1),1])*catchscale
	newse[i] <- sqrt(log(1+(((seout)/(newabund[[i]]))^2)))
}
paramsout[[i]] <-  results_savev$OutLogit[2:(kk+1),1]
trueout[[i]] <- t(betavarin)
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

abundout$diffperc = (abundout$TrueCPUE - abundout$EstCPUE)/abundout$TrueCPUE

abundout <- abundout[order(abundout$year),]

write.table(abundout, file=abundfile, sep=",", row.names=FALSE, quote = FALSE)
write.table(cbind(do.call(rbind, startout), do.call(rbind, choiceout)), 
	file=paste0(getwd(), "/", abundtitle, "/sampout-", 
	gsub("/", "-", abundtitle), ".csv"), 
	sep=",", row.names=FALSE, quote = FALSE)
write.table(cbind(do.call(rbind, trueout), do.call(rbind, paramsout)), 
	file=paste0(getwd(), "/", abundtitle, "/paramsout-", 
	gsub("/", "-", abundtitle), ".csv"), 
	sep=",", row.names=FALSE, quote = FALSE)
    
} else {

abundout <- read.table(abundfile, sep=",", header=TRUE)

abundout <- abundout[order(abundout$year),]

# dattemp$CPUE$obs <- abundout$EstCPUE
# dattemp$CPUE$index <- 3
# dattemp$CPUE <- dattemp$CPUE[is.na(dattemp$CPUE$obs) == FALSE, ]


}

# if (is.null(abundse) == TRUE) {
    # dattemp$CPUE$se_log <- sqrt(log(1+((sd(dattemp$CPUE$obs)/
        # mean(dattemp$CPUE$obs))^2)))
# } else {
    # dattemp$CPUE$se_log <- abundse
# }

dat_list$CPUE <- rbind(dat_list$CPUE, 
    dattemp$CPUE[is.na(dattemp$CPUE$obs)==FALSE,])

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)

}
