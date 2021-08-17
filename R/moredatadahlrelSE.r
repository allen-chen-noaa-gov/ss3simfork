#' moredatadahlrelSE
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatadahlrelSE <- function(datfile.origsave,dat_list){
dattemp <- sample_index(dat_list        = datfile.origsave,
						outfile         = NULL,
						fleets          = 2,
						years           = list(74:100),
						sds_obs         = list(0.2),
						write_file      = FALSE)
						
realcpue <- (datfile.origsave$CPUE[74:100,])

scaleq <- mean(realcpue$obs)/4.5
		
scalecatch <- realcpue$obs/mean(realcpue$obs)

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)

if (file.exists(paste0("C:\\Users\\Allen.Chen\\SS3SIM_SCRATCH\\080721_mortup4\\",
    "abund_indices\\rel-se-abund-",gsub("/", "-", abundtitle),".csv")) == 
    FALSE) {

newabund <- list()
for (i in 1:length(scalecatch)) { 

library(barebones.FishSET)

kk <- 4
ii <- rep(500, kk)

alpha <- 3
betac <- -(1)

betavar <- as.matrix(c(1.50, 1.25, 1.00, 0.75))*scalecatch[i]

##########CATCH DIST HERE##########
distance <- rbind(t(as.matrix(c(0.0, 0.5, 0.5, 0.707))), 
			t(as.matrix(c(0.5, 0.0, 0.707, 0.5))), 
			t(as.matrix(c(0.5, 0.707, 0, 0.5))), 
			t(as.matrix(c(0.707, 0.5, 0.5, 0))))

distance <- distance*3

si <- list()
zi <- list()
for (l in 1:kk) {
si[[l]] <- matrix(sample(1:5,ii[l]*1,replace=TRUE),ii[l],1)
zi[[l]] <- matrix(sample(1:10,ii[l]*1,replace=TRUE),ii[l],1)
}

bik <- list()	
for (l in 1:kk) { #can vectorize?
bik[[l]] <- matrix(rnorm(ii[l]*kk,0,3),ii[l],kk)
}

wijk <- list()
for (l in 1:kk) { 
wijk[[l]] <- matrix((-log(rexp(ii[l]*kk,1))),ii[l],kk) 
#-log(exp(1)) is standard type 1 extreme value i.e. gumbel beta=1 mu=0
}

choice <- list()
yikchosen <- list()
siout <- list()
siout2 <- list()
ziout <- list()
ziout2 <- list()
startlocout <- list()
distanceout <- list()
predyik <- list()
predyik2 <- list()
for (j in 1:kk) {

Vijk <- list()
yik <- list()

for (k in 1:kk) {

###################################Here for multiple params
tijk <- betac*distance[j,k]*zi[[j]] + wijk[[j]][,k]

###################################Here for catch error
yik[[k]] <- betavar[k,]*si[[j]] + bik[[j]][,k]

###################################Here for multiple params
Vijk[[k]] <- alpha*yik[[k]] + tijk

}

choice[[j]] <- as.matrix(which(t(apply(matrix(unlist(Vijk),ii[j],kk),1,max)==matrix(unlist(Vijk),ii[j],kk)))-(((1:ii[j])-1)*kk))
yikchosen[[j]] <- as.matrix(diag(matrix(unlist(yik),ii[j],kk)[,choice[[j]]]))

siout[[j]] <- si[[j]]
ziout[[j]] <- zi[[j]]
startlocout[[j]] <- rep(j,ii[j])

distanceout[[j]] <- t(matrix(rep(distance[j,],ii[j]),kk,ii[j]))

}

zifin <- data.frame(V1 = as.numeric(unlist(ziout)))
startlocfin <- data.frame(V1 = as.numeric(unlist(startlocout)))

sifin <- data.frame(V1 = as.numeric(unlist(siout)),V2 = as.numeric(unlist(siout)),V3 = as.numeric(unlist(siout)),V4 = as.numeric(unlist(siout)))
sifin <- matrix(as.numeric(unlist(siout)),sum(ii),kk)

catchfin <- data.frame(V1 = unlist(yikchosen))
choicefin <- data.frame(V1 = unlist(choice))

distancefin <- data.frame(do.call(rbind,distanceout))
colnames(distancefin) <- c("V1","V2","V3","V4")
				
###################################Here for multiple params

intdatfin <- list(zi=zifin)

griddatfin <- list(si=sifin)

startlocdatfin <- list(startloc=startlocfin)

polyn <- 3
polyintnum <- 1
regconstant <- 0
polyconstant <- 1
singlecor <- 0

otherdatfin <- list(griddat=list(as.matrix(sifin)), noCgriddat = NA,
    intdat=list(as.matrix(zifin)), startloc=as.matrix(startlocfin), 
    polyn=polyn, polyintnum = polyintnum, regconstant = regconstant, 
    polyconstant = polyconstant, singlecor = singlecor)

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

optimOpt <- c(100000,1.00000000000000e-08,1,0) #Optimization options for the maximum number of
					   #function evaluations, maximum iterations, and the relative tolerance of x.
					   #Then, how often to report output, and whether to report output.
methodname = "BFGS"

bw <- -1

otherdatfin$distance <- as.matrix(distancefin)
otherdatfin$bw <- bw

func <- logit_correction_polyint_estscale

results_savev <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)
    
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
    logit_correction_polyint, catchfin, choicefin, distancefin, otherdatfin,
    changevec)

LLmat <- data.frame(cbind(1:searchspace,unlist(results$saveLLstarts)))
LLmatorder <- LLmat[order(LLmat$X2),]

initparamssave <- results$savestarts[LLmatorder$X1[1:100]]

while ((any(is.na(as.numeric(results_savev$OutLogit[,2]))) == TRUE ||
    results_savev$OutLogit[1,1] < 0) & initcount < 20) {

initcount <- initcount + 1
    
initparams <- initparamssave[[initcount]]


results_savev <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)

}

newabund[i] <- sum(results_savev$OutLogit[2:5,1])*1000000

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
# dattemp$CPUE$se_log <- sqrt(log(1+((sd(dattemp$CPUE$obs)/mean(dattemp$CPUE$obs))^2)))

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
