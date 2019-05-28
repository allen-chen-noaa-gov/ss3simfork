#' moredatabiasrel
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatabiasrel <- function(datfile.origsave,dat_list){
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

newabund <- list()
sumcatches <- list()
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
bikreal <- list()	
for (l in 1:kk) { #can vectorize?
bik[[l]] <- matrix(rnorm(ii[l]*kk,0,3),ii[l],kk)
bikreal[[l]] <- matrix(rnorm(ii[l]*kk,0,1),ii[l],kk)
}

wijk <- list()
for (l in 1:kk) { 
wijk[[l]] <- matrix((-log(rexp(ii[l]*kk,1))),ii[l],kk) 
#-log(exp(1)) is standard type 1 extreme value i.e. gumbel beta=1 mu=0
}

choice <- list()
yikchosen <- list()
yikreal <- list()
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

yikreal[[k]] <- betavar[k,]*si[[j]] + bik[[j]][,k] + bikreal[[j]][,k]

###################################Here for multiple params
Vijk[[k]] <- alpha*yik[[k]] + tijk

}

choice[[j]] <- as.matrix(which(t(apply(matrix(unlist(Vijk),ii[j],kk),1,max)==matrix(unlist(Vijk),ii[j],kk)))-(((1:ii[j])-1)*kk))
yikchosen[[j]] <- as.matrix(diag(matrix(unlist(yikreal),ii[j],kk)[,choice[[j]]]))

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

otherdatfin <- list(griddat=as.matrix(sifin),intdat=as.matrix(zifin),startloc=as.matrix(startlocfin),
				polyn=polyn,gridnum=1,intnum=1)

###################################Here for multiple params
initparams <- c(3, betavar, rep(0, (((polyn+1)*2) + 2)*kk), 3, 3) #Initial paramters for revenue then cost.

optimOpt <- c(100000,1.00000000000000e-08,1,0) #Optimization options for the maximum number of
					   #function evaluations, maximum iterations, and the relative tolerance of x.
					   #Then, how often to report output, and whether to report output.
methodname = "BFGS"

func <- logit_correction_v
otherdatfin$distance <- as.matrix(distancefin)

XX <- model.matrix(~as.factor(V1)-1, choicefin)*sifin
YY <- catchfin$V1

results_savev <- lm(YY~XX-1)

# results_savev$catches = cbind(sifin,catchfin)

newabund[i] <- sum(results_savev$coef)*1000000 #was scaleq
sumcatches[i] <- sum(catchfin)
}

dattemp$CPUE$obs <- unlist(newabund)
dattemp$CPUE$index <- 3

abundout <- merge(realcpue, dattemp$CPUE,  by = c("year", "seas"))

abundout <- abundout[order(abundout$year),]

abundout$index.x <- NULL
abundout$se_log.x <- NULL
abundout$index.y <- NULL
names(abundout)[names(abundout) == "obs.x"] <- "TrueCPUE"
names(abundout)[names(abundout) == "obs.y"] <- "BiasCPUE"

abundout$diffperc = (abundout$TrueCPUE - abundout$BiasCPUE)/abundout$TrueCPUE

abundout$sumcatches <- unlist(sumcatches)

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)
write.table(abundout, 
file=paste0("C:\\Users\\allen.chen\\SS3SIM_SCRATCH\\abund_indices\\flat_makebias\\biasabund-",gsub("/", "-", abundtitle),".csv"), 
sep=",", row.names=FALSE, quote = FALSE)

# dattemp$CPUE$se_log <- mean(abs(abundout$diffperc))
# dattemp$CPUE$se_log <- sqrt(log(1+((sd(dattemp$CPUE$obs)/mean(dattemp$CPUE$obs))^2)))
dattemp$CPUE$se_log <- 0.2

dat_list$CPUE <- rbind(dat_list$CPUE, dattemp$CPUE)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)
}