#' moredatasamplerel
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatasamplerel <- function(datfile.origsave,dat_list){
dattemp <- sample_index(dat_list        = datfile.origsave,
						outfile         = NULL,
						fleets          = 2,
						years           = list(74:100),
						sds_obs         = list(0.2),
						write_file      = FALSE)
						
realcpue <- (datfile.origsave$CPUE[74:100,])

scaleq <- mean(realcpue$obs)/6.75
		
scalecatch <- realcpue$obs/mean(realcpue$obs)

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

newabund <- list()
for (i in 1:length(scalecatch)) { 

library(barebones.FishSET)
library(MASS)

kk <- 6
# ii <- rep(500, kk)
ii <- as.numeric(table(sample(1:kk,2000,replace=TRUE)))

alpha <- 3
betac <- -(1)

# betavar <- as.matrix(c(1.50, 1.25, 1.00, 0.75))*scalecatch[i]
			
distance = matrix(scan("C:/Users/Allen/Dropbox/Allen/Work/NMFS/packages_docs/testcodefor_barebones.FishSET/dahlMCs/sixdistances.csv", sep=",",quiet=TRUE),kk,kk)

distance <- distance*3

sigmaspace <- 1

sspace <- 5

covarspace <- (sigmaspace^2)*exp(-(distance/sspace)^2)

checkzero <- 1

while (checkzero>0 || is.na(checkzero)==TRUE){
Xspace <- mvrnorm(1,mu=rep(0,kk),Sigma=covarspace)
deltaspace <- mvrnorm(1,mu=rep(0,kk),Sigma=covarspace)

betanotspace <- 2

betaonespace <- runif(kk, -0.5, 0.5)

# exp(betanotspace + betaonespace*Xspace)

# saveabund[[i]] <- rpois(4,exp(betanotspace + betaonespace*Xspace))
betavar <- rpois(kk,exp(betanotspace + betaonespace*Xspace + deltaspace))
betavar <- betavar/(sum(betavar)/6.75)

checkzero <- (min(betavar)<0.5)
}

betavar <- as.matrix(betavar)*scalecatch[i]

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

# choice[[j]] <- as.matrix(which(t(apply(matrix(unlist(Vijk),ii[j],kk),1,max)==matrix(unlist(Vijk),ii[j],kk)))-(((1:ii[j])-1)*kk))
choice[[j]] <- matrix(as.numeric(sample(1:kk,ii[j]*1,replace=TRUE),ii[j],1))
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

newabund[i] <- sum(results_savev$coef)*1000000
}

dattemp$CPUE$obs <- unlist(newabund)
dattemp$CPUE$index <- 3

abundout <- merge(realcpue, dattemp$CPUE,  by = c("year", "seas"))
abundout$index.x <- NULL
abundout$se_log.x <- NULL
abundout$index.y <- NULL
names(abundout)[names(abundout) == "obs.x"] <- "TrueCPUE"
names(abundout)[names(abundout) == "obs.y"] <- "BiasCPUE"

abundout$diffperc = (abundout$TrueCPUE - abundout$BiasCPUE)/abundout$TrueCPUE

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)
write.table(abundout, 
file=paste0("C:\\Users\\Allen\\Desktop\\abund_indices\\flatalpha6\\sampleabund-",gsub("/", "-", abundtitle),".csv"), 
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