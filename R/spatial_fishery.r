#' spatial_fishery
#'
#'
#' @param locnum description
#' @param obsnum description
#'
#' @export

spatial_fishery <- function(locnum,obsnum,betavar,uparams){

alpha <- uparams$alpha
betac <- uparams$betac

# betavar <- as.matrix(c(1.50, 1.25, 1.00, 0.75))*scalecatch[i]
# betavar <- as.matrix(betavar)*scalecatch[i]

##########CATCH DIST HERE##########
# distance <- rbind(t(as.matrix(c(0.0, 0.5, 0.5, 0.707))), 
			# t(as.matrix(c(0.5, 0.0, 0.707, 0.5))), 
			# t(as.matrix(c(0.5, 0.707, 0, 0.5))), 
			# t(as.matrix(c(0.707, 0.5, 0.5, 0))))
# distance <- distance*3

distance <- locnum

# kk <- 4
kk <- dim(distance)[1]

# ii <- rep(500, kk)
if (length(obsnum) == 1) {
    ii <- rep(obsnum, kk)
} else {
    ii <- obsnum
}

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

choice[[j]] <- as.matrix(which(t(apply(matrix(unlist(Vijk),ii[j],kk),1,max) == 
    matrix(unlist(Vijk),ii[j],kk)))-(((1:ii[j])-1)*kk))
yikchosen[[j]] <- as.matrix(diag(matrix(unlist(yik),ii[j],kk)[,choice[[j]]]))

siout[[j]] <- si[[j]]
ziout[[j]] <- zi[[j]]
startlocout[[j]] <- rep(j,ii[j])

distanceout[[j]] <- t(matrix(rep(distance[j,],ii[j]),kk,ii[j]))

}

zifin <- data.frame(V1 = as.numeric(unlist(ziout)))
startlocfin <- data.frame(V1 = as.numeric(unlist(startlocout)))

sifin <- data.frame(V1 = as.numeric(unlist(siout)), 
    V2 = as.numeric(unlist(siout)), V3 = as.numeric(unlist(siout)), 
    V4 = as.numeric(unlist(siout)))
sifin <- matrix(as.numeric(unlist(siout)),sum(ii),kk)

catchfin <- data.frame(V1 = unlist(yikchosen))
choicefin <- data.frame(V1 = unlist(choice))

distancefin <- data.frame(do.call(rbind,distanceout))
colnames(distancefin) <- c("V1","V2","V3","V4")
				
###################################Here for multiple params

intdatfin <- list(zi=zifin)

griddatfin <- list(si=sifin)

startlocdatfin <- list(startloc=startlocfin)

otherdatfin <- list(griddat=list(as.matrix(sifin)), noCgriddat = NA,
    intdat=list(as.matrix(zifin)), startloc=as.matrix(startlocfin),
    catchfin=catchfin, choicefin=choicefin)
otherdatfin$distance <- as.matrix(distancefin)

return(otherdatfin)

}
