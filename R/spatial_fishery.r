#' spatial_fishery
#'
#'
#' @param locnum description
#' @param obsnum description
#'
#' @export

spatial_fishery <- function(locnum,obsnum,betavar,uparams,totcatches,year,
    random=FALSE,avghauls,catchscale){

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

betavar <- betavar
distance <- locnum

# kk <- 4
kk <- dim(distance)[1]

# ii <- rep(500, kk)
if (length(obsnum) == 1) {
    ii <- rep(obsnum, kk)
} else {
    ii <- obsnum
}

ii <- ii/sum(ii)

yikchosen <- list()
choice <- list()
siout <- list()
ziout <- list()
startlocout <- list()
distanceout <- list()
triplength <- list()
counter <- 1
haulcounter <- 1
#this is calibrated for about 3000 hauls in a year
while (sum(unlist(yikchosen))*(catchscale) < 
    totcatches$Fishery[year]) {

si <- sample(1:5,1)
zi <- sample(1:10,1)

bik <- rnorm(kk,0,3)

wijk <- -log(rexp(kk,1))
#-log(exp(1)) is standard type 1 extreme value i.e. gumbel beta=1 mu=0

if (haulcounter == 1) {
j <- sample(1:kk, size = 1, prob = ii)
}

yik <- list()
Vijk <- list()
for (k in 1:kk) {

###################################Here for multiple params
tijk <- betac*distance[j,k]*zi + wijk[k]

###################################Here for catch error
yik[[k]] <- betavar[k,]*si + bik[k]

###################################Here for multiple params
Vijk[[k]] <- alpha*yik[[k]] + tijk

}

if (random == TRUE) {
choice[[counter]] <- sample(1:kk,1)
} else {
choice[[counter]] <- which(max(unlist(Vijk)) == Vijk)
}

yikchosen[[counter]] <- yik[[choice[[counter]]]]

siout[[counter]] <- si
ziout[[counter]] <- zi
startlocout[[counter]] <- j

distanceout[[counter]] <- distance[j,]

if (rpois(1, haulcounter) >= avghauls) {
triplength[[counter]] <- haulcounter
haulcounter <- 1 
} else {
j <- choice[[counter]]
haulcounter <- haulcounter + 1
}

counter <- counter+1

}

zifin <- data.frame(V1 = as.numeric(unlist(ziout)))
startlocfin <- data.frame(V1 = as.numeric(unlist(startlocout)))

# sifin <- data.frame(V1 = as.numeric(unlist(siout)), 
    # V2 = as.numeric(unlist(siout)), V3 = as.numeric(unlist(siout)), 
    # V4 = as.numeric(unlist(siout)))
sifin <- matrix(as.numeric(unlist(siout)), length(unlist(siout)),max(unlist(choice)))

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
