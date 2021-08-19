#' moredatasample
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatasample <- function(datfile.origsave,dat_list,locnum,obsnum,betavar,
    uparams,abundse=NULL,abundscale=NULL,filename=NULL,...){
dattemp <- sample_index(dat_list        = datfile.origsave,
						outfile         = NULL,
						fleets          = 2,
						years           = list(74:100),
						sds_obs         = list(0.2),
						write_file      = FALSE)
                        
realcpue <- (datfile.origsave$CPUE[74:100,])

if (is.null(abundscale) == TRUE) {
    scaleq <- mean(realcpue$obs)/sum(betavar)
} else {
    scaleq <- abundscale
}

scalecatch <- realcpue$obs/mean(realcpue$obs)

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)

if (is.null(filename) == TRUE) {
    abundfile <- paste0(getwd(), "/", abundtitle, "/samp-abund-", 
        gsub("/", "-", abundtitle), ".csv")
} else {
    
    abundfile <- paste0(filename, gsub("/", "-", gsub("^[^/]*/","",abundtitle)),
    ".csv")
}

newabund <- list()
for (i in 1:length(scalecatch)) { 

library(barebones.FishSET)

betavarin <- as.matrix(betavar)*scalecatch[i]

alpha <- uparams$alpha
betac <- uparams$betac

distance <- locnum

kk <- dim(distance)[1]

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

choice <- list()
yikchosen <- list()
siout <- list()
for (j in 1:kk) {

yik <- list()

for (k in 1:kk) {

###################################Here for catch error
yik[[k]] <- betavarin[k,]*si[[j]] + bik[[j]][,k]

}

choice[[j]] <- matrix(as.numeric(sample(1:kk,ii[j]*1,replace=TRUE),ii[j],1))
yikchosen[[j]] <- as.matrix(diag(matrix(unlist(yik),ii[j],kk)[,choice[[j]]]))

siout[[j]] <- si[[j]]

}

sifin <- data.frame(V1 = as.numeric(unlist(siout)),
    V2 = as.numeric(unlist(siout)), V3 = as.numeric(unlist(siout)),
    V4 = as.numeric(unlist(siout)))
sifin <- matrix(as.numeric(unlist(siout)),sum(ii),kk)

catchfin <- data.frame(V1 = unlist(yikchosen))
choicefin <- data.frame(V1 = unlist(choice))

XX <- model.matrix(~as.factor(V1)-1, choicefin)*sifin
YY <- catchfin$V1

results_savev <- lm(YY~XX-1)

newabund[i] <- sum(results_savev$coef)*scaleq

}

dattemp$CPUE$obs <- unlist(newabund)
dattemp$CPUE[is.na(dattemp$CPUE$obs)==FALSE,]
dattemp$CPUE$index <- 3

abundout <- merge(realcpue, dattemp$CPUE,  by = c("year", "seas"))
abundout$index.x <- NULL
abundout$se_log.x <- NULL
abundout$index.y <- NULL
names(abundout)[names(abundout) == "obs.x"] <- "TrueCPUE"
names(abundout)[names(abundout) == "obs.y"] <- "DahlCPUE"

abundout$diffperc = (abundout$TrueCPUE - abundout$DahlCPUE)/abundout$TrueCPUE

abundout <- abundout[order(abundout$year),]

write.table(abundout, file=abundfile, sep=",", row.names=FALSE, quote = FALSE)

if (is.null(abundse) == TRUE) {
    dattemp$CPUE$se_log <- sqrt(log(1+((sd(dattemp$CPUE$obs)/
        mean(dattemp$CPUE$obs))^2)))
} else {
    dattemp$CPUE$se_log <- abundse
}

dat_list$CPUE <- rbind(dat_list$CPUE, dattemp$CPUE)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)

}
