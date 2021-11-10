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
    abundfile <- paste0(getwd(), "/", abundtitle, "/samp-abund-", 
        gsub("/", "-", abundtitle), ".csv")
} else {
    
    abundfile <- paste0(filename, gsub("/", "-", gsub("^[^/]*/","",abundtitle)),
    ".csv")
}

newabund <- list()
for (i in 1:length(scalecatch)) { 
set.seed(i)

betavarin <- as.matrix(betavar)*scalecatch[i]

alpha <- uparams$alpha
betac <- uparams$betac

distance <- locnum

kk <- dim(distance)[1]

otherdatfin <- spatial_fishery(locnum,obsnum,betavarin,uparams,
    datfile.origsave$catch[(100-list(...)$obsyears):100,],year=i,random=TRUE,
    list(...)$avghauls)

choicefin <- otherdatfin$choicefin
sifin <- do.call(cbind,otherdatfin$griddat)
catchfin <- otherdatfin$catchfin

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
names(abundout)[names(abundout) == "obs.y"] <- "EstCPUE"

abundout$diffperc = (abundout$TrueCPUE - abundout$EstCPUE)/abundout$TrueCPUE

abundout <- abundout[order(abundout$year),]

write.table(abundout, file=abundfile, sep=",", row.names=FALSE, quote = FALSE)

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
