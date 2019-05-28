#' moredatabiasreadrelSE
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatabiasreadrelSE <- function(datfile.origsave,dat_list){

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)
dahldirout <- "C:/Users/allen.chen/SS3SIM_SCRATCH/abund_indices/"

abundtitle2 <- gsub("/", "-", abundtitle)
abundtitle3 <- gsub(".*-", "", abundtitle2)

dattemp <- read.csv(paste(dahldirout,"newbiasabund-D225-E39-F2-M0-cod-",abundtitle3,".csv",sep=""))

#year seas index obs se_log

dattemp$index <- 3
dattemp$obs <- dattemp$BiasCPUE
# dattemp$se_log = 0.2
dattemp$se_log <- sqrt(log(1+((sd(dattemp$obs)/mean(dattemp$obs))^2)))

dattemp$sumcatches <- NULL
dattemp$TrueCPUE <- NULL
dattemp$BiasCPUE <- NULL
dattemp$se_log.y <- NULL
dattemp$diffperc <- NULL

# dattemp$obs <- dattemp$obs/1000

# dattemp$se_log <- 0.2

dat_list$CPUE <- rbind(dat_list$CPUE, dattemp)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)

}