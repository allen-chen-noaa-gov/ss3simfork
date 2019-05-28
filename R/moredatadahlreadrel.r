#' moredata
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatadahlreadrel <- function(datfile.origsave,dat_list){

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)
dahldirout <- "C:/Users/Allen/Desktop/abund_indices/dahldat/flatalpha6/"

abundtitle2 <- gsub("/", "-", abundtitle)
abundtitle3 <- gsub(".*-", "", abundtitle2)

dattemp <- read.csv(paste(dahldirout,"dahl-D225-E39-F2-M0-cod-",abundtitle3,".csv",sep=""))

dattemp$sumcatches <- NULL

# dattemp$obs <- dattemp$obs/1000

# dattemp$se_log <- sqrt(log(1+((sd(dattemp$obs)/mean(dattemp$obs))^2)))
dattemp$se_log <- 0.2

dat_list$CPUE <- rbind(dat_list$CPUE, dattemp)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)

}