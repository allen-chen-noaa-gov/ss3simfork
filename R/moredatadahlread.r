#' moredata
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatadahlread <- function(datfile.origsave,dat_list){

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)
dahldirout <- "C:/Users/allen.chen/Work/SS3SIM/dahldat/"

abundtitle2 <- gsub("/", "-", abundtitle)
abundtitle3 <- gsub(".*-", "", abundtitle2)

dattemp <- read.csv(paste(dahldirout,"dahl-D33-E33-F0-M0-cod-",abundtitle3,".csv",sep=""))

dattemp$se_log <- sqrt(log(1+((sd(dattemp$obs)/mean(dattemp$obs))^2)))

dat_list$CPUE <- rbind(dat_list$CPUE, dattemp)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)

}