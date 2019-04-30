moredatadahlread <- function(datfile.origsave,dat_list){

abundtitle <- sub("/\\s*em\\b.*", "", dat_list$`sourcefile`)
dahldirout <- "C:/Users/allen.chen/Work/SS3SIM/dahldat/"
dattemp <- read.csv(paste(dahldirout,"dahl-",gsub("/", "-", abundtitle),".csv",sep=""))

dat_list$CPUE <- rbind(dat_list$CPUE, dattemp)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)

}