#' moredatatest
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredatatest <- function(datfile.origsave,dat_list){
dattemp <- sample_index(dat_list        = datfile.origsave,
						outfile         = NULL,
						fleets          = 3,
						years           = list(74:100),
						sds_obs         = list(0.2),
						write_file      = FALSE)
		
newabund <- dattemp$CPUE
newabund$obs <- 100000000
newabund$se_log <- 0.2

dat_list$CPUE <- rbind(dat_list$CPUE, newabund)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)
}