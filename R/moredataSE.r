#' moredataSE
#'
#'
#' @param datfile.origsave description
#' @param dat_list description
#'
#' @export

moredataSE <- function(datfile.origsave,dat_list){
dattemp <- sample_index(dat_list        = datfile.origsave,
						outfile         = NULL,
						fleets          = 2,
						years           = list(74:100),
						sds_obs         = list(0.2),
						write_file      = FALSE)
		
dattemp$CPUE$index <- 3
newabund <- dattemp$CPUE

newabund$se_log <- sqrt(log(1+((sd(newabund$obs)/mean(newabund$obs))^2)))
# newabund$se_log <- 0.2

dat_list$CPUE <- rbind(dat_list$CPUE, newabund)

rownames(dat_list$CPUE) <- seq(length=nrow(dat_list$CPUE))

dat_list$N_cpue <- dim(dat_list$CPUE)[1]
dat_list$NCPUEObs[3] <- dim(dat_list$CPUE[dat_list$CPUE$index==3,])[1]

return(dat_list)
}