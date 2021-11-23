rm(list=ls())
gc()
options(scipen = 99)

#for running the model
library(ss3simfork)
library(doParallel)
library(foreach)

#for making the plots
library(plyr)
library(ggplot2)
library(RColorBrewer)

d <- system.file("extdata", package = "ss3simfork")
om <- paste0(d, "/models/cod-om")
em <- paste0(d, "/models/cod-em")

#rawoutputdir for files (iterations and scenarios)
rawoutputdir <- #dir

#nicedir for cases
nicedir <- #dir
case_folder <- paste0(nicedir, "\\eg-cases")

###############################################################################

setwd(rawoutputdir) #This is directory where RAW output get written

registerDoParallel(cores = 8) 
#The current directory when you register parallel cores is the directory output 
#gets written

print(Sys.time())
begtime <- Sys.time()

run_ss3sim(iterations = 1:100, scenarios =
     c(
    "D4-E1-F2-M0-cod",
    "D1-E1-F2-M0-cod", 
    "D2-E1-F2-M0-cod", 
    "D3-E1-F2-M0-cod"
     ),
  case_files = list(F = "F", D = c("index", "lcomp", "agecomp", "econ"),
    E = "E", M = "M"),
  case_folder = case_folder, om_dir = om,
  em_dir = em, bias_adjust = FALSE, parallel = TRUE, parallel_iterations = TRUE)
#I thought I remembered something about a problem with bias_adjust?
  
endtime <- Sys.time()
print(Sys.time())

endtime - begtime

###############################################################################

get_results_all(directory = rawoutputdir, 
    overwrite_files = TRUE, 						
    #this writes to rawoutputdir even if you set a nicedir before call
    user_scenarios = 
    c(    
    "D1-E1-F2-M0-cod", 
    "D2-E1-F2-M0-cod", 
    "D3-E1-F2-M0-cod", 
    "D4-E1-F2-M0-cod"
    ), 
    parallel = TRUE)
		
scalar_dat <- read.csv("ss3sim_scalar.csv")
ts_dat <- read.csv("ss3sim_ts.csv")

ts_dat$SpawnBio <- (ts_dat$SpawnBio_em - ts_dat$SpawnBio_om)/
    (ts_dat$SpawnBio_om)
  
ts_dat <- merge(ts_dat, scalar_dat[,c("scenario", "replicate",
    "max_grad")])

subsetlist <- c("D1","D2","D3","D4")

ts_dat_sto <- subset(ts_dat, D %in% subsetlist)

ts_dat_sto$D <- as.factor(ts_dat_sto$D)

ts_dat_sto$D <- revalue(ts_dat_sto$D, c(
    "D1"="Quadrennial fishery-independent survey",
    "D2"="+ Annual randomly-sampled fishery data",
    "D3"="+ Annual uncorrected fishery data",
    "D4"="+ Annual corrected fishery data"))
    
#Small values of the maximum gradient (approximately 0.001 or less) indicate 
#that model convergence is likely. 
#Larger values (greater than 1) indicate that model convergence is unlikely. 
ts_dat_sto <- ts_dat_sto[ts_dat_sto$max_grad<1,] 

test.50 <- ddply(ts_dat_sto, .(D,year), summarise, 
    five0 = quantile(SpawnBio, probs = c(0.5), na.rm = TRUE))

test.05 <- ddply(ts_dat_sto, .(D,year), summarise, 
    zero5 = quantile(SpawnBio, probs = c(0.025), na.rm = TRUE))
test.95 <- ddply(ts_dat_sto, .(D,year), summarise, 
    nine5 = quantile(SpawnBio, probs = c(0.975), na.rm = TRUE))
test.25 <- ddply(ts_dat_sto, .(D,year), summarise, 
    two5 = quantile(SpawnBio, probs = c(0.25), na.rm = TRUE))
test.75 <- ddply(ts_dat_sto, .(D,year), summarise, 
    seven5 = quantile(SpawnBio, probs = c(0.75), na.rm = TRUE))    
test.05$nine5 <- test.95$nine5
test.05$two5 <- test.25$two5
test.05$seven5 <- test.75$seven5

p <- ggplot(test.05, aes(x = year)) + xlab("Year") +
    theme_bw() + 
	geom_ribbon(data = test.05, aes(ymin = zero5, ymax = nine5), alpha = 0.50, 
        fill = "#4682b4") +
	geom_ribbon(data = test.05, aes(ymin = two5, ymax = seven5), alpha = 0.75, 
        fill = "#4682b4") +
	theme(text = element_text(size=30), axis.title.y=element_text(vjust=1.5)) + 
	geom_line(data=test.50, aes(y=five0, x = year), colour="black") +
	ylim(-0.55,0.55) + ylab("Relative Error in Spawning Biomass") + 
    xlab("Year") +
	# theme(axis.text.x = element_text(angle=45, hjust=0.1, vjust=0.1)) + 
	facet_wrap(~D, ncol=2)

ggsave(p, file=paste0(nicedir, "\\SB", ".png"), 
    width = 16, height = 9)

#############################################################################

abundout <- list()
for (it in 1:100) {

#The default naming convention for each abundance index .csv is below, if 
#the file name is left NULL, but can be specified if desired.

casename <- paste0(rawoutputdir, "\\D3-E1-F2-M0-cod\\", it, 
    "\\bias-abund-D3-E1-F2-M0-cod-")
abundtrue <- readabund(filename = paste(casename,it,".csv",sep=""),
    type="True", istruecpue = TRUE)

casename <- paste0(rawoutputdir, "\\D3-E1-F2-M0-cod\\", it, 
    "\\bias-abund-D3-E1-F2-M0-cod-")
abundbias <- readabund(filename = paste(casename,it,".csv",sep=""), 
    type="Uncorrected")

casename <- paste0(rawoutputdir, "\\D4-E1-F2-M0-cod\\", it, 
    "\\cor-abund-D4-E1-F2-M0-cod-")
abundcor <- readabund(filename = paste(casename,it,".csv",sep=""), 
    type="Corrected")

casename <- paste0(rawoutputdir, "\\D2-E1-F2-M0-cod\\", it, 
    "\\samp-abund-D2-E1-F2-M0-cod-")    
abundsamp <- readabund(filename = paste(casename,it,".csv",sep=""), 
    type="Randomly sampled")
    
tempabundout <- rbind(abundtrue, abundbias, abundcor, abundsamp)

tempabundout$iter <- it

abundout[[it]] <- tempabundout

}

totabundout <- do.call(rbind, abundout)

totabundout <- totabundout[totabundout$CPUE > 0, ]

linegraph = aggregate(totabundout$RelDiff, 
    list(Year = totabundout$year, Index = totabundout$Index), median)
names(linegraph)[names(linegraph) == 'x'] = 'RelDiff'

linegraphse = aggregate(totabundout$RelDiff, 
    list(Year = totabundout$year, Index = totabundout$Index), sd)
names(linegraphse)[names(linegraphse) == 'x'] = 'sd'

linegraphlen = aggregate(totabundout$RelDiff, 
    list(Year = totabundout$year, Index = totabundout$Index), length)
names(linegraphlen)[names(linegraphlen) == 'x'] = 'len'

linegraph = merge(linegraph, linegraphse, by=c("Year","Index"))
linegraph = merge(linegraph, linegraphlen, by=c("Year","Index"))

linegraph$se = linegraph$sd/sqrt(linegraph$len)

linegraph$Index <- revalue(linegraph$Index , c("Corrected"="Corrected     ",
    "Randomly sampled"="Randomly sampled     ",
    "True"="True     ",
    "Uncorrected"="Uncorrected     "))

pd <- position_dodge(0.05)

lineout = ggplot(data = linegraph, 
		aes(x=Year, y=RelDiff)) + 
        geom_line(aes(colour=Index, group=Index),size=2, position=pd) +
		geom_point(aes(shape=Index), size=4, fill="white", stroke = 1.5, 
            position=pd) + 
		xlab("Year") + ylab("Relative Differences in Abundance") +
		geom_hline(yintercept=1, linetype="dashed", color = "orange", size=2) +
		guides(colour = guide_legend(override.aes = list(size=6))) +
		theme(legend.direction = 'horizontal', legend.position = 'top', 
            legend.key = element_rect(size = 4), 
            legend.text = element_text(size=22),
            legend.key.size = unit(2.5, 'lines'), legend.title=element_blank(), 
            axis.text.x = element_text(angle = 45, hjust = 1),
            text = element_text(size=30), 
            axis.title.y=element_text(vjust=1.5)) +
		scale_shape_manual(values=c(21,22,25,24)) + 
        scale_colour_manual(values=brewer.pal(4,"Set1"))
        
ggsave(lineout, file=paste0(nicedir, "\\abundline", ".png"), 
    width = 16, height = 9)
    