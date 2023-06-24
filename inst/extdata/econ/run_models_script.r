rm(list = ls())
gc()
options(scipen = 99)

# for running the model
library(ss3simfork)
library(doParallel)
library(foreach)

# for making the plots
library(plyr)
library(ggplot2)
library(RColorBrewer)

d <- system.file("extdata", package = "ss3simfork")
om <- paste0(d, "/models/cod-om")
em <- paste0(d, "/models/cod-em")

# rawoutputdir for files (iterations and scenarios). make a folder for each
# scenario
# rawoutputdir <- paste0(ENTER YOUR DIRECTORY FOR OUTPUT HERE)

dirsin <- c("mortup_9loc",
  "mortup_9loc_singlestart",
  "mortup_8loc_random_20",
  "mortup_8loc_trend_20")

print(Sys.time())
begtime <- Sys.time()

for (i in 1:4) {

# nicedir for cases
case_folder <- paste0(d, "/econ/", dirsin[i], "/eg-cases")

if (file.exists(paste0(rawoutputdir, "/", dirsin[i]))) {
  setwd(paste0(rawoutputdir, "/", dirsin[i]))
} else {
  dir.create(file.path(paste0(rawoutputdir, "/", dirsin[i])),
    showWarnings = FALSE)
  setwd(paste0(rawoutputdir, "/", dirsin[i]))
}

# This is directory where raw output gets written

# The current directory when you register parallel cores is the directory output
# gets written
# registerDoParallel(cores = CHECK THE NUMBER OF WORKERS ON YOUR MACHINE)

run_ss3sim(iterations = 1:200, scenarios =
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

}

endtime <- Sys.time()
print(Sys.time())

endtime - begtime
