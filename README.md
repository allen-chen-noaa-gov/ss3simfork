<!-- README.md is generated from README.Rmd. Please edit that file -->

ss3simfork
==========

This is a fork of the ss3sim package that incorporates a spatial
economic model, testing the use of corrected fishery-dependent abundance
indices.

Please refer to the [ss3sim Github](https://github.com/ss3sim/ss3sim)
for documentation and the original repository. This readme outlines
changes to that original repository.

Please also scroll to the end of this readme for instructions on
replication of the paper “Economists Counting Fish”.

The files added or modified from the original fork should be:

-   /README.Rmd
-   /inst/extdata/econ
    -   /inst/extdata/econ/paper\_script.r
    -   /inst/extdata/econ/make\_fig5.r
    -   …
-   /R/moredatabias.r
-   /R/moredatadahl.r
-   /R/moredatasample.r
-   /R/no\_data.r
-   /R/readabund.r
-   /R/case\_econ.r
-   /R/sample\_econ.r
-   /R/spatial\_fishery.r
-   /R/ss3sim\_base.r (line 478)

Where the economic model interacts with ss3sim
----------------------------------------------

The primary interaction the economic model has with ss3sim is to
potentially add an additional abundance index to the estimating model.
During a standard run, ss3sim produces a “true” abundance index by year
in its operating model. The economic model then takes that true
abundance and distributes it across a spatial grid as defined by the
analyst. Then in each year fishers sample across the spatial grid
depending on some specified sampling pattern, and those samples are used
to construct an estimate of abundance in that year. This index can then
be used in the estimating model.

Running the model
-----------------

Running the model is largely the same as before. A reproducible script
is included below. The ggplot2, RColorBrewer, plyr, doParallel, and
foreach packages can be found on CRAN. Please get ss3simfork and
[barebones.FishSET](https://github.com/allen-chen-noaa-gov/barebones.FishSET)
from their respective Github repositories.

``` r
library(devtools)
# this assumes you've set the directory to wherever you've downloaded the
# packages
install("barebones.FishSET")
install("ss3simfork")
```

The first three packages run the model, and the latter three are just
for making the plots.

``` r
# for running the model
library(ss3simfork)
library(doParallel)
library(foreach)

# for making the plots
library(plyr)
library(ggplot2)
library(RColorBrewer)
```

We use the “cod-like” species that is in the ss3sim package. We also set
two directories, one where we put the cases, and the other where the raw
output is written. You could choose the same directory if you wanted.

``` r
d <- system.file("extdata", package = "ss3simfork")
om <- paste0(d, "/models/cod-om")
em <- paste0(d, "/models/cod-em")

# rawoutputdir for files (iterations and scenarios)
rawoutputdir <- # enter your preferred directory here

# nicedir for cases
nicedir <- # enter your preferred directory here; a folder eg-cases should be
  # in this directory containing the case files.
case_folder <- paste0(nicedir, "/eg-cases")
```

The primary difference is that there is an extra “econ” case file in the
eg-cases folder now. The example cases below can be found in the
`/inst/extdata/econ/eg-cases_readme` folder in the package, to be copied
into the directory you chose above.

An example econ case file is shown below. The econ case file specifies
the economic model and sampling process in `functionname`. There are
currently four options: 1) `nodata` adds no additional data, 2)
`moredatabias` allows fishers to sample based on their preferences, and
uses those samples to construct abundance indices without correction, 3)
`moredatasample` has fishers randomly sample a location, not based on
their preferences, and uses those samples to construct abundance
indices, and 4) `moredatadahl` has fishers sample based on their
preferences, but uses those samples in a process correcting for
selection bias, creating abundance indices.

A matrix describing the distances to each respective location in the
spatial grid is specified in `locnum`. Then, `obsnum` describes the
number of observations simulated at each location, `betavar` the
parameters dictating average catches at each location, and `uparams` the
utility function parameters describing the fisher’s preferences. If you
specify a file name, abundance indices will be written to the specified
folder and name, ending in the ss3sim iteration number. Otherwise it
will be written to each ss3sim iteration folder.

``` r
functionname; nodata
locnum; rbind(c(0.0, 0.5, 0.5, 0.707), c(0.5, 0.0, 0.707, 0.5), c(0.5, 0.707, 0, 0.5), c(0.707, 0.5, 0.5, 0))*3
obsnum; 500
betavar; c(1.50, 1.25, 1.00, 0.75)
uparams; list(alpha = 3, betac = -1)
abundse; NULL
abundscale; NULL
filename; NULL
```

This case file then gets used in the run\_ss3sim call.

You should check how many workers your machine can use. With 8 workers
it took about 2 days.

``` r
# This is directory where RAW output get written
setwd(rawoutputdir)

# The current directory when you register parallel cores is the directory output
# gets written
# registerDoParallel(cores = CHECK THE NUMBER OF WORKERS ON YOUR MACHINE)

print(Sys.time())
begtime <- Sys.time()

run_ss3sim(iterations = 1:100, scenarios =
  c(
  "D1-E1-F2-M0-cod",
  "D2-E1-F2-M0-cod",
  "D3-E1-F2-M0-cod",
  "D4-E1-F2-M0-cod"
  ),
  case_files = list(F = "F", D = c("index", "lcomp", "agecomp", "econ"),
    E = "E", M = "M"),
  case_folder = case_folder, om_dir = om,
  em_dir = em, bias_adjust = FALSE, parallel = TRUE, parallel_iterations = TRUE)

endtime <- Sys.time()
print(Sys.time())

endtime - begtime
```

The operating model runs the same as before. Now, during the estimating
model, in `ss3sim_base`, there is an additional `sample_econ` call that
runs the function specified in the econ case file. As part of this
process the `spatial_fishery` function takes the abundance in each year
and distributes it according to the `betavar` parameters across the
designated spatial grid. Then, depending on the sampling process
specified, an index of abundance is estimated.

Some results
------------

Using the function from the ss3sim package, we can grab aggregated
results `scalar_dat` and `ts_dat`.

``` r
get_results_all(directory = rawoutputdir,
  overwrite_files = TRUE,
  # this writes to rawoutputdir even if you set a nicedir before call
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

ts_dat <- merge(ts_dat, scalar_dat[, c("scenario", "replicate",
  "max_grad")])

subsetlist <- c("D1", "D2", "D3", "D4")

ts_dat_sto <- subset(ts_dat, D %in% subsetlist)

ts_dat_sto$D <- as.factor(ts_dat_sto$D)

ts_dat_sto$D <- revalue(ts_dat_sto$D, c(
  "D1" = "Quadrennial fishery-independent survey",
  "D2" = "+ Annual randomly-sampled fishery data",
  "D3" = "+ Annual uncorrected fishery data",
  "D4" = "+ Annual corrected fishery data"))

# Small values of the maximum gradient (approximately 0.001 or less) indicate
# that model convergence is likely.
# Larger values (greater than 1) indicate that model convergence is unlikely.
ts_dat_sto <- ts_dat_sto[ts_dat_sto$max_grad < 1, ]

test.50 <- ddply(ts_dat_sto, .(D, year), summarise,
  five0 = quantile(SpawnBio, probs = c(0.5), na.rm = TRUE))

test.05 <- ddply(ts_dat_sto, .(D, year), summarise,
  zero5 = quantile(SpawnBio, probs = c(0.025), na.rm = TRUE))
test.95 <- ddply(ts_dat_sto, .(D, year), summarise,
  nine5 = quantile(SpawnBio, probs = c(0.975), na.rm = TRUE))
test.25 <- ddply(ts_dat_sto, .(D, year), summarise,
  two5 = quantile(SpawnBio, probs = c(0.25), na.rm = TRUE))
test.75 <- ddply(ts_dat_sto, .(D, year), summarise,
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
  theme(text = element_text(size = 30),
    axis.title.y = element_text(vjust = 1.5)) +
  geom_line(data = test.50, aes(y = five0, x = year), colour = "black") +
  ylim(-0.55, 0.55) + ylab("Relative Error in Spawning Biomass") +
  xlab("Year") +
  facet_wrap(~D, ncol = 2)
print(p)
```

![Spawning biomass output
example](https://github.com/allen-chen-noaa-gov/ss3simfork/blob/master/inst/extdata/econ/readme_SB_081821_mortup4.png?raw=true)

We can also look at the relative abundance indices across different
sampling patterns.

``` r
abundout <- list()
for (it in 1:100) {

# The default naming convention for each abundance index .csv is below, if
# the file name is left NULL, but can be specified if desired.

casename <- paste0(rawoutputdir, "/D3-E1-F2-M0-cod/", it,
  "/bias-abund-D3-E1-F2-M0-cod-")
abundtrue <- readabund(filename = paste(casename, it, ".csv", sep = ""),
  type = "True", istruecpue = TRUE)

casename <- paste0(rawoutputdir, "/D3-E1-F2-M0-cod/", it,
  "/bias-abund-D3-E1-F2-M0-cod-")
abundbias <- readabund(filename = paste(casename, it, ".csv", sep = ""),
  type = "Uncorrected")

casename <- paste0(rawoutputdir, "/D4-E1-F2-M0-cod/", it,
  "/cor-abund-D4-E1-F2-M0-cod-")
abundcor <- readabund(filename = paste(casename, it, ".csv", sep = ""),
  type = "Corrected")

casename <- paste0(rawoutputdir, "/D2-E1-F2-M0-cod/", it,
"/samp-abund-D2-E1-F2-M0-cod-")
abundsamp <- readabund(filename = paste(casename, it, ".csv", sep = ""),
    type = "Randomly sampled")

tempabundout <- rbind(abundtrue, abundbias, abundcor, abundsamp)

tempabundout$iter <- it

abundout[[it]] <- tempabundout

}

totabundout <- do.call(rbind, abundout)

totabundout <- totabundout[totabundout$CPUE > 0, ]

linegraph <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), median)
names(linegraph)[names(linegraph) == "x"] <- "RelDiff"

linegraphse <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), sd)
names(linegraphse)[names(linegraphse) == "x"] <- "sd"

linegraphlen <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), length)
names(linegraphlen)[names(linegraphlen) == "x"] <- "len"

linegraph <- merge(linegraph, linegraphse, by = c("Year", "Index"))
linegraph <- merge(linegraph, linegraphlen, by = c("Year", "Index"))

linegraph$se <- linegraph$sd/sqrt(linegraph$len)

linegraph$Index <- revalue(linegraph$Index, c("Corrected" = "Corrected     ",
  "Randomly sampled" = "Randomly sampled     ",
  "True" = "True     ",
  "Uncorrected" = "Uncorrected     "))

pd <- position_dodge(0.05)

lineout <- ggplot(data = linegraph,
  aes(x = Year, y = RelDiff)) +
  geom_line(aes(colour = Index, group = Index), size = 2, position = pd) +
  geom_point(aes(shape = Index), size = 4, fill = "white", stroke = 1.5,
    position = pd) +
  xlab("Year") + ylab("Relative Differences in Abundance") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "orange",
    size = 2) +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme(legend.direction = "horizontal", legend.position = "top",
    legend.key = element_rect(size = 4),
    legend.text = element_text(size = 22),
    legend.key.size = unit(2.5, "lines"), legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 30),
    axis.title.y = element_text(vjust = 1.5)) +
  scale_shape_manual(values = c(21, 22, 25, 24)) +
  scale_colour_manual(values = brewer.pal(4, "Set1"))
print(lineout)
```

![Abundance indices output
example](https://github.com/allen-chen-noaa-gov/ss3simfork/blob/master/inst/extdata/econ/readme_abundline_081821_mortup4.png?raw=true)

Economists Counting Fish
------------------------

We also include the code to reproduce the results from the paper
“Economists Counting Fish”. A script can be found in the
`/inst/extdata/econ/` folder in the package, titled `paper_script.r`.
Note this script requires you to set a directory to write the raw output
into, as the raw output is expected to be a few GB.

Then, the `make_fig5.r` script will reproduce Figure 5 from the paper,
which outlines the various relative abundance indices and errors. With
8 workers (Intel i7-8700), it took a little over 16 days.

DISCLAIMER
------------------------
“This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.”
