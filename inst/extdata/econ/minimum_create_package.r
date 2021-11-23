# install.packages("devtools")
# install.packages("roxygen2")
# install.packages("rmarkdown")

library(devtools)

setwd("U:\\R_packages")

setwd("./ss3simfork")
library(rmarkdown)
library(roxygen2)
document()
rmarkdown::render('README.Rmd')
setwd("..")

install("ss3simfork", build_vignettes = FALSE)
