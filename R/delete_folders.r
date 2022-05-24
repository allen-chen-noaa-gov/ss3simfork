setwd(rawoutputdir) #This is directory where RAW output get written

foldernames <- c(
    "D1-E1-F2-M0-cod", 
    "D2-E1-F2-M0-cod", 
    "D3-E1-F2-M0-cod",
    "D4-E1-F2-M0-cod")
    
for (i in 1:100) {

for (j in 1:length(foldernames)) {

unlink(paste0(rawoutputdir, "\\", foldernames[j], "\\", i, "\\em"), recursive=TRUE)
unlink(paste0(rawoutputdir, "\\", foldernames[j], "\\", i, "\\om"), recursive=TRUE)
file.remove(paste0(rawoutputdir, "\\", foldernames[j], "\\", i, "\\log.txt"))

}

}
