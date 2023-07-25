library(png)
library(sjPlot)

png("tiled_fig1.png", width = 400, height = 225, units='mm', res = 300)

plot(NA, xlim = c(0, 2), ylim = c(0, 2), type = "n", xaxt = "n", yaxt = "n",
  xlab = "", ylab = "")

rasterImage(readPNG(source = "fmortup.png"), 0, 1, 1, 2)
rasterImage(readPNG(source = "fmortup_abundex.png"), 1, 1, 2, 2)
rasterImage(readPNG(source = "fishery_selex.png"), 0, 0, 1, 1)
rasterImage(readPNG(source = "survey_selex.png"), 1, 0, 2, 1)

dev.off()
