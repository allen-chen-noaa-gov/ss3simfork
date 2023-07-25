library(ggplot2)

fvals <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0.00432478632478633, 0.00864957264957265, 0.012974358974359, 0.0172991452991453,
0.0216239316239316, 0.025948717948718, 0.0302735042735043, 0.0345982905982906,
0.0389230769230769, 0.0432478632478633, 0.0475726495726496, 0.0518974358974359,
0.0562222222222222, 0.0605470085470085, 0.0648717948717949, 0.0691965811965812,
0.0735213675213675, 0.0778461538461539, 0.0821709401709402, 0.0864957264957265,
0.0908205128205128, 0.0951452991452992, 0.0994700854700855, 0.103794871794872,
0.108119658119658, 0.112444444444444, 0.116769230769231, 0.121094017094017,
0.125418803418803, 0.12974358974359, 0.134068376068376, 0.138393162393162,
0.142717948717949, 0.147042735042735, 0.151367521367521, 0.155692307692308,
0.160017094017094, 0.16434188034188, 0.168666666666667)
years <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
100)

plotdata <- data.frame(cbind(years, fvals))
plotdata <- plotdata[plotdata$years > 40, ]
lineout <- ggplot(data = plotdata,
    aes(x = years, y = fvals)) + 
    geom_point(shape = 21, colour = "black", fill = "white", size = 2.5,
        stroke = 2.5) +
    xlab("Year") + ylab("Fishing mortality") +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    theme(legend.direction = "horizontal", legend.position = "top",
        legend.key = element_rect(size = 4),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2.5, "lines"), legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 30),
        axis.title.y = element_text(vjust = 1.5))

ggsave(lineout, file = paste0(getwd(), "/inst/extdata/econ/fmortup",
  ".png"), width = 16, height = 9)

#this is from iteration one of base scenario
abundex <- read.csv(paste0(getwd(), "/inst/extdata/econ/abundex.csv"))

abundex <- subset(abundex, abundex$index == 2)
abundex <- abundex[abundex$year > 0, ]
abundex$obs <- abundex$obs/(10^9)

lineout <- ggplot(data = abundex,
    aes(x = year, y = obs)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 2.5,
        stroke = 2.5) +
    xlab("Year") + ylab("Abundance (10^9)") +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    ylim(0, max(abundex$obs)) +
    theme(legend.direction = "horizontal", legend.position = "top",
        legend.key = element_rect(size = 4),
        legend.text = element_text(size = 22),
        legend.key.size = unit(2.5, "lines"), legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 30),
        axis.title.y = element_text(vjust = 1.5))

ggsave(lineout, file = paste0(getwd(),
  "/inst/extdata/econ/fmortup_abundex", ".png"), width = 16, height = 9)
