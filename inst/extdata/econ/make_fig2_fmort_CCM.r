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

## ## plotdata <- plotdata[plotdata$years > 40, ]
## lineout <- ggplot(data = plotdata,
##     aes(x = years, y = fvals)) +
##     geom_point(shape = 21, colour = "black", fill = "white", size = 2.5,
##         stroke = 2.5) +
##     xlab("Year") + ylab("Fishing mortality") +
##     guides(colour = guide_legend(override.aes = list(size = 6))) +
##     theme(legend.direction = "horizontal", legend.position = "top",
##         legend.key = element_rect(size = 4),
##         legend.text = element_text(size = 22),
##         legend.key.size = unit(2.5, "lines"), legend.title = element_blank(),
##         axis.text.x = element_text(angle = 45, hjust = 1),
##         text = element_text(size = 30),
##         axis.title.y = element_text(vjust = 1.5))
## ggsave(lineout, file = paste0(getwd(), "/inst/extdata/econ/fmortup",
##   ".png"), width = 16, height = 9)

abundex <- read.csv(paste0(getwd(), "/inst/extdata/econ/abundex.csv"))
abundex <- subset(abundex, abundex$index == 2)
abundex <- abundex[abundex$year > 0, ]
abundex$obs <- abundex$obs/(10^9)

## lineout <- ggplot(data = abundex,
##     aes(x = year, y = obs)) +
##     geom_point(shape = 21, colour = "black", fill = "white", size = 2.5,
##         stroke = 2.5) +
##     xlab("Year") + ylab("Abundance (10^9)") +
##     guides(colour = guide_legend(override.aes = list(size = 6))) +
##     ylim(0, max(abundex$obs)) +
##     theme(legend.direction = "horizontal", legend.position = "top",
##         legend.key = element_rect(size = 4),
##         legend.text = element_text(size = 22),
##         legend.key.size = unit(2.5, "lines"), legend.title = element_blank(),
##         axis.text.x = element_text(angle = 45, hjust = 1),
##         text = element_text(size = 30),
##         axis.title.y = element_text(vjust = 1.5))
## ggsave(lineout, file = paste0(getwd(),
##   "/inst/extdata/econ/fmortup_abundex", ".png"), width = 16, height = 9)




#' Calculate the double-normal selectivity curve
#'
#' @param x Vector of lengths or ages
#' @param sp Vector of 6 representing the parameters of the curve
#'
#' @return Vector of selectivity values
#'
#' @details Modified from r4ss function
doublenormal <- function(x, sp) {
  ## Stole this from r4ss and modified to return selex. Cole 4/25/2017.
  sel <- rep(NA, length(x))
  startbin <- 1
  peak <- sp[1]
  upselex <- exp(sp[3])
  downselex <- exp(sp[4])
  final <- sp[6]
  if(sp[5] < -1000) {
    j1 <-  -1001 - round(sp[5])
    sel[1:j1] <- 1.0e-06
  }
  if(sp[5] >= -1000) {
    j1 <- startbin - 1
    if(sp[5] > -999) {
      point1 <- 1.0/(1.0+exp(-sp[5]))
      t1min <- exp(-(x[startbin]-peak)^2 / upselex)
    }
  }
  if(sp[6] < -1000) j2 <- -1000- round(sp[6])
  if(sp[6] >= -1000) j2 <- length(x)
  peak2 <- peak + 2 + (0.99*x[j2]- peak - 2)/(1.+exp(-sp[2]))
  if(sp[6] > -999) {
    point2 <- 1.0/(1.0 + exp(-final))
    t2min <- exp(-(x[j2]-peak2)^2 / downselex)
  }
  t1 <- x - peak
  t2 <- x - peak2
  join1 <- 1.0/(1.0 + exp(-(20./(1.0 + abs(t1)))*t1))
  join2 <- 1.0/(1.0 + exp(-(20./(1.0 + abs(t2)))*t2))
  if(sp[5] > -999) asc <- point1 + (1.0-point1) * (exp(-t1^2 / upselex)-t1min)/(1.0-t1min)
  if(sp[5] <= -999) asc <- exp(-t1^2 / upselex)
  if(sp[6] > -999) dsc <- 1.0 + (point2 - 1.0) * (exp(-t2^2 / downselex)-1.0) / (t2min-1.0)
  if(sp[6] <= -999) dsc <- exp(-(t2)^2/downselex)
  sel[(j1+1):j2] <- asc*(1.0-join1)+join1*(1.0-join2+dsc*join2)
  if(startbin > 1 && sp[5] >= -1000) {
    sel[1:startbin] <- (x[1:startbin] / x[startbin])^2 * sel[startbin]
  }
  if(j2 < length(x)) sel[(j2+1):length(x)] <- sel[j2]
  return(sel)
} # end sel.line function

x <- 0:180

pf <- c(50.8, -3,5.1,15,-999,-999)
ps <- c(41.8,-4,5.2,14,-99,99)
self <- doublenormal(x, pf)
sels <- doublenormal(x, ps)

## journal guidelines are
width1 <- 170/ 25.4 # mm to in
width2 <- 85/ 25.4 # mm to in

print.letter <- function(label="(a)",xy=c(0.1,0.925),...) {
    tmp <- par("usr")
    text.x <- tmp[1]+xy[1]*diff(tmp[1:2])   #x position, diff=difference
    text.y <- tmp[3]+xy[2]*diff(tmp[3:4])   #y position
    text(x=text.x, y=text.y, labels=label, ...)
}

xy <- c(.04,.95)
png(paste0(getwd(), "/inst/extdata/econ/fig2.png"), width = width2,
    height = 6.5, units='in', res=500)
## pdf(paste0(getwd(), "/inst/extdata/econ/fig2.pdf"), width = width2,
##     height = 6.5)
par(mfrow=c(3,1), mar=c(3,3,.5,.5), mgp=c(1.5,.4,0), tck=-.02)
plot(plotdata$years, plotdata$fvals, type='l', xlab='Year',
     ylab='Fishing effort (F)', xlim=c(40,100))
print.letter('(a)', xy)
plot(abundex$year, abundex$obs, xlim=c(40,100), type='l',
     xlab='Year', ylab='Biomass (M t)', ylim=c(0,5))
print.letter('(b)', xy)
plot(x, self, type='l', xlab='Length (cm)', ylab='Selectivity', xlim=c(0,100))
lines(x,sels, col=1, lty=2)
legend('bottomright', legend=c('Fishery', 'Survey'), lty=c(1,2), bty='n')
print.letter('(c)', xy)
dev.off()
