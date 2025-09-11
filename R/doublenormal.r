
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
} 