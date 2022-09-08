#' bland.altman
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param ylim Y axis limits
#' @param pch pch as in plot
#' @param main Title of the plot
#' @param xlab X axis label
#' @param ylab Y axis label
#'
#' @return
#' @export
#'
#' @examples
bland.altman = function(x, y, ylim=c(), pch=c(), main=c(), xlab=c(), ylab=c()){
  means = (x + y) / 2
  diffs = x - y

  mean_bias = mean(x - y, na.rm=T)
  loa = 1.96 * sd(x - y, na.rm=T)
  loa = c(mean_bias - loa, mean_bias + loa)

  if(length(ylim)==0){
    ylim = c(min(diffs, loa, na.rm = T), max(diffs, loa, na.rm = T))
  }

  plot(means, diffs,
       pch = pch, ylim = ylim, main = main,
       xlab = xlab, ylab = ylab)
  abline(h = c(mean_bias, loa), lty = c(1,2,2))
}
