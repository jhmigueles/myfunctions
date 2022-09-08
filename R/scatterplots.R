#' scatterplots
#'
#' @param x Numeric vector (variable x axis)
#' @param y Numeric vector (variable y axis)
#' @param ylim Y axis limits
#' @param pch pch as in plot
#' @param main Title of the plot
#' @param xlab X axis label
#' @param ylab Y axis label
#' @param data Data frame
#'
#' @return
#' @export
#'
#' @examples
scatterplots = function(x, y, ylim=c(), pch=c(), main=c(), xlab=c(), ylab=c(), data = c()){

  # Pearson
  correlation = cor.test(x, y)

  # Spearman
  rank.x=data.frame(rank=1:(length(x)),id=data$ID[order(x)])
  rank.y=data.frame(rank=1:(length(y)),id=data$ID[order(y)])
  rankdat.xy=data.frame(rank.x=rank.x[order(rank.x$id),'rank'],
                        rank.y=rank.y[order(rank.y$id),'rank'])
  spearman = cor.test(x=rankdat.xy$rank.x,y=rankdat.xy$rank.y,alternative='greater',method='spearman')


  plot(x, y,
       pch = pch, ylim = ylim, main = main,
       xlab = xlab, ylab = ylab)
  abline(0, 1, col = "grey")
  legend("topleft", paste0("Pearson r = ", round(correlation$estimate, 3),
                           "\nSpearman rank r = ", round(spearman$estimate, 3)), bty = "n")

}
