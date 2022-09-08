#' LR_ILR
#' @description Calculates the isometric log ratio coordinates of a composition for Analysis of Compositional Data.
#'
#' @param compo Matrix or data frame with a composition
#'
#' @return
#' @export
#'
#' @examples
LR_ILR = function(compo){
  D = ncol(compo) #n part composition
  z = rep(list(NA), times = D)
  d = 1 #rotate through the basis
  while (d <= D) {
    rest = which(1:ncol(compo) != d)
    data = compo[,c(d, rest)]
    tmp = 1
    while (tmp < D) {
      if (length(z[[d]]) == 1) z[[d]] = data[, 1:(D - 1)]
      cols = which(1:D > tmp)
      for (i in 1:nrow(z[[d]])) {
        z[[d]][i, tmp] = sqrt((D - tmp)/(D - tmp + 1))*log(data[i,tmp]/prod(data[i,cols])^(1/(D - tmp)))
      }
      colnames(z[[d]])[tmp] = paste0(colnames(data)[tmp], "_", paste(colnames(data)[cols], collapse = "."))
      tmp = tmp + 1
    }
    d = d + 1 #next rotation
  }
  return(z)
}
