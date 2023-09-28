#' LR_ILR
#' @description Calculates the isometric log ratio coordinates of a composition for Analysis of Compositional Data.
#'
#' @param compo Matrix or data frame with a composition
#'
#' @return Matrix with isometric log ratios calculated
#' @export
#' @importFrom compositions ilr
LR_ILR = function(compo){
  D = ncol(compo) #n part composition
  z = rep(list(NA), times = D)
  d = 1 #rotate through the basis
  
  # define balance
  balance = matrix(data = NA, 
                   nrow = D - 1, 
                   ncol = D)
  
  for (i in 1:nrow(balance)) {
    if (i == 1) {
      balance[i, ] = c(1, rep(-1, ncol(balance) - 1))
    } else {
      zeroes = sum(balance[i - 1, ] == 0) + 1
      balance[i, 1:zeroes] = 0
      balance[i, zeroes + 1] = 1
      balance[i, is.na(balance[i, ])] = -1
    }
  }
  psi = compositions::gsi.buildilrBase(t(balance))
  while (d <= D) {
    rest = which(1:ncol(compo) != d)
    data = compo[,c(d, rest)]
    tmp = 1
    z[[d]] = as.data.frame(compositions::ilr(data, V = psi))
    # colnames
    while (tmp < D) {
      cols = which(1:D > tmp)
      colnames(z[[d]])[tmp] = paste0(colnames(data)[tmp], "_", paste(colnames(data)[cols], collapse = "."))
      tmp = tmp + 1
    }
    d = d + 1 #next rotation
  }
  return(z)
}
