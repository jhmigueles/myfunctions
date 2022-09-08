#' covariance
#' @description Calculates the covariance across a set of covariates for Analysis of Compositional Data
#'
#' @param x data frame or matrix with the variables to calculate the covariance
#'
#' @return Matrix with covariance values
#' @export
#'
#' @importFrom gtools permutations
#'
#' @examples
covariance = function(x){
  x = x[complete.cases(x),]
  D = ncol(x)
  perms = permutations(n = D, r = 2, v = 1:D, repeats.allowed = T)

  # covariance matrix
  cov = matrix(data = NA, ncol = D, nrow = D)
  colnames(cov) = rownames(cov) = colnames(x)
  for (i in 1:nrow(perms)) {
    p1 = perms[i,1]; p2 = perms[i,2]
    cov[p1,p2] = var(log(x[,p1] / x[,p2]))
  }
  return(cov)
}
