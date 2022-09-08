#' outliers
#' @description Identify outliers in a given numeric vector
#'
#' @param x Numeric vector
#' @param k Constant to consider outliers (default = 1.5)
#'
#' @return Number of outliers in variable
#' @export
#'
outliers = function(x, k = 1.5){
  iqr = IQR(x, na.rm = T)
  q1 = quantile(x, probs = 0.25, na.rm = T)
  q3 = quantile(x, probs = 0.75, na.rm = T)
  high = q3 + k*iqr
  low = q1 - k*iqr
  outliers_l = which(x < low)
  outliers_h = which(x > high)

  return(length(outliers_l) + length(outliers_h))
}
