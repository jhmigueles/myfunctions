#' gm
#' @description Calculates the geometric mean of a numeric vector
#'
#' @param x Numeric vector
#' @param na.rm Logical. Whether to ignore NA values.
#'
#' @return Geometric mean of the vector
#' @export
#'
gm = function(x, na.rm = T) {
  exp(mean(log(x), na.rm = na.rm))
}
