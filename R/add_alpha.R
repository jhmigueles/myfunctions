#' add_alpha
#' @description A function to add alpha levels to RGB colors
#'
#' @param col Color
#' @param alpha Alpha level (between 0 and 1)
#'
#' @return
#' @export
#'
#' @importFrom grDevices rgb
#'
#' @examples
add_alpha <- function(col, alpha = 1){
  if (missing(col)) stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha = alpha))
}