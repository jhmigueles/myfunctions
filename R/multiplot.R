#' multiplot
#' @description Place multiple \code{ggplot} plots on one page.
#'
#' @param ... ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' @param plotlist ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' @param cols numeric scalar giving the number of columns in the layout
#' @param layout Matrix with layout for plots. If present, cols is ignored.
#'
#' @return
#' @export
#' @importFrom grid grid.newpage pushViewport grid.layout
#' @examples
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
