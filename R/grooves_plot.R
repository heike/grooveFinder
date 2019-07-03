#' Internal function to plot crosscut + grooves
#'
#' @param land data.frame with columns x and value
#' @param grooves numeric vector of length 2 identifying both grooves.
#'        If only one groove is identified, the other should be NA
#' @importFrom assertthat assert_that
#' @import ggplot2
#' @return a ggplot2 object
grooves_plot <- function(land, grooves) {
  assert_that(has_name(land, "x"))
  assert_that(has_name(land, "value"))
  assert_that(is.numeric(grooves), sum(!is.na(grooves)) >= 1)
  x <- value <- NULL

  ggplot(aes(x = x, y = value), data = land) + geom_line(size = .5) +
    theme_bw() +
    geom_vline(xintercept = grooves[1], colour = "blue") +
    geom_vline(xintercept = grooves[2], colour = "blue")
}

