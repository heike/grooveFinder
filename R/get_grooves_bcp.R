#' Use Bayesian Changepoint Detection to identify groove locations
#'
#' @inheritParams get_grooves_quadratic
#' @importFrom bulletcp get_grooves_bcp
#' @export
get_grooves_bcp <- function(x, value, adjust = 10, return_plot = FALSE) {
  bcp_out <- bulletcp::get_grooves_bcp(x = x, value = value, adjust = adjust, ...)
  if (return_plot == T) {
    grooves <- list(
      groove = bcp_out$groove,
      plot = grooves_plot(land = land, grooves = bcp_out$groove)
    )
  } else {
    grooves <- bcp_out
  }

grooves
}
