#' Quadratic fit to find groove locations
#'
#' Use a robust fit of a quadratic curve to find groove locations
#' @param x numeric vector of locations (in microns)
#' @param value numeric values of surface measurements in microns
#' @param adjust positive number to adjust the grooves
#' @param return_plot return plot of grooves?
#' @return list of groove vector and plot of crosscut with shoulder locations
#' @importFrom MASS rlm
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' \dontrun{
#' # Set the data up to be read in, cleaned, etc.
#' library(bulletxtrctr)
#' library(x3ptools)
#'
#' example_data <- bullet_pipeline(
#'   location = list(Bullet1 = c(hamby252demo$bullet1[3])),
#'   stop_at_step = "crosscut",
#'   x3p_clean = function(x) x %>%
#'       x3p_scale_unit(scale_by = 10^6) %>%
#'       rotate_x3p(angle = -90) %>%
#'       y_flip_x3p()
#' )
#'
#' get_grooves_quadratic(example_data$ccdata[[1]]$x,
#'   example_data$ccdata[[1]]$value,
#'   adjust = 30, return_plot = T
#' )
#' cc_locate_grooves(example_data$ccdata[[1]],
#'   method = "quadratic",
#'   adjust = 30, return_plot = T
#' )
#' }
get_grooves_quadratic <- function(x, value, adjust, return_plot = F) {
  assert_that(
    is.numeric(x), is.numeric(value), is.numeric(adjust),
    is.logical(return_plot)
  )

  land <- data.frame(x = x, value = value)

  lm0 <- rlm(value ~ poly(x, 2), data = land, maxit = 100)
  land$pred <- predict(lm0, newdata = land)

  land$absresid <- with(land, abs(value - pred))
  absresid90 <- NULL
  land$absresid90 <- with(
    land, absresid > 4 * median(land$absresid, na.rm = TRUE)
  )

  groove <- range(filter(land, !absresid90)$x) + c(adjust, -adjust)

  if (return_plot) {
    return(list(
      groove = groove,
      plot = grooves_plot(land = land, grooves = groove)
    ))
  } else {
    return(list(groove = groove))
  }
}
