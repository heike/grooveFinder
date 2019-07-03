#' Use the center of a crosscut
#'
#' @param x numeric vector of locations in microns
#' @param value numeric vector of surface measurements in microns
#' @param middle middle percent to use for the identification
#' @param return_plot return plot?
#' @return list of groove vector and plot of crosscut, if return_plot is true
#' @importFrom assertthat assert_that
#' @importFrom dplyr filter mutate group_by summarize count
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
#'       x3p_scale_unit(scale_by=10^6) %>%
#'       rotate_x3p(angle = -90) %>%
#'       y_flip_x3p()
#' )
#'
#' get_grooves_middle(example_data$ccdata[[1]]$x,
#'   example_data$ccdata[[1]]$value,
#'   return_plot = T
#' )
#' cc_locate_grooves(example_data$ccdata[[1]],
#'   method = "middle",
#'   return_plot = T
#' )
#' }
get_grooves_middle <- function(x, value, middle = 75, return_plot = F) {
  assert_that(
    is.numeric(x), is.numeric(value), is.numeric(middle),
    is.logical(return_plot)
  )
  assert_that(middle >= 0, middle <= 100)

  land <- data.frame(x = x, value = value)
  # summarize values for each x:
  ns <- land %>% count(x)
  if (max(ns$n) > 1) message(sprintf("summarizing across %d profiles ...", max(ns$n)))

  land <- land %>% group_by(x) %>% summarize(value = mean(value, na.rm = TRUE))
  groove <- quantile(land$x,
                     probs = c((100 - middle) / 200, (100 + middle) / 200)
  )

  if (return_plot) {
    return(list(
      groove = groove,
      plot = grooves_plot(land = land, grooves = groove)
    ))
  } else {
    return(list(groove = groove))
  }
}
