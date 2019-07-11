#' Using rollapply to find grooves in a crosscut
#'
#' @param x numeric vector of locations (in microns)
#' @param value numeric values of surface measurements in microns
#' @param smoothfactor The smoothing window to use
#' @param adjust positive number to adjust the grooves
#' @param groove_cutoff The index at which a groove cannot exist past
#' @param mean_left If provided, the location of the average left groove
#' @param mean_right If provided, the location of the average right groove
#' @param mean_window The window around the means to use
#' @param second_smooth Whether or not to smooth a second time
#' @param which_fun Which function to use in the rollapply statement
#' @param return_plot return plot of grooves?
#' @export
#' @importFrom assertthat assert_that
#' @importFrom zoo rollapply na.fill
#' @importFrom utils head tail
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
#' get_grooves_rollapply(example_data$ccdata[[1]]$x,
#'   example_data$ccdata[[1]]$value,
#'   adjust = 30, return_plot = T
#' )
#' cc_locate_grooves(example_data$ccdata[[1]],
#'   method = "rollapply",
#'   adjust = 30, return_plot = T
#' )
#' }
get_grooves_rollapply <- function(x, value, smoothfactor = 15, adjust = 10,
                                  groove_cutoff = 400, mean_left = NULL,
                                  mean_right = NULL, mean_window = 100,
                                  second_smooth = T, which_fun = mean,
                                  return_plot = F) {
  . <- NULL

  assert_that(
    is.numeric(x), is.numeric(value), is.numeric(adjust),
    is.numeric(smoothfactor), is.numeric(groove_cutoff),
    is.logical(second_smooth), is.logical(return_plot),
    "function" %in% class(which_fun)
  )

  land <- data.frame(x = x, value = value)
  original_land <- land

  if (!is.null(mean_left) && !is.null(mean_right)) {
    mean.left.ind <- which.min(abs(land$x - mean_left))
    mean.right.ind <- which.min(abs(land$x - mean_right))

    window.left.left <- max(1, mean.left.ind - mean_window)
    window.left.right <- mean.left.ind + mean_window

    window.right.left <- mean.right.ind - mean_window
    window.right.right <- min(length(land$x), mean.right.ind + mean_window)

    land <- land[c(
      window.left.left:window.left.right,
      window.right.left:window.right.right
    ), ]

    groove_cutoff <- Inf
  }

  value_filled <- zoo::na.fill(land$value, "extend")
  smoothed <- rollapply(value_filled, smoothfactor, function(x) which_fun(x))

  # Add in an if statement, to only do the first smoothing if the second_smooth
  # parameter is equal to FALSE
  if (second_smooth) {
    smoothed_truefalse <- rollapply(smoothed, smoothfactor,
      function(x) which_fun(x),
      partial = TRUE
    )
  } else {
    smoothed_truefalse <- smoothed
  }

  lengthdiff <- length(land$value) - length(smoothed_truefalse)

  peak_ind_smoothed <- rollapply(
    smoothed_truefalse, 3,
    function(x) which.max(x) == 2
  ) %>%
    which() %>%
    head(n = 1)
  peak_ind <- peak_ind_smoothed + floor(lengthdiff / 2)
  if (length(peak_ind) == 0) {
    groove_ind <- peak_ind
  } else {
    groove_ind <- tail(smoothed_truefalse, n = -peak_ind_smoothed) %>%
      rollapply(., 3, function(x) which.min(x) == 2) %>%
      which() %>%
      head(n = 1)
    groove_ind <- groove_ind + peak_ind
  }

  peak_ind2_smoothed_temp <- smoothed_truefalse %>%
    rev() %>%
    rollapply(., 3, function(x) which.max(x) == 2) %>%
    which() %>%
    head(n = 1)

  peak_ind2_temp <- peak_ind2_smoothed_temp + floor(lengthdiff / 2)
  if (length(peak_ind2_temp) == 0) {
    groove_ind2_temp <- peak_ind2_temp
  } else {
    groove_ind2_temp <- rev(smoothed_truefalse) %>%
      tail(., n = -peak_ind2_smoothed_temp) %>%
      rollapply(., 3, function(x) which.min(x) == 2) %>%
      which() %>%
      head(n = 1)
    groove_ind2_temp <- groove_ind2_temp + peak_ind2_temp
  }

  # peak_ind2 <- length(land$value) - peak_ind2_temp + 1
  groove_ind2 <- length(land$value) - groove_ind2_temp + 1

  ## Check that it actually FOUND a groove...
  if (length(groove_ind) == 0 || groove_ind > groove_cutoff) {
    groove_ind <- 1
  }
  if (length(groove_ind2) == 0 ||
    groove_ind2 < length(land$value) - groove_cutoff) {
    groove_ind2 <- length(land$value)
  }

  xvals <- original_land$x
  # yvals <- original_land$value

  # plot_peak_ind <- which(original_land$x == land$x[peak_ind])
  plot_groove_ind <- which(original_land$x == land$x[groove_ind])
  # plot_peak_ind2 <- which(original_land$x == land$x[peak_ind2])
  plot_groove_ind2 <- which(original_land$x == land$x[groove_ind2])

  center <- which.min(abs(xvals - mean(xvals)))

  # I can't figure out how to test this if statement...
  if (plot_groove_ind > center) {
    plot_groove_ind2 <- plot_groove_ind
    plot_groove_ind <- 0
  }

  if (plot_groove_ind2 < center) {
    plot_groove_ind <- plot_groove_ind2
    plot_groove_ind2 <- length(xvals)
  }

  # smoothed_diff <- floor(lengthdiff/2)

  groove <- c(
    original_land$x[plot_groove_ind + adjust],
    original_land$x[plot_groove_ind2 - adjust]
  )

  if (return_plot) {
    return(list(
      groove = groove,
      plot = grooves_plot(land = original_land, grooves = groove)
    ))
  } else {
    return(list(groove = groove))
  }
}
