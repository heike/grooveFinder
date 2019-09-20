

#' Re-parameterize the Hough transform output for pixel identification
#'
#' Helper function, that should not be used outside the package.
#' Computes x and y intercepts for a line given in polar coordinate representation.
#' A line in polar coordinate representation is given by angle theta (with respect to the x-axis) and
#' distance of the line to the origin rho.
#' @param rho Numeric vector containing the shortest distance from the line to the origin
#' @param theta Numeric vector containing the angle of the line from the positive x axis
#' @param df Data frame containing output from a Hough transformation
#' @return data frame with variables rho, theta, score (original data frame) expanded by yintercept and xintercept.


rho_to_ab <- function(rho = NULL, theta = NULL, df = NULL) {
  if (is.null(df)) {
    df <- data.frame(rho = rho, theta = theta)
  }
  stopifnot(c("rho", "theta") %in% names(df))

  df <- df %>%
    mutate(
      yintercept = ifelse(theta == 0, NA, rho / sin(theta)),
  #    slope = -cos(theta) / sin(theta), # don't compute slope at this point - it's ambiguous
      xintercept = rho / cos(theta)
    ) # cos(theta) == 0 for theta = ± pi/2 but realistically we don't get there
  df
}



#' Use Hough transformation to identify groove locations.
#'
#' Hough transformations are used to identify the location of left and right groove.
#' Generates Hough lines
#' and selects lines with angle corresponding to (close to) vertical lines.
#'
#' @param land dataframe of surface measurements in microns in the x, y, and x direction. Use `x3p_to_df` to access the data from an x3p scan.
#' @param qu quantile (between 0 and 1) to specify score quantile for which vertical lines are considered. Defaults to qu = 0.999
#' If groove are not strongly expressed, lower this threshold.
#' @param adjust (generally) positive number in micron used to adjust the grooves inward
#' @param return_plot boolean value - should a plot of the crosscut with the grooves be returned? defaults to FALSE
#' @return list object consisting of functions to describe the left and right groove.
#' Parameters for the functions are given in microns and return results in microns.
#'
#' @importFrom x3ptools df_to_x3p
#' @importFrom imager as.cimg width height
#' @importFrom imager imgradient
#' @importFrom imager hough_line
#' @importFrom assertthat assert_that
#' @importFrom assertthat has_name
#' @importFrom stats quantile median sd na.omit
#' @importFrom dplyr filter mutate group_by summarize count
#' @importFrom x3ptools df_to_x3p x3p_to_df x3p_get_scale
#'
#' @examples
#' library(x3ptools)
#' library(ggplot2)
#' data("br411")
#' x3p <- br411
#'
#' # Get grooves fit for left and right groove from x3p
#' grooves <- get_grooves_hough(x3p_to_df(x3p), qu = 0.999, adjust = 0)
#'
#' # Find optimized crosscut location, this may take some time
#' # if (require(bulletxtrctr)) {
#' #  crosscut <- x3p %>% bulletxtrctr::x3p_crosscut_optimize()
#' # } else {
#'   crosscut <- 125
#' # }
#'
#' \dontrun{
#' a <- get_mask_hough(x3p, grooves)
#' a <- a %>% x3p_add_hline(yintercept=crosscut)
#' x3ptools::image_x3p(a)
#' }
#'
#' # Find groove locations for specified crosscut
#' grooves$left.groove.fit(crosscut)
#' grooves$right.groove.fit(crosscut)
#'
#' # Plot profile
#' ccdata <- x3p %>% x3p_to_df() %>% dplyr::filter(y == crosscut)
#'
#' ccdata %>%
#'   ggplot(aes(x = x, y = value))+
#'   geom_line()+
#'   geom_vline(xintercept = grooves$left.groove.fit(crosscut), color = "red") +
#'   geom_vline(xintercept = grooves$right.groove.fit(crosscut), color = "blue")
#'
#' # grooves at a different crosscut:
#' x3p %>% x3p_to_df() %>% dplyr::filter(y == 150) %>%
#'   ggplot(aes(x = x, y = value))+
#'   geom_line()+
#'   geom_vline(xintercept = grooves$left.groove.fit(crosscut), color = "red") +
#'   geom_vline(xintercept = grooves$right.groove.fit(crosscut), color = "blue")
#' @export
get_grooves_hough <- function(land, qu = 0.999, adjust = 10, return_plot = FALSE) {
  assert_that(
    has_name(land, "x"), has_name(land, "y"), has_name(land, "value"),
    is.numeric(land$x), is.numeric(land$y), is.numeric(land$value)
  )

  # Convert to cimage
  land.x3p <- df_to_x3p(land)

  # helper function to convert from pixel to microns
  pix_to_micron <- function(x) {
    (x-1) * x3p_get_scale(land.x3p)
  }


  # visible bindings problem
  theta <- score <- rho <- xintercept <- yintercept <- 0
  slope <- x <- value <- 0

  # HH: this if condition distinguishes between lands based on the aspect ratio.
  # please assume that we have lands that are wider than high.
  # If they are not, spit out a warning rather than different code.
  # It has bitten us badly into rear elements before to try to think for the user.
  sizes <- dim(land.x3p$surface.matrix)
  width <- sizes[1] # faster to use than function call to width(strong)
  height <- sizes[2] # same
  if (width < height) {
    # cimg <- as.cimg(t(land.x3p$surface.matrix))
    warning("This scan seems to not be rotated correctly. Proceed with caution. ")
  }

  # else{
  cimg <- as.cimg(land.x3p$surface.matrix)
  # }

  # Create image gradient
  dx <- imgradient(cimg, "x")
  dy <- imgradient(cimg, "y")

  grad.mag <- sqrt(dx^2 + dy^2)

  strong <- grad.mag > quantile(grad.mag, .99, na.rm = TRUE)
  # create the hough transform
  hough.df <- hough_line(strong, data.frame = TRUE, shift = FALSE) # we want to get values with respect to (0,0) not (1,1)
  assert_that(
    has_name(hough.df, "theta"),
    has_name(hough.df, "rho"),
    has_name(hough.df, "score")
  )

  # Subset based on score and angle rotation
  hough.df <- hough.df %>%
    dplyr::mutate(theta = ifelse(theta <= pi, theta, theta - 2 * pi)) %>%
    dplyr::filter(
      score > quantile(score, qu),
      theta > (-pi / 16), # identify only vertical(ish) lines
      theta < (pi / 16),
      (rho < abs(width) * 1 / 6 | rho > width * 5 / 6) # at either end of the LEA
    )
  hough.df <- unique(hough.df) # get rid of duplicates

  # get x and y intercepts  (in pixel dimensions)
  segments <- rho_to_ab(df = hough.df)
  assert_that(
    has_name(segments, "theta"),
    has_name(segments, "rho"),
    has_name(segments, "score"),
    has_name(segments, "xintercept"),
    has_name(segments, "yintercept") #,
  #  has_name(segments, "slope") # don't have a slope any more
  )

  if (nrow(segments) == 0) stop(sprintf("No results found. Try to reduce the quantile qu, currently at %f.", qu))

  # browser()
  segments <- segments %>%
    dplyr::mutate(
      pixset.intercept = ifelse(theta == 0, xintercept, (height - yintercept) / slope),
      xaverage = ifelse(theta == 0, xintercept,
        ((0 - yintercept) / slope + (height - yintercept) / slope) / 2
      ),
      xbottom = ifelse(theta == 0, xintercept, (height - yintercept) / slope),
      xtop = xintercept
    )

  # Find the middle 2/3rds
  lthird <- width / 6
  uthird <- 5 * width / 6

  # Find best bottom and top index of groove for both sides
  top.left <- segments$xintercept[which.min(abs(segments$xintercept - lthird))]
  bottom.left <- segments$pixset.intercept[which.min(abs(segments$pixset.intercept - lthird))]


  top.right <- segments$xintercept[which.min(abs(segments$xintercept - uthird))]
  bottom.right <- segments$pixset.intercept[which.min(abs(segments$pixset.intercept - uthird))]

  # Calculate equation of line for each side
  slope.left <-  (top.left-bottom.left)/height
  slope.right <- (top.right - bottom.right)/height

  # Crate two functions to calculate the x output for each y input
  left_groove_fit <- function(yinput) {
    assert_that(is.numeric(yinput))

    bottom.micron <- pix_to_micron(bottom.left) # scale bottom.left to microns

    if (is.infinite(slope.left)) return(NA) # hough didn't find a groove

    if (length(slope.left) == 0) { # straight vertical line
      left.groove <- rep(bottom.micron, length(yinput)) + adjust
    } else {
      left.groove <- (bottom.micron + slope.left*yinput) + adjust
    }
    return(left.groove)
  }

  right_groove_fit <- function(yinput) {
    assert_that(is.numeric(yinput))
    bottom.micron <- pix_to_micron(bottom.right) # scale bottom.right to microns

    if (is.infinite(slope.right)) return(NA) # hough didn't find a groove

    if (length(slope.right) == 0) { # straight vertical line
      right.groove <- rep(bottom.micron, length(yinput)) - adjust
    } else {
      right.groove <- (bottom.micron + slope.right*yinput) - adjust
    }
    return(right.groove)
  }

  # summarize the land before visualizing
  land.summary <- dplyr::summarize(dplyr::group_by(land, x), value = median(value, na.rm = TRUE))


  return(list(left.groove.fit = left_groove_fit, right.groove.fit = right_groove_fit))
}
