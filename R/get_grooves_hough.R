

#' Re-parameterize the Hough transform output for pixel identification
#'
#' Helper function, that should not be used outside the package.
#' Computes x and y intercepts for a line given in polar coordinate representation.
#' A line in polar coordinate representation is given by angle theta (with respect to the x-axis) and
#' distance of the line to the origin rho.
#' @param rho Numeric vector containing the shortest distance from the line to the origin
#' @param theta Numeric vector containing the angle of the line from the positive x axis
#' @param df Data frame containing output from a Hough transformation
#' @return data frame with variables rho, theta, score (original data frame) expanded by yintercept, xintercept and slope.
rho_to_ab <- function(rho = NULL, theta = NULL, df = NULL) {
  if (is.null(df)) {
    df <- data.frame(rho = rho, theta = theta)
  }
  stopifnot(c("rho", "theta") %in% names(df))

  df <- df %>%
    mutate(
      yintercept = ifelse(theta == 0, NA, rho/sin(theta)),
      slope = -cos(theta)/sin(theta),
      xintercept = rho/cos(theta)) # cos(theta) == 0 for theta = ± pi/2 but realistically we don't get there
  df
}



#' Use Hough transformation to identify groove locations.
#'
#' Choose strong edges based on whether scores exceed the 99.75% percentile of scores and if the angle of line is less than
#' Pi/4.
#' @param land dataframe of surface measurements in microns in the x, y, and x direction
#' @param qu quantile (between 0 and 1) to specify score quantile for which vertical lines are considered. If groove are not strongly expressed, lower this threshold.
#' @param adjust positive number to adjust the grooves inward
#' @param return_plot return plot of grooves
#' @return list object consisting of a vector of groove values (left and right) and, if return_plot is TRUE, a plot of the profile with the groove locations
#'
#' @importFrom x3ptools df_to_x3p
#' @importFrom imager as.cimg width height
#' @importFrom imager imgradient
#' @importFrom imager hough_line
#' @importFrom assertthat assert_that
#' @importFrom assertthat has_name
#' @importFrom stats quantile median sd na.omit
#' @importFrom dplyr filter mutate group_by summarize count
#' @importFrom x3ptools x3p_get_scale df_to_x3p
#' @export
#' @examples
get_grooves_hough <- function(land, qu = 0.999, adjust=10, return_plot=F){
  assert_that(has_name(land, "x"), has_name(land, "y"), has_name(land, "value"),
              is.numeric(land$x), is.numeric(land$y), is.numeric(land$value))
  # Convert to cimage
  land.x3p <- df_to_x3p(land)

  # visible bindings problem
  theta <- score <- rho <- xintercept <- yintercept <- 0
  slope <- x <- value <- 0

  # HH: this if condition distinguishes between lands based on the aspect ratio.
  # please assume that we have lands that are wider than high.
  # If they are not, spit out a warning rather than different code.
  # It has bitten us badly into rear elements before to try to think for the user.
  sizes <- dim(land.x3p$surface.matrix)
  width <- sizes[1]
  height <- sizes[2]
  if(width < height){
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

  strong <- grad.mag > quantile(grad.mag,.99, na.rm=TRUE)
  # create the hough transform
  hough.df <- hough_line(strong, data.frame = TRUE, shift = FALSE) # we want to get values with respect to (0,0) not (1,1)
  assert_that(
    has_name(hough.df, "theta"),
    has_name(hough.df, "rho"),
    has_name(hough.df, "score"))

  # Subset based on score and angle rotation
  hough.df <- hough.df %>%
    dplyr::mutate(theta = ifelse(theta <= pi, theta, theta - 2*pi)) %>%
    dplyr::filter(score > quantile(score, qu),
           theta > (-pi/16), # identify only vertical(ish) lines
           theta < (pi/16),
           (rho < abs(width(strong))*1/6 | rho > width(strong)*5/6) # at either end of the LEA
           )

  summary.save <- summary(hough.df)

  # get x and y intercepts  (in pixel dimensions)
  segments <- rho_to_ab(df = hough.df)
  assert_that(
    has_name(segments, "theta"),
    has_name(segments, "rho"),
    has_name(segments, "score"),
    has_name(segments, "xintercept"),
    has_name(segments, "yintercept"),
    has_name(segments, "slope")
  )

# browser()
  segments <- segments %>%
    dplyr::mutate(
      pixset.intercept = ifelse(theta==0, xintercept, (height(strong) - yintercept)/slope),
      xaverage = ifelse(theta==0, xintercept, ((0-yintercept)/slope + (height(strong) - yintercept)/slope)/2))

  good_vertical_segs <- segments$xaverage

  lthird <- width(strong)/6
  uthird <- 5*width(strong)/6

  # Find hough line index where
  closelthird <- good_vertical_segs[which.min(abs(good_vertical_segs - lthird))]
  closeuthird <- good_vertical_segs[which.min(abs(good_vertical_segs - uthird))]

  groove <- c(closelthird, closeuthird) + adjust*c(1,-1) # adjust locations inward from steep drop-off
  groove <- groove *x3p_get_scale(land.x3p) # change from image width to locations in microns

  # summarize the land before visualizing
  land.summary <- dplyr::summarize(dplyr::group_by(land, x), value = median(value, na.rm=TRUE))

  if(return_plot){
    return(
      list(
        groove,
        plot = grooves_plot(land = land.summary, grooves = groove)
      )
    )
  }
  else{
    return(list(groove = groove))
  }

}

