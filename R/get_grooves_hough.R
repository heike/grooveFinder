

#' Re-parameterize the Hough transform output for pixel identification
#'
#' @param rho Numeric vector containing the shortest distance from the line to the origin
#' @param theta Numeric vector containing the angle of the line from the positive x axis
#' @param df Data frame containing output from a Hough transformation
#' @export


rho_to_ab <- function(rho = NULL, theta = NULL, df = NULL) {
  if (is.null(df)) {
    df <- data.frame(rho = rho, theta = theta)
  }
  stopifnot(c("rho", "theta") %in% names(df))
  idx <- df$theta ==0
  df <- df %>%
    mutate(yintercept = ifelse(idx, NA, rho/sin(theta)),
           slope = -cos(theta)/sin(theta),
           xintercept = ifelse(idx, rho, rho/cos(theta)))
  df
}



#' Use Hough transformation to identify groove locations.
#' Choose strong edges based on whether scores exceed the 99.75% percentile of scores and if the angle of line is less than
#' Pi/4.
#'
#'
#' @param land dataframe of surface measurements in microns in the x, y, and x direction
#' @param return_plot return plot of grooves
#'
#' @importFrom x3ptools df_to_x3p
#' @importFrom imager as.cimg
#' @importFrom imager imgradient
#' @importFrom imager hough_line
#' @importFrom assertthat assert_that
#' @importFrom assertthat has_name
#' @export

get_grooves_hough <- function(land, return_plot=F){
  assert_that(has_name(land, "x"), has_name(land, "y"), has_name(land, "value"),
              is.numeric(land$x), is.numeric(land$y), is.numeric(land$value))
  # Convert to cimage
  land.x3p <- df_to_x3p(land)

  # HH: this if condition distinguishes between lands based on the aspect ratio.
  # please assume that we have lands that are wider than high.
  # If they are not, spit out a warning rather than different code.
  # It has bitten us badly into rear elements before to try to think for the user.
  if(width(land.x3p$surface.matrix) < height(land.x3p$surface.matrix)){
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
  hough.df <- hough_line(strong, data.frame = TRUE)

  # Subset based on score and angle rotation

  hough.df <- hough.df %>%
    mutate(theta = ifelse(theta <= pi, theta, theta - 2*pi)) %>%
    filter(score > quantile(score, .999),
           theta > (-pi/4),
           theta < (pi/4))

  summary.save <- summary(hough.df)

  # Find only the good vertical segments
  segments <- rho_to_ab(df = hough.df)

  # Calculate the intercept of where each Hough line

  segments <- segments %>%
    mutate(pixset.intercept = ifelse(theta==0, xintercept, (height(strong) - yintercept)/slope),
           xaverage = ifelse(theta==0, xintercept, ((0-yintercept)/slope + (height(strong) - yintercept)/slope)/2))

  good_vertical_segs <- segments %>%
    extract2("xaverage")

  lthird <- width(strong)/6
  uthird <- 5*width(strong)/6

  # Find hough line index where
  closelthird <- good_vertical_segs[which.min(abs(good_vertical_segs - lthird))]
  closeuthird <- good_vertical_segs[which.min(abs(good_vertical_segs - uthird))]

  groove <- c(closelthird, closeuthird)

  # summarize the land before visualizing
  land.summary <- summarize(group_by(land, x), value = median(value, na.rm=TRUE))

  if(return_plot){
    return(
      list(
        groove,
        plot = grooves_plot(land = land, grooves = groove)
      )
    )
  }
  else{
    return(list(groove = groove))
  }

}

#' Return a list of errors and results from get_groove_hough
#'
#' @param land dataframe of surface measurements in microns in the x, y, and x direction
#' @param return_plot return plot of grooves
#'
#' @importFrom purrr safely
#' @export


safely_get_grooves_hough <- purrr::safely(get_grooves_hough)

