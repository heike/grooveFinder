

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
#' Choose strong edges based on whether scores exceed the 99.75 percentile of scores and if the angle of line is less than
#' Pi/4.
#' @param land dataframe of surface measurements in microns in the x, y, and x direction
#' @param qu quantile (between 0 and 1) to specify score quantile for which vertical lines are considered. If groove are not strongly expressed, lower this threshold.
#' @param adjust positive number to adjust the grooves inward
#' @return list object consisting of functions to describe the left and right groove
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

get_grooves_hough <- function(land, qu = 0.999, adjust=10){
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
      xaverage = ifelse(theta==0, xintercept,
                        ((0-yintercept)/slope + (height(strong) - yintercept)/slope)/2),
      xbottom = ifelse(theta==0, xintercept, (height(strong) - yintercept)/slope),
      xtop = xintercept
    )

  # Find the middle 2/3rds

  lthird <- width(strong)/6
  uthird <- 5*width(strong)/6

  # Find best bottom and top index of groove for both sides
  top.left <- top.left <- segments$xintercept[which.min(abs(segments$xintercept - lthird))]
  bottom.left <- segments$pixset.intercept[which.min(abs(segments$pixset.intercept - lthird))]

  top.right <- segments$xintercept[which.min(abs(segments$xintercept - uthird))]
  bottom.right <- segments$pixset.intercept[which.min(abs(segments$pixset.intercept - uthird))]

  # Calculate equation of line for each side
  slope.left <- -height(strong)/(top.left - bottom.left)
  yint.left <- (-(slope.left*top.left))*x3p_get_scale(land.x3p)

  slope.right <- -height(strong)/(top.right - bottom.right)
  yint.right <- (-(slope.right*top.right))*x3p_get_scale(land.x3p)

  #Crate two functions to calculate the x output for each y input
  left_groove_fit <- function(yinput){
    assert_that(is.numeric(yinput))

    if(is.infinite(slope.left)){
      left.groove <- bottom.left
    }

    else{
      # Do I need the scaled +1 to adjust location of the groove?
      left.groove <- ((yinput - yint.left)/slope.left)
    }
    return(left.groove)
  }

  right_groove_fit <- function(yinput){
    assert_that(is.numeric(yinput))

    if(is.infinite(slope.right)){
      right.groove <- bottom.right
    }

    else{
      right.groove <- ((yinput - yint.right)/slope.right)
    }
    return(right.groove)
  }

  # summarize the land before visualizing
  land.summary <- dplyr::summarize(dplyr::group_by(land, x), value = median(value, na.rm=TRUE))


  return(list(left.groove.fit = left_groove_fit, right.groove.fit = right_groove_fit))
}





