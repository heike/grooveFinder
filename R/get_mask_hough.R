#' Create mask layer for grooves identified with hough
#'
#' Function that creates a mask layer for an x3p of a bullet land based on the predictions
#' of the Hough process.
#'
#' @param land.x3p an x3p file containing a bullet land, including surface matrix
#' @param grooves a list object containing the output of the get_grooves_hough function. Should include a left and right groove fit equations
#' @return an x3p file with the added mask that estimates the GEAs from the Hough method.
#' @export
#' @importFrom x3ptools x3p_get_scale x3p_add_mask
#' @examples
#' library(x3ptools)
#' data("br411")
#' x3p <- br411
#' grooves <- get_grooves_hough(x3p_to_df(x3p), qu = 0.995)
#'
#' a <- get_mask_hough(x3p, grooves)
#' \dontrun{x3ptools::image_x3p(a)}
get_mask_hough <- function(land.x3p, grooves){
  left <- grooves$left.groove.fit(0:(ncol(land.x3p$surface.matrix)-1)*x3p_get_scale(land.x3p))
  left <- floor(left/x3p_get_scale(land.x3p) + 1)
  right <- grooves$right.groove.fit(0:(ncol(land.x3p$surface.matrix)-1)*x3p_get_scale(land.x3p))
  right <- floor(right/x3p_get_scale(land.x3p) + 1)

  if (is.null(land.x3p$mask)){
     land.x3p <-  x3p_add_mask(land.x3p)
  }
  # mask <- land.x3p$mask

  # create mask. Masks have transposed size as surface matrix. :(
  mask <- matrix(data = FALSE,
                 nrow = ncol(land.x3p$surface.matrix),
                 ncol = nrow(land.x3p$surface.matrix))

  leftgroove <- sapply(1:length(left), FUN = function(i){
    mask[ i,1:floor(left[i])] <- TRUE
    mask[i,]
  }) %>%t()

  land.x3p <- x3ptools::x3p_add_mask_layer(land.x3p, mask = leftgroove, color="#c4221a", annotation = "Left groove")

  rightgroove <- sapply(1:length(right), FUN = function(i){
      mask[i, floor(right[i]):ncol(mask)] <- TRUE
      mask[i,]
  }) %>% t()

  land.x3p <- x3ptools::x3p_add_mask_layer(land.x3p, mask = rightgroove, color="#3279a8", annotation = "Right groove")

  return(land.x3p)
}

