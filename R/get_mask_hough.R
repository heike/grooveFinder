#' @examples
#'
#' x3p <- read_x3p("")
#' grooves <-get_grooves_hough(x3p_to_df(x3p))
#'
#' get_mask_hough(x3p, grooves)
#'
get_mask_hough <- function(land.x3p, grooves){
  left <- grooves$left.groove.fit(1:ncol(land.x3p$surface.matrix))
  right <- grooves$right.groove.fit(1:ncol(land.x3p$surface.matrix))


  if(is.null(land.x3p$mask)){
    land.x3p <- land.x3p %>%
      x3p_add_mask()
  }

  mask <- raster::as.matrix(land.x3p$mask)
  mask <- sapply(1:length(left), FUN = function(i){
    mask[ i,1:floor(left[i])] <- 1
    mask[i,]
  })


  land.x3p <- land.x3p %>% x3p_add_mask_layer(mask, color = "red", annotation = "Left groove")


  browser()

}

