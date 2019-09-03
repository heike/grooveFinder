#' @examples
#'
#' x3p <- read_x3p("")
#' grooves <-get_grooves_hough(x3p_to_df(x3p), qu = 0.99)
#'
#' get_mask_hough(x3p, grooves)
#'
get_mask_hough <- function(land.x3p, grooves){
  left <- grooves$left.groove.fit(1:ncol(land.x3p$surface.matrix))/x3p_get_scale(land.x3p)
  right <- grooves$right.groove.fit(1:ncol(land.x3p$surface.matrix))/x3p_get_scale(land.x3p)


  if(is.null(land.x3p$mask)){
    land.x3p <- land.x3p %>%
      x3p_add_mask()
  }

  mask <- raster::as.matrix(land.x3p$mask)
  mask <- sapply(1:length(left), FUN = function(i){
    mask[ i,1:floor(left[i])] <- "#c4221a"
    mask[i,]
  }) %>%t()

  mask <- sapply(1:length(right), FUN = function(i){
      mask[i, floor(right[i]):ncol(mask)] <- "#3279a8"
      mask[i,]
  }) %>% t()

  land.x3p <- land.x3p %>% x3p_add_mask(mask = mask)


}

