#' Check grooves for correctness
#'
#' @param x output from cc_locate_grooves
#' @return TRUE if ok, error otherwise
#' @importFrom assertthat assert_that
check_grooves <- function(x) {
  assert_that(has_name(x, "groove"))
  assert_that(is.numeric(x$groove))
  assert_that(length(x$groove) <= 2)
  return(TRUE)
}
