#' @export
get_grooves_manual <- function() {
  appDir <- system.file("shiny-examples", "detectGrooves", package = "mypackage")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
