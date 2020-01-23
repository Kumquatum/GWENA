#' Muting a function
#'
#' Prevent a function to output multiple message.
#' Source https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
#'
#' @param func Function who need to be muted.
#'
#' @return Nothing, just mute the called function

quiet <- function(func) {
  sink(tempfile())
  on.exit(sink())
  suppressMessages(invisible(force(func)))
}
