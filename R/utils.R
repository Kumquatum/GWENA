#' Muting a function
#'
#' Prevent a function to output multiple message.
#' Source
#'
#' @param func Function who need to be muted.
#'
#' @return Nothing, just mute the called function

quiet <- function(func) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(func))
}
