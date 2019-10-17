#' Match a correlation function based on a name
#'
#' Translate a function name into a usable function in code.
#' TODO finish
#'
#' @param cor_name string of the name of the correlation to be use
#'
#' @return A function corresponding to the correlation required
#' @author  Gwenaëlle Lemoine <lemoine.gwenaelle@@gmail.com>
#' @examples
#'
#' @import WGCNA

cor_func_match <- function(cor_func = c("pearson", "spearman", "bicor")){
  # Checks
  if (is.null(cor_func)) stop("cor_func must be a character vector representing one of the supported correlation functions. See ?cor_func_match for more information.")
  cor_func <- match.arg(cor_func)

  # Function assignation
  if (cor_func == "pearson") {
    cor_func <- function(x) WGCNA::cor(x, method = "pearson", use = "e")
  } else if (cor_func == "spearman") {
    cor_func <- function(x) stats::cor(x, method = "spearman", use = "e")
  } else if (cor_func == "bicor") {
    cor_func <- function(x) WGCNA::bicor(x, use = "e")
  } else {
    stop("Should never be triggered")
  }
}
