# ==== Format ====

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


# ==== Checks ====

#' Determine if an object is a network
#'
#' Check content of a given object to determine if it's a network, meaning a squared matrix of similarity
#' score between genes
#'
#' @param network matrix or data.frame, object to test to be a network
#'
#' @return list, a boolean as first element and in second element NULL or the reason why boolean is
#' set to FALSE
#'
#' @export

is_network <- function(network) {
  if (!(is.data.frame(network) || is.matrix(network))) return(list(bool = FALSE, reason = "network must be a data.frame or a matrix"))
  if (is.null(colnames(network)) || is.null(rownames(network))) return(list(bool = FALSE, reason = "network must have colnames and rownames"))
  if (!all(colnames(network) %in% rownames(network))) return(list(bool = FALSE, reason = "colnames and rownames form network doesn't match"))
  if (ncol(network) != nrow(network)) return(list(bool = FALSE, reason = "network must be a squared matrix"))
  if ((any(network > 1) || any(network < -1)) && !any(is.na(network))) return(list(bool = FALSE, reason = "network should be filled with value in the [-1,1] range"))
  #else
  return(list(bool = TRUE, reason = NULL))
}

#' Check if an object is a network
#'
#' Check content of a given object to determine if it's a network, meaning a squared matrix of similarity
#' score between genes.
#'
#' @param network matrix or data.frame, object to test to be a network
#'
#' @return Throw an error if doesn't correspond

.check_is_network <- function(network) {
  check <- is_network(network)
  if (!check$bool) {
    stop(check$reason)
  }
}


#' Determine if an object is a module or a list of modules
#'
#' Check content of a given object to determine if it's a module or a list of modules, meaning a single
#' vector of characters which are gene names, or a named list of these vectors
#'
#' @param module vector or list, object to test to be a module or list of modules
#' @param is_list boolean, indicate if module must be tested as a single module or a list of modules
#'
#' @return list, a boolean as first element and in second element NULL or the reason why boolean is
#' set to FALSE
#'
#' @export

is_module <- function(module, is_list = FALSE) {
  if (is_list) {
    if (!is.list(module)) return(list(bool = FALSE, reason = "module must be a list"))
    if (is.null(names(module))) return(list(bool = FALSE, reason = "module list must have names"))
    if (any(names(module) == "" %>% unlist)) return(list(bool = FALSE, reason = "module list must have names for all elements"))
    if (any(lapply(module, is.vector, "character") %>% unlist %>% `!`)) return(list(bool = FALSE, reason = "module list element must be vector of gene names"))
  } else {
    if (!is.vector(module, "character")) return(list(bool = FALSE, reason = "module must be vector of gene names"))
  }
  return(list(bool = TRUE, reason = NULL))
}


#' Check if an object is a module or a list of modules
#'
#' Check content of a given object to determine if it's a module or a list of modules, meaning a single
#' vector of characters which are gene names, or a named list of these vectors
#'
#' @param module vector or list, object to test to be a module or list of modules
#' @param is_list boolean, indicate if module must be tested as a single module or a list of modules
#'
#' @return Throw an error if doesn't correspond

.check_is_module <- function(module, is_list = FALSE) {
  check <- is_module(module, is_list)
  if (!check$bool) {
    stop(check$reason)
  }
}

