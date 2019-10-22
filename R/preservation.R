#' Checking preservation of modules between one or more conditions
#'
#' Detect which modules are preserved from a condition to another
#'
#' @param conditions List form modules_phenotype function output
#' TODO finish
#'
#' @details
#' TODO
#' @return
#' TODO
#'
#' @examples
#' TODO
#'
#' @import NetRep
#' @import dplyr
#'
#' @export

modules_preservation <- function(data_expr_list, net_list, modules_list, reference, test) {
  # Checks

  #
  preservation <- NetRep::modulePreservation(
    network = lapply(net_list, function(condition) condition$tom %>% as.matrix),
    data = lapply(data_expr_list, function(condition) condition %>% as.matrix),
    correlation = lapply(net_list, function(condition) condition$adjacency %>% as.matrix),
    moduleAssignments = modules_list[[1]]$modules, discovery = reference, test = test,
    nPerm = 10000
  ) # TODO Check if it's adj/tom/cor for network and/or cor
}
