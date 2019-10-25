#' Full GWENA pipeline analysis
#'
#' Transform a network matrix into a graph vertices and edges
#'
#' @param modules_phenotype List form modules_phenotype function output
#' TODO finish
#'
#' @details
#' TODO
#' @return
#' TODO
#'
#' @examples
#' #TODO
#'
#' @import
#'
#' @export

GWENA_full <- function(data_expr,
                  # net_building
                  cor_func = c("pearson", "spearman", "bicor", "other"), your_func = NULL,
                  power_value = NULL, fit_cut_off = 0.90, network_type = c("unsigned", "signed", "signed hybrid", "none"),
                  tom_type = c("unsigned", "signed", "signed Nowick"), save_adjacency = FALSE,
                  # modules_detection
                  min_module_size = min(20, ncol(data_expr) / 2), merge_cut_height = 0.25,
                  # modules enrichment

                  #

                  # exec and report params
                  n_threads = 0, detailled_result = FALSE, ...) {

}

GWENA_simplified <- function() {

}
