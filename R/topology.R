#' Return graph from squared matrix network
#'
#' Takes a squared matrix containing the pairwise similarity scores for each gene and return a igraph object.
#'
#' @param sq_mat matrix or data.frame, squared matrix representing
#'
#' @export

build_graph_from_sq_mat <- function(sq_mat) {
  # Checks
  if (!(is.matrix(sq_mat) || is.data.frame((sq_mat)))) stop("sq_mat should be a matrix or a data.frame")
  if (any(is.na(sq_mat))) warning("sq_mat should not contain any missing value")
  if ((any(sq_mat > 1) || any(sq_mat < -1)) && !any(is.na(sq_mat))) stop("sq_mat should be filled with value in the [-1,1] range")
  if (nrow(sq_mat) != ncol(sq_mat)) stop("sq_mat must be a squared matrix")

  # From matrix to graph
  sq_mat %>%
    as.data.frame %>%
    tibble::rownames_to_column("gene_from") %>%
    tidyr::pivot_longer(-gene_from, names_to = "gene_to", values_to = "weight") %>%
    .[!duplicated(t(apply(.[, c("gene_from", "gene_to")], 1, sort))), ] %>% # removing pairwise duplicates
    as.data.frame %>%
    dplyr::filter(gene_from != gene_to) %>% #removing self edge
    igraph::graph_from_data_frame(directed = FALSE)
}


#' Determine hub genes based on connectivity
#'
#' Compute connectivity of each gene by module if provided or for whole network if not, and return the top_n highest connected ones.
#'
#' @param network matrix or data.frame, square table representing connectivity between each genes as returned by
#' build_net. Can be whole network or a single module.
#' @param modules list, modules defined as list of gene vectors. If null, network is supposed to be the whole network
#' or an already split module
#' @param top_n integer, number of genes to be considered as hub genes
#'
#' @return list of vectors, or single vector of gene names
#'
#' @importFrom magrittr %>%
#'
#' @export

get_hub_high_co <- function(network, modules = NULL, top_n = 5) {
  # Checks
  .check_is_network(network)
  if (!is.null(modules)) {
    .check_is_module(modules, is.list(modules))}
  if (length(top_n) > 1) stop("top_n must be a single numeric value")
  if (!is.numeric(top_n)) stop("top_n must be a numeric value")
  if (top_n < 1 || top_n %% 1 != 0) stop("If not NULL, block_size must be a whole number >= 1")

  # Highly connected genes selection
  if (is.null(modules)) { # considered network is already a single module, or looking for hubs independently of modules split
    if (top_n > ncol(network)) stop("top_n > to number of genes into network")
    hubs <- rowSums(network) %>% sort(decreasing = TRUE) %>% .[1:top_n] %>% names
  } else {
    hubs <- lapply(modules, function(x){
      net <- network[x,x]
      if (top_n > ncol(net)) {
        warning("top_n > to number of genes into one of the modules. Taking all genes in this case.")
        top_n <- ncol(net)
      } else {
        hubs_name <- rowSums(net) %>% sort(decreasing = TRUE) %>% .[1:top_n] %>% names
      }
    })
  }
  return(hubs)
}

#' Determine hub genes based on degree
#'
#' Remove edges from the graph which value is under weight_th then compute degree of each node (gene).
#' Hub gene are genes whose degree value is above average degree value of the thresholded network.
#'
#' @param network matrix or data.frame, square table representing connectivity between each genes as returned by
#' build_net. Can be whole network or a single module.
#' @param modules list, modules defined as list of gene vectors. If null, network is supposed to be the whole network
#' or an already split module
#' @param weight_th decimal, weight threshold under or equal to which edges will be removed
#'
#' @detail GWENA natively build networks using WGCNA. These networks are complete in a graph theory sens, meaning all nodes
#' are connected to each other. Therefore a threshold need to be applied so degree of all nodes isn't the same ()
#'
#' @importFrom igraph
#'
#' @return list of vectors, or single vector of gene names
#'
#' @export

get_hub_degree <- function(network, modules, weight_th = 0.2) {
  # Checks
  if (!(is.data.frame(network) || is.matrix(network))) stop("network must be a data.frame or a matrix")
  if (is.null(colnames(network)) || is.null(rownames(network))) stop("network must have colnames and rownames")
  if (!is.null(modules)) {
    if (!is.list(modules)) stop("modules must be a list")
    if (is.null(names(modules))) stop("modules list must have names")
    if (any(names(modules) == "" %>% unlist)) stop("modules list must have names for all elements")
    if (any(lapply(modules, is.vector, "character") %>% unlist)) stop("modules list element must be vector of gene names")
  }

  # TODO : check if can avoid transforming to igraph object. It takes a lot of time...


  # Above average degree genes
  if (is.null(modules)) { # considered network is already a single module, or looking for hubs independently of modules split
    # Graph whole network
    net_degree <- get_graph_from_sq_mat(network) %>%
      igraph::delete.edges(which(igraph::E(.)$weight < weight_th)) %>%
      igraph::degree()
    hubs <- net_degree[which(net_degree > (net_degree %>% mean))] %>% names
  } else {
    hubs <- lapply(modules, function(x){
      net_degree <- network[x,x] %>%
        get_graph_from_sq_mat() %>%
        igraph::delete.edges(which(igraph::E(.)$weight < weight_th)) %>%
        igraph::degree()
      x_hubs <- net_degree[which(net_degree > (net_degree %>% mean))] %>% names
    })
  }
  return(hubs)
}


#' Determine hub genes inside each module
#'
#' Return genes considered as hub genes inside each module of a network
#'
#' @param network
#' @param modules
#' @param method
#'
#' @detail
#' \describe{
#'   \item{highest connectivity}{Select the top n (n depending on parameter given) highest connected genes. Similar to WGCNA::chooseTopHubInEachModule.}
#'   \item{superior degree}{Select genes which degree is greater than average connection degree of the network. Definition from network theory.}
#'   \item{Kleinberg's score}{Select genes which Kleinberg's score superior to provided threshold.}
#' }
#'
#' @return list of vectors representing hub genes, by module
#'
#' @export

get_hub_genes <- function(network, modules, method) {

}

# Idee pour fournir method avec nom method + param : des listes à passer avec un element nommé "name" puis les autres elements avec le
# nom des parametres et leur valeur
# Ex : method = list(name = "highest connectivity", n = 10)
