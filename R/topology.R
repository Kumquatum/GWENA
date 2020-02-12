#' Return graph from squared matrix network
#'
#' Takes a squared matrix containing the pairwise similarity scores for each gene and return a igraph object.
#'
#' @param sq_mat matrix or data.frame, squared matrix representing
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
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
    modules <- list(colnames(network))
  }

  hubs <- lapply(modules, function(x){
    net <- network[x,x]
    if (top_n > ncol(net)) {
      warning("top_n > to number of genes into one of the modules. Taking all genes in this case.")
      top_n <- ncol(net)
    } else {
      hubs_name <- rowSums(net) %>% sort(decreasing = TRUE) %>% .[1:top_n]
    }
  })
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
#' @importFrom igraph delete.edges degree E
#'
#' @return list of vectors, or single vector of gene names
#'
#' @export

get_hub_degree <- function(network, modules = NULL, weight_th = 0.2) {
  # Checks
  .check_is_network(network)
  if (!is.null(modules)) {
    .check_is_module(modules, is.list(modules))}
  if (!is.null(weight_th)) {
    if (length(weight_th) > 1) stop("weight_th must be a single numeric value")
    if (!is.numeric(weight_th)) stop("weight_th must be a numeric value")
    if (weight_th < 0 || weight_th >= 1) stop("weight_th must be a in [0;1[")
  }
  # TODO : check if can avoid transforming to igraph object. It takes a lot of time...

  # Above average degree genes
  if (is.null(modules)) { # considered network is already a single module, or looking for hubs independently of modules split
    modules <- list(colnames(network))
  }
  hubs <- lapply(modules, function(x){
    net_degree <- network[x,x] %>%
      build_graph_from_sq_mat() %>%
      igraph::delete.edges(which(igraph::E(.)$weight <= weight_th)) %>%
      igraph::degree()
    if (table(net_degree) %>% length == 1) {
      x_hubs <- c()
    } else {
      x_hubs <- net_degree[which(net_degree > (net_degree %>% mean))] }

  })

  return(hubs)
}


#' Determine hub genes based on Kleinberg's score
#'
#' Compute Kleinberg's score (defined as the principal eigenvector of A*t(A), where A is the similarity matrix of the graph)
#' of each gene by module if provided or for whole network if not, and return the top_n highest ones.
#'
#' @param network matrix or data.frame, square table representing connectivity between each genes as returned by
#' build_net. Can be whole network or a single module.
#' @param modules list, modules defined as list of gene vectors. If null, network is supposed to be the whole network
#' or an already split module
#' @param top_n integer, number genes to be considered as hub genes
#' @param k_th decimal, Kleinberg's score threshold above or equal to which genes are considered as hubs
#'
#' @details If you provide a top_n value, you can't provide a k_th value and vice versa. If none of them is provided
#' top_n = 5.
#' For more information on Kleinberg's score, look at \code{\link[igraph{hub_score}]} from igraph.
#'
#' @importFrom igraph hub_score
#'
#' @return list of vectors, or single vector of gene names
#'
#' @export

get_hub_kleinberg <- function(network, modules = NULL, top_n = NULL, k_th = NULL) {
  # Checks
  .check_is_network(network)
  if (!is.null(modules)) {
    .check_is_module(modules, is.list(modules))}
  if (is.null(top_n) && is.null(k_th)) {
    top_n <- 5
    warning("No top_n or k_th value provided. Default: top_n = ", top_n) }
  if (!is.null(top_n) && !is.null(k_th)) {
    k_th <- NULL
    warning("Conflict: top_n and k_th value provided. Keeping only top_n value") }
  if (!is.null(k_th)) {
    if (length(k_th) > 1) stop("k_th must be a single numeric value")
    if (!is.numeric(k_th)) stop("k_th must be a numeric value")
    if (k_th < 0 || k_th >= 1) stop("k_th must be a in [0;1[")
  }
  if (!is.null(top_n)) {
    if (length(top_n) > 1) stop("top_n must be a single numeric value")
    if (!is.numeric(top_n)) stop("top_n must be a numeric value")
    if (top_n < 1 || top_n %% 1 != 0) stop("If not NULL, block_size must be a whole number >= 1")
  }

  if (is.null(modules)) { # considered network is already a single module, or looking for hubs independently of modules split
    modules <- list(module = colnames(network))
  }
  hubs <- lapply(modules, function(x){
    net_hub_score <- network[x,x] %>%
      build_graph_from_sq_mat() %>%
      igraph::hub_score(scale = TRUE) %>%
      .$vector

    # With top_n
    if (!is.null(top_n)) {
      x_hubs <- net_hub_score %>% sort(decreasing = TRUE) %>% .[1:top_n]
    } else { # With k_th
      x_hubs <- net_hub_score[which(net_hub_score >= k_th)]
    }
    return(x_hubs)
  })

  return(hubs)
}


#' Determine hub genes inside each module
#'
#' Return genes considered as hub genes inside each module of a network following the selected method. Method will be lauched with default
#' parameters. If specific parameters desired, please use directly the function \code{get_hub_...} itself.
#'
#' @param network matrix or data.frame, square table representing connectivity between each genes as returned by
#' build_net. Can be whole network or a single module.
#' @param modules list, modules defined as list of gene vectors. If null, network is supposed to be the whole network
#' or an already split module
#' @param methodstring, name of the method to be used for hub gene detection. See details.
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

get_hub_genes <- function(network, modules, method = c("highest connectivity", "superior degree", "Kleinberg's score")) {
  # Checks
  method = match.arg(method)

  # Call
  if (method == "highest connectivity") { hubs <- get_hub_high_co(network, modules)
  } else if (method == "superior degree") { hubs <- get_hub_degree(network, modules)
  } else if (method == "Kleinberg's score") { hubs <- get_hub_kleinberg(network, modules) }

  return(hubs)

}







