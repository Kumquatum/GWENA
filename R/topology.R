#' Determine hub genes based on connectivity
#'
#' Compute connectivity of each gene by module if provided or for whole network if not and return the top_n highest connected ones.
#'
#' @param network matrix or data.frame, square table representing connectivity between each genes as returned by
#' build_net. Can be whole network or a single module.
#' @param modules list, modules defined as list of gene vectors. If null, network is supposed to be the whole network
#' or an already split module
#' @param top_n integer, number of highest connected genes to be considered as hub genes
#'
#' @return list of vectors, or single vector of gene names

get_hub_gene_high_co <- function(network, modules = NULL, top_n = 5) {
  # Checks
  if (!(is.data.frame(network) || is.matrix(network))) stop("network must be a data.frame or a matrix")
  if (is.null(colnames(network)) || is.null(rownames(network))) stop("network must have colnames and rownames")
  if (!is.null(modules)) {
    if (!is.list(modules)) stop("modules must be a list")
    if (is.null(names(modules))) stop("modules list must have names")
    if (any(names(modules) == "" %>% unlist)) stop("modules list must have names for all elements")
    if (any(lapply(modules, is.vector, "character") %>% unlist)) stop("modules list element must be vector of gene names")
  }
  if (!is.numeric(top_n)) stop("top_n must be an integer")
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
#'   \item{local and connector}{}
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
