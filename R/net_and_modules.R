#' Match a correlation function based on a name
#'
#' Translate a function name into an R function.
#'
#' @param cor_func string of the name of the correlation to be use
#'
#' @return A function corresponding to the correlation required
#'
#' @importFrom WGCNA bicor cor

.cor_func_match <- function(cor_func = c("pearson", "spearman", "bicor")){
  # Checks
  if (is.null(cor_func))
    stop("cor_func must be a character vector representing one of the",
         " supported correlation functions.")
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


# Removing errors about dplyr data-variables
utils::globalVariables(c("Power", "SFT.R.sq"))

#' Calculating best fit of a power low on correlation matrix computed on
#' expression data
#'
#' Adjust a correlation matrix depending of the type of network, then try to
#' parameter a power law for best fit
#'
#' @param cor_mat matrix or data.frame of genes correlation.
#' @param fit_cut_off float, cut off by which R^2 (coefficient of
#' determination) will be thresholded. Must be in ]0;1[.
#' @param network_type string giving type of network to be used. Either
#' "unsigned", "signed", "signed hybrid". See details.
#' @param block_size integer giving size of blocks by which operations can be
#' proceed. Helping if working with low capacity computers. If null, will be
#' estimated.
#' @param ... any other parameter compatible with
#' \code{\link{pickSoftThreshold.fromSimilarity}}
#'
#' @details
#' network_type indicate which transformation will be applied on the
#' correlation matrix to return the similarity score.
#' \describe{
#'   \item{signed}{will modify the range [-1;1] to [0.5;1.5] (because of log10
#'   beeing used for scale free index computation)}
#'   \item{unsigned}{will return absolute value (moving from [-1;1] to [0;1])}
#'   \item{signed hybrid}{will replace all negative values by 0 (moving from
#'   [-1;1] to [0;1])}
#' }
#'
#' @return A list containing power of the law for best fit, fit table, and
#' metadata about the arguments used.
#' @examples
#' get_fit.cor(cor_mat = cor(kuehne_expr[, seq_len(100)]))
#'
#' @importFrom WGCNA pickSoftThreshold.fromSimilarity
#' @importFrom magrittr %>%
#' @importFrom dplyr select top_n
#'
#' @export

get_fit.cor <- function(cor_mat, fit_cut_off = 0.90, network_type =
                          c("unsigned", "signed", "signed hybrid"),
                        block_size = NULL, ...){
  # Checking args
  if (!is.matrix(cor_mat)) stop("cor_mat should be a matrix")
  if (any(is.na(cor_mat)))
    warning("cor_mat should not contain any missing value. It may be due to",
            " non variating probe/transcript when you computed your",
            " correlation matrix.")
  if ((any(cor_mat > 1) | any(cor_mat < -1)) & !any(is.na(cor_mat)))
    stop("cor_mat should be filled with value in the [-1,1] range")
  if (nrow(cor_mat) != ncol(cor_mat))
    stop("cor_mat should be a squared matrix")
  if (length(fit_cut_off) != 1 | !is.numeric(fit_cut_off))
    stop("power_cut_off should be a single number")
  if (fit_cut_off < 0 | fit_cut_off > 1) stop("power_cut_off should be a",
                                              " number between 0 and 1")
  network_type <- match.arg(network_type)
  if (!is.null(block_size)) {
    if (block_size < 2 | block_size %% 1 != 0)
      stop("If not NULL, block_size must be a whole number > 1")}

  # Calculating similarity
  similarity <- matrix()
  if (network_type == "unsigned") {
    similarity <- abs(cor_mat)
  } else if (network_type == "signed") {
    similarity <- (1 + cor_mat) / 2
  } else if (network_type == "signed hybrid") {
    cor_mat[cor_mat < 0] <- 0
    similarity <- cor_mat
  }

  # Getting fit
  sft_fit <- quiet(
    WGCNA::pickSoftThreshold.fromSimilarity(similarity = similarity,
                                            RsquaredCut = fit_cut_off,
                                            blockSize = block_size, ...))
  fit_above_cut_off <- TRUE
  if (is.na(sft_fit$powerEstimate)) { # If no fit, taking maximum fitting power
    warning("No fitting power could be found for provided fit_cut_off. Taking",
            " power for maximum fit. See FAQ for known causes.")
    sft_fit$powerEstimate <- sft_fit$fitIndice %>%
      dplyr::top_n(1, SFT.R.sq) %>% dplyr::select(Power) %>% as.numeric()
    fit_above_cut_off <- FALSE
  }

  # Final list with all infos
  fit <- list(
    power_value = sft_fit$powerEstimate,
    fit_table = sft_fit$fitIndice,
    fit_above_cut_off = fit_above_cut_off,
    metadata = list(
      network_type = network_type
    )
  )

  return(fit)
}


#' Calculating best fit of a power low on expression data
#'
#' Computes correlation matrix of the gene expression data, adjust it depending
#' of the type of network, then try to parameter a power law for best fit
#'
#' @param data_expr matrix or data.frame or SummarizedExperiment, expression
#' data with genes as column and samples as row.
#' @param fit_cut_off float, cut off by which R^2 (coefficient of
#' determination) will be thresholded. Must be in ]0;1[.
#' @param cor_func string specifying correlation function to be used. Must be
#' one of "pearson", "spearman", "bicor", "other". If "other", your_func must
#' be provided
#' @param your_func function returning correlation values. Final values must be
#' in [-1;1]
#' @param network_type string giving type of network to be used. Either
#' "unsigned", "signed", "signed hybrid". See details.
#' @param block_size integer giving size of blocks by which operations can be
#' proceed. Helping if working with low capacity computers. If null, will be
#' estimated.
#' @param ... any other parameter compatible with
#' \code{\link{pickSoftThreshold.fromSimilarity}}
#'
#' @details
#' network_type indicate which transformation will be applied on the
#' correlation matrix to return the similarity score.
#' \describe{
#'   \item{signed}{will modify the range [-1;1] to [0.5;1.5] (because of log10
#'   beeing used for scale free index computation)}
#'   \item{unsigned}{will return absolute value (moving from [-1;1] to [0;1])}
#'   \item{signed hybrid}{will replace all negative values by 0 (moving from
#'   [-1;1] to [0;1])}
#' }
#'
#' @return A list containing power of the law for best fit, fit table, and
#' metadata about the arguments used.
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#'
#' @examples
#' get_fit.expr(kuehne_expr[, seq_len(100)])
#'
#' @export

get_fit.expr <- function(data_expr, fit_cut_off = 0.90,
                         cor_func = c("pearson", "spearman", "bicor", "other"),
                         your_func = NULL, network_type =
                           c("unsigned", "signed", "signed hybrid"),
                         block_size = NULL, ...){
  # Checking args
  if (methods::is(data_expr, "SummarizedExperiment")) {
    data_expr <- t(SummarizedExperiment::assay(data_expr))
  } else .check_data_expr(data_expr)
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" & (is.null(your_func) | !is.function(your_func)))
    stop("If you specify other, your_func must be a function.")

  # Cor selection
  if (cor_func == "other") {
    cor_to_use <- your_func
  } else {
    cor_to_use <- .cor_func_match(cor_func)
  }

  # Calculating correlation matrix
  cor_mat <- cor_to_use(data_expr %>% as.matrix)

  # If personnal function, checking it returns a matrix with values in [-1;1]
  if (cor_func == "other"){
    if (min(cor_mat) < -1 | max(cor_mat) > 1)
      stop("your_func must be a function which returns values in [-1;1]")
  }

  # Getting fit
  fit <- get_fit.cor(cor_mat = cor_mat, fit_cut_off = fit_cut_off,
                     network_type = network_type, block_size = block_size, ...)

  # Final list with all infos
  fit$metadata$cor_func = ifelse(
    cor_func == "other",
    paste("other:", deparse(substitute(your_func))),
    cor_func)

  return(fit)
}


#' Network building by co-expression score computation
#'
#' Compute the adjacency matrix, then the TOM to build the network. Than detect
#' the modules by hierarchical clustering and thresholding
#'
#' @param data_expr matrix or data.frame or SummarizedExperiment, expression
#' data with genes as column and samples as row.
#' @param fit_cut_off float, cut off by which R^2 (coefficient of
#' determination) will be thresholded. Must be in ]0;1[.
#' @param cor_func string, name of the correlation function to be used. Must be
#' one of "pearson", "spearman", "bicor", "other". If "other", your_func must
#' be provided
#' @param your_func function returning correlation values. Final values must be
#' in [-1;1]
#' @param keep_matrices string, matrices to keep in final object. Can be one of
#' "none", "cor", "adja", "both". It is usefull to keep both if you plant to use
#' \code{\link{compare_conditions}}.
#' @param power_value integer, power to be applied to the adjacency matrix. If
#' NULL, will be estimated by trying different power law fitting.
#' @param block_size integer, size of blocks by which operations can be
#' proceed. Helping if working with low capacity computers. If null, will be
#' estimated.
#' @param stop_if_no_fit boolean, does not finding a fit above fit_cut_off
#' should stop process, or just print a warning and return the highest fitting
#' power.
#' @param network_type string, type of network to be used. Either "unsigned",
#' "signed", "signed hybrid". See details.
#' @param tom_type string, type of the topological overlap matrix to be
#' computed. Either "none", "unsigned", "signed", "signed Nowick",
#' "unsigned 2", "signed 2" and "signed Nowick 2". See detail at
#' \code{\link[WGCNA]{TOMsimilarityFromExpr}}.
#' @param n_threads integer, number of threads that can be used to paralellise
#' the computing
#' @param ... any other parameter compatible with
#' \code{\link{adjacency.fromSimilarity}}
#'
#' @return list containing network matrix, metadata of input parameters and
#' power fitting information.
#'
#' @importFrom WGCNA adjacency.fromSimilarity TOMsimilarity
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#'
#' @examples
#' net <- build_net(kuehne_expr[, seq_len(350)], n_threads = 1)
#'
#' @export

build_net <- function(data_expr, fit_cut_off = 0.90, cor_func =
                        c("pearson", "spearman", "bicor", "other"),
                      your_func = NULL, power_value = NULL, block_size = NULL,
                      stop_if_no_fit = FALSE,
                      network_type = c("unsigned", "signed", "signed hybrid"),
                      tom_type = c("unsigned", "signed", "signed Nowick",
                                   "unsigned 2", "signed 2", "none"),
                      keep_matrices = c("none", "cor", "adja", "both"),
                      n_threads = NULL, ...)
                      # TODO program the mclapply version
{
  # Checking
  if (methods::is(data_expr, "SummarizedExperiment")) {
    data_expr <- t(SummarizedExperiment::assay(data_expr))
  } else .check_data_expr(data_expr)
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" & (is.null(your_func) | !is.function(your_func)))
    stop("If you specify other, your_func must be a function.")
  if (!is.null(power_value)) {
    if (power_value < 1 | power_value %% 1 != 0)
      stop("If not NULL, power_value must be a whole number >= 1.")}
  if (!is.logical(stop_if_no_fit))
    stop("stop_if_no_fit must be a boolean.")
  network_type <- match.arg(network_type)
  tom_type <- match.arg(tom_type)
  keep_matrices <- match.arg(keep_matrices)
  if (!is.null(n_threads)) {
    if (!is.numeric(n_threads) | n_threads < 1 | n_threads %% 1 != 0)
      stop("n_threads must be a whole number >= 1")}

  # WGCNA's multi-threading
  if (n_threads == 1) {
    quiet(WGCNA::disableWGCNAThreads())
  } else quiet(WGCNA::enableWGCNAThreads(n_threads))
  # TODO : change this function for an internal one because WGCNA's function
  # isn't prefixed, causing warnings

  # Correlation selection and correlation matrix computation
  if (cor_func == "other") {
    cor_to_use <- your_func
  } else {
    cor_to_use <- .cor_func_match(cor_func)
  }
  cor_mat <- cor_to_use(data_expr %>% as.matrix)
  if (min(cor_mat < -1) | max(cor_mat) > 1)
    stop("Provided correlation function returned values outside [-1,1].")

  # Getting power
  if (is.null(power_value)) {
    fit <- get_fit.cor(cor_mat = cor_mat, fit_cut_off = fit_cut_off,
                       network_type = network_type, ...)
    if (stop_if_no_fit & fit$fit_above_cut_off == FALSE)
      stop("No fitting power could be found for provided fit_cut_off.",
           "You should verify your data (or lower fit_cut_off). See FAQ.")
  } else {
    fit <- list(
      power_value = power_value,
      fit_table = "None. Custom power_value provided")
    }

  # Adjacency
  adj = WGCNA::adjacency.fromSimilarity(similarity = cor_mat,
                                        type = network_type,
                                        power = fit$power_value)
  # Outputed matrix class is AsIs. Changing it to avoid later side effects
  class(adj) <- c("matrix", "array")

  # Topological overlap matrix
  if (tom_type != "none") {
    tom <- 1 - quiet(WGCNA::TOMsimilarity(adj, TOMType = tom_type)) %>%
      magrittr::set_colnames(colnames(adj)) %>%
      magrittr::set_rownames(rownames(adj))
  } else { tom <- adj }

  net = list(
    network = tom,
    metadata = list(
      cor_func = cor_func,
      network_type = network_type,
      tom_type = tom_type,
      power = fit$power_value,
      fit_power_table = fit$fit_table
    )
  )

  # Adding matrices to keep to the list
  if (keep_matrices %in% c("cor", "both"))
    net$cor_mat <- cor_mat
  if (keep_matrices %in% c("adja", "both"))
    net$adja_mat <- adj

  return(net)
}


# Removing errors about dplyr data-variables
utils::globalVariables(c("", ""))

#' Modules detection in a network
#'
#' Detect the modules by hierarchical clustering .
#'
#' @param data_expr matrix or data.frame or SummarizedExperiment, expression
#' data with genes as column and samples as row.
#' @param network matrix or data.frame, strengh of gene co-expression (edge
#' values).
#' @param min_module_size integer, lowest number of gene allowed in a module.
#' If none provided, estimated.
#' @param clustering_th float, threshold to be used by the clustering method.
#' For now \code{\link[dynamicTreeCut]{cutreeDynamic}}.
#' @param merge_close_modules boolean, does closest modules (based on
#' eigengene) should be merged together.
#' @param merge_threshold float, eigengenes correlation value over which
#' close modules will be merged. Must be in ]0;1[. See
#' \code{\link[WGCNA]{mergeCloseModules}}
#' @param detailled_result boolean, does pre-merge modules (if applicable)
#' and dendrogram included in output.
#' @param pam_respects_dendro boolean, If TRUE, the Partitioning Around Medoids
#' (PAM) stage will respect the dendrogram in the
#' sense that objects and small clusters will only be assigned to clusters that
#'  belong to the same branch that the objects or
#' small clusters being assigned belong to.
#' @param ... any other parameter compatible with
#' \code{\link[WGCNA]{mergeCloseModules}}
#'
#' @return list containing modules detected, modules_eigengenes, and if asked
#' for, modules pre-merge and dendrogram
#'
#' @importFrom WGCNA mergeCloseModules
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats setNames as.dist hclust
#' @importFrom methods is
#' @importFrom stringr str_sort
#'
#' @examples
#' df <- kuehne_expr[1:24, 1:350]
#' net <- build_net(df, n_threads = 1)
#' detect_modules(df, net$network)
#'
#' @export

detect_modules <- function(data_expr, network, min_module_size =
                             min(20, ncol(data_expr) / 2),
                           clustering_th = NULL, merge_close_modules = TRUE,
                           merge_threshold = 0.75, detailled_result = TRUE,
                           pam_respects_dendro = FALSE, ...) {
  # Checks
  if (methods::is(data_expr, "SummarizedExperiment")) {
    data_expr <- t(SummarizedExperiment::assay(data_expr))
  } else .check_data_expr(data_expr)
  if (!(is.data.frame(network) | is.matrix(network)))
    stop("network should be a data.frame or a matrix.")
  if (ncol(network) != nrow(network))
    stop("network should be squarred")
  if (is.null(rownames(network)) | !all(colnames(network) %in%
                                        rownames(network)))
    stop("network should have the same genes names as colnames and rownames")
  # TODO : finish checks

  # Order network matrix in the same order as data_expr
  network <- network[colnames(data_expr), colnames(data_expr)]

  # Hierarchical clustering
  gene_tree = stats::hclust(stats::as.dist(network), method = "average")
  # Tree cut
  dynamicMods = quiet(dynamicTreeCut::cutreeDynamic(
    dendro = gene_tree, distM = network, cutHeight = clustering_th,
    deepSplit = 2, pamRespectsDendro = pam_respects_dendro,
    minClusterSize = min_module_size))
  # TODO manage to divide ellipsis (...) to allow users to add args
  # for cutreeDynamic

  # Re-assign gene names
  dynamicMods <- stats::setNames(dynamicMods, colnames(data_expr))

  if (length(table(dynamicMods)) == 1)
    stop("No modules detected")
  if (length(table(dynamicMods)) == 2)
    warning("Only one module found, plus 'non-classified genes' module")

  # Merging closest modules
  if (merge_close_modules) {
    if (length(table(dynamicMods)) == 1)
      stop("All genes were classified in the same module")
    merge <- quiet(WGCNA::mergeCloseModules(
      data_expr, colors = dynamicMods, cutHeight = 1 - merge_threshold,
      relabel = TRUE, ...))

    # Reordering module eigengenes
    merge$newMEs <- merge$newMEs[, stringr::str_sort(colnames(merge$newMEs),
                                                     numeric = TRUE)]
    merge$newMEs <- merge$newMEs[, stringr::str_sort(colnames(merge$newMEs),
                                                     numeric = TRUE)]

  } else {
    merge <- list(
      colors = dynamicMods,
      newMEs = WGCNA::moduleEigengenes(data_expr, dynamicMods)$eigengenes)
  }

  # Re-formating
  modules_list <- stats::setNames(merge$colors,
                           colnames(data_expr)) %>% split(names(.), .)
  modules_list_premerge <- stats::setNames(dynamicMods,
                                    colnames(data_expr)) %>% split(names(.), .)

  # Return
  if (detailled_result) {
    detection <- list(
      modules = modules_list,
      modules_premerge = modules_list_premerge,
      modules_eigengenes = merge$newMEs,
      dendrograms = stats::hclust(stats::as.dist(1 - cor(merge$newMEs)),
                                  method = "average")
    )
  } else {
    detection <- list(
      modules = modules_list,
      modules_eigengenes = merge$newMEs
    )
  }

  return(detection)
}


# Removing errors about dplyr data-variables
utils::globalVariables(c("before"))

#' Modules merge plot
#'
#' Plot a bipartite graph to see in which modules all modules have been merged
#'
#' @param modules_premerge vector, id (whole number or string) of module before
#' merge associated to each gene.
#' @param modules_merged vector, id (whole number or string) of module after
#' merge associated to each gene.
#'
#' @details Both vectors must be in the same gene order before passing them to
#' the function. No check is applied on this.
#'
#' @return A bipartite graph
#'
#' @importFrom igraph graph_from_data_frame V add_layout_ as_bipartite
#' @importFrom magrittr %>% set_colnames
#' @importFrom dplyr left_join mutate_if distinct arrange
#' @importFrom utils stack
#'
#' @examples
#' df <- kuehne_expr[1:24, 1:350]
#' net <- build_net(df, n_threads = 1)
#' detection <- detect_modules(df, net$network, detailled_result = TRUE)
#' plot_modules_merge(modules_premerge = detection$modules_premerge,
#'                    modules_merged = detection$modules)
#'
#' @export

plot_modules_merge <- function(modules_premerge, modules_merged) {
  # Checks
  if (!is.list(modules_premerge))
    stop("modules_premerge must be a named list with modules id as names, and",
    " vectors of gene names as content")
  if (!is.list(modules_merged))
    stop("modules_merged must be a named list with modules id as names, and",
    " vectors of gene names as content")
  if (length(modules_premerge) < 2)
    stop("modules_premerge must contain at least 2 modules")
  if (length(modules_merged) < 1 |
      length(modules_merged) > length(modules_premerge))
    stop("modules_premerge must contain at least 1 module and less or equal",
    "number of module")
  if (is.null(names(modules_premerge)))
    stop("modules_premerge must be named with modules id as names ")
  if (is.null(names(modules_merged)))
    stop("modules_merged must be named with modules id as names ")
  if (!all(lapply(modules_merged, is.vector, "character") %>% unlist))
    stop("module_merged values for each module must be a vector of gene names")

  # data.frame indicating which module got merge into another
  g <- dplyr::left_join(utils::stack(modules_premerge),
                        utils::stack(modules_merged), by = "values") %>%
    tibble::column_to_rownames("values") %>%
    magrittr::set_colnames(c("before", "after")) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate_if(is.character, as.numeric) %>%
    dplyr::distinct() %>%
    dplyr::arrange(before) %>%
    igraph::graph_from_data_frame()

  # Adding property used by bipartite plot
  igraph::V(g)$type <- lapply(igraph::V(g) %>% names %>% as.numeric,
                              function(x) ifelse(x %in% (modules_merged %>%
                                                           names %>%
                                                           unique),
                                                 TRUE, FALSE)) %>% unlist

  # Bipartite plot
  g %>% igraph::add_layout_(igraph::as_bipartite()) %>%
    plot(vertex.label.color = "gray20",
         vertex.label.family = "Helvetica",
         vertex.label.cex = 0.9,
         vertex.color = "lightskyblue",
         vertex.frame.color = "white")
}


# Removing errors about dplyr data-variables
utils::globalVariables(c("module", "gene", "gene_gene", "expression_gene",
                         "expression_eigengene", "cor_sign"))

#' Modules expression profiles
#'
#' Plot expression profiles for all modules with eigengene highlighted
#'
#' @param data_expr matrix or data.frame or SummarizedExperiment, expression
#' data with genes as column and samples as row.
#' @param modules vector, id (whole number or string) of modules associated to
#' each gene.
#'
#' @return A ggplot representing expression profile and eigengene by module
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_grid theme element_text
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate left_join group_by rename select
#' @importFrom stats setNames prcomp
#' @importFrom utils stack
#' @importFrom methods is
#'
#' @examples
#' df <- kuehne_expr[1:24, 1:350]
#' net <- build_net(df, n_threads = 1)
#' detection <- detect_modules(df, net$network, detailled_result = TRUE)
#' plot_expression_profiles(df, detection$modules)
#'
#' @export

plot_expression_profiles <- function(data_expr, modules) {
  # Check
  if (methods::is(data_expr, "SummarizedExperiment")) {
    data_expr <- t(SummarizedExperiment::assay(data_expr))
  } else .check_data_expr(data_expr)
  if (is.list(modules)) {
    if (any(!unlist(lapply(modules, is.vector, "character")))) {
      stop("If modules is a list of modules, all elements of the list must be",
      " vectors of gene names") }
    if (is.null(names(modules)))
      warning("No name provided for the list of modules.")
  } else if (!is.vector(modules, "character")) {
    stop("modules must be either a list of modules or a single module",
    " represented by a vector of gene names")
  } else if (!is.list(modules) & length(modules) < 2) {
    warning("modules represent a single modules and only contains one",
    " gene name")
  }

  # Tables preparation for ggplot
  if (is.list(modules)) {
    df <- modules %>%
      utils::stack() %>%
      stats::setNames(c("gene", "module")) %>%
      dplyr::mutate(module = as.character(module))
  } else {
    df <- data.frame(gene = modules,
                     module = "module",
                     stringsAsFactors = FALSE) %>%
      dplyr::mutate(module = as.character(module))
  }

  eigengenes <- df %>%
    dplyr::left_join(data_expr %>%
                       t %>%
                       as.data.frame %>%
                       tibble::rownames_to_column(var = "gene"),
                     by = "gene") %>%
    split.data.frame(.$module) %>%
    lapply(function(y) y[,c(-1,-2)] %>%
             t %>%
             stats::prcomp(center = TRUE, scale.=TRUE) %>%
             .$x %>%
             .[, "PC1"] %>%
             scale) %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column(var="sample") %>%
    dplyr::mutate(gene = "eigengene") %>%
    tidyr::pivot_longer(c(-gene, -sample),
                        names_to = "module",
                        values_to = "expression")

  plot_table <- df %>%
    dplyr::left_join(data_expr %>%
                       scale %>%
                       as.data.frame %>%
                       t %>%
                       as.data.frame %>%
                tibble::rownames_to_column(var = "gene"),
                by = "gene") %>%
    tidyr::pivot_longer(c(-gene, -module),
                        names_to = "sample",
                        values_to = "expression") %>%
    dplyr::left_join(eigengenes,
                     c("module", "sample"),
                     suffix = c("_gene", "_eigengene")) %>%
    dplyr::group_by(gene_gene) %>%
    dplyr::mutate(cor_sign = ifelse(
      cor(expression_gene, expression_eigengene) >= 0, "+", "-")) %>%
    dplyr::mutate(expression_eigengene = ifelse(
      cor_sign == "-", -(expression_eigengene), expression_eigengene)) %>%
    dplyr::rename(gene = gene_gene) %>%
    dplyr::rename(expression = expression_gene) %>%
    dplyr::select(gene, module, sample, expression, cor_sign,
                  expression_eigengene)

  ggplot2::ggplot(plot_table,
                  ggplot2::aes(x=sample, y=expression, group=gene)) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::facet_grid(cor_sign ~ module) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_line(ggplot2::aes(x = sample, y = expression_eigengene),
                       color = "red")
}
