#' Match a correlation function based on a name
#'
#' Translate a function name into an R function.
#'
#' @param cor_func string of the name of the correlation to be use
#'
#' @return A function corresponding to the correlation required
#' @examples
#' cor_func_match("pearson")
#'
#' @importFrom WGCNA bicor cor

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

#' Calculating best fit of a power low on correlation matrix computed on expression data
#'
#' Adjust a correlation matrix depending of the type of network, then try to parameter a power law for best fit
#'
#' @param cor_mat matrix or data.frame of genes correlation.
#' @param fit_cut_off integer by which R^2 (coefficient of determination) will be thresholded.
#' @param network_type string giving type of network to be used. Either "unsigned", "signed", "signed hybrid". See details.
#' @param block_size integer giving size of blocks by which operations can be proceed. Helping if working with low capacity computers. If null, will be estimated.
#'
#' @details
#' network_type indicate which transformation will be applied on the correlation matrix to return the similarity score.
#' \describe{
#'   \item{signed}{will modify the range [-1;1] to [0.5;1.5] (because of log10 beeing used for scale free index computation)}
#'   \item{unsigned}{will return absolute value (moving from [-1;1] to [0;1])}
#'   \item{signed hybrid}{will replace all negative values by 0 (moving from [-1;1] to [0;1])}
#' }
#'
#' @return A list containing power of the law for best fit, fit table, and metadata about the arguments used.
#' @examples
#' get_fit.cor(cor_mat = cor(kuehne_expr[, 1:100]))
#'
#' @importFrom WGCNA pickSoftThreshold.fromSimilarity
#' @importFrom magrittr %>%
#' @importFrom dplyr select top_n
#'
#' @export

get_fit.cor <- function(cor_mat, fit_cut_off = 0.90, network_type = c("unsigned", "signed", "signed hybrid"), block_size = NULL, ...){
  # Checking args
  if (!is.matrix(cor_mat)) stop("cor_mat should be a matrix")
  if (any(is.na(cor_mat))) warning("cor_mat should not contain any missing value. It may be due to non variating probe/transcript when you computed your correlation matrix.")
  if ((any(cor_mat > 1) || any(cor_mat < -1)) && !any(is.na(cor_mat))) stop("cor_mat should be filled with value in the [-1,1] range")
  if (nrow(cor_mat) != ncol(cor_mat)) stop("cor_mat should be a squared matrix")
  if (length(fit_cut_off) != 1 | !is.numeric(fit_cut_off)) stop("power_cut_off should be a single number")
  if (fit_cut_off < 0 | fit_cut_off > 1) stop("power_cut_off should be a number between 0 and 1")
  network_type <- match.arg(network_type)
  if (!is.null(block_size) && (block_size < 2 || block_size %% 1 != 0)) stop("If not NULL, block_size must be a whole number > 1")


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
  sft_fit <- quiet(WGCNA::pickSoftThreshold.fromSimilarity(similarity = similarity, RsquaredCut = fit_cut_off, blockSize = block_size, ...))
  fit_above_cut_off <- TRUE
  if (is.na(sft_fit$powerEstimate)) { # If no fit, taking maximum fitting power
    warning("No fitting power could be found for provided fit_cut_off. Taking power for maximum fit. See FAQ for known causes.")
    sft_fit$powerEstimate <- sft_fit$fitIndice %>% dplyr::top_n(1, SFT.R.sq) %>% dplyr::select(Power) %>% as.numeric()
    fit_above_cut_off <- FALSE
  }

  # Final list with all infos
  fit <- list(
    power_value = sft_fit$powerEstimate,
    fit_table = sft_fit$fitIndice,
    fit_above_cut_off = fit_above_cut_off,
    metadata = list(
      net_type = network_type
    )
  )

  return(fit)
}



#' Calculating best fit of a power low on expression data
#'
#' Computes correlation matrix of the gene expression data, adjust it depending of the type of network, then try to
#' parameter a power law for best fit
#'
#' @param data_expr matrix of data only normalized for constructor specificities, with genes as column and samples as row.
#' @param fit_cut_off integer by which R^2 (coefficient of determination) will be thresholded.
#' @param cor_func string specifying correlation function to be used. Must be one of "pearson", "spearman", "bicor", "other". If "other", your_func must be provided
#' @param your_func function returning correlation values. Final values must be in [-1;1]
#' @param network_type string giving type of network to be used. Either "unsigned", "signed", "signed hybrid". See details.
#' @param block_size integer giving size of blocks by which operations can be proceed. Helping if working with low capacity computers. If null, will be estimated.
#'
#' @details
#' network_type indicate which transformation will be applied on the correlation matrix to return the similarity score.
#' \describe{
#'   \item{signed}{will modify the range [-1;1] to [0.5;1.5] (because of log10 beeing used for scale free index computation)}
#'   \item{unsigned}{will return absolute value (moving from [-1;1] to [0;1])}
#'   \item{signed hybrid}{will replace all negative values by 0 (moving from [-1;1] to [0;1])}
#' }
#'
#' @return A list containing power of the law for best fit, fit table, and metadata about the arguments used.
#'
#' @examples
#' get_fit.expr(data_expr = kuehne_expr[, 1:100])
#'
#' @importFrom magrittr %>%
#'
#' @export

get_fit.expr <- function(data_expr, fit_cut_off = 0.90, cor_func = c("pearson", "spearman", "bicor", "other"),
                     your_func = NULL, network_type = c("unsigned", "signed", "signed hybrid"), block_size = NULL, ...){
  # Checking args
  if (!(is.data.frame(data_expr) || is.matrix(data_expr))) stop("data_expr should be a data.frame or a matrix.")
  if (any(is.na(data_expr))) stop("data_expr cannot contain any missing value. To approximate them, see FAQ answer on this subject.")
  if (min(data_expr) < 0) stop("data_expr cannot contain any negative value.")
  if (ncol(data_expr) < nrow(data_expr)) warning("Number of columns inferior to number of rows. Check if columns are the genes name.")
  # No need to check in theory because checked in get_fit.cor
  # if (length(fit_cut_off) != 1 | !is.numeric(fit_cut_off)) stop("power_cut_off should be a single number.")
  # if (fit_cut_off < 0 | fit_cut_off > 1) stop("power_cut_off should be a number between 0 and 1.")
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" && (is.null(your_func) || !is.function(your_func))) stop("If you specify other, your_func must be a function.")
  # network_type <- match.arg(network_type)

  # Cor selection
  if (cor_func == "other") {
    cor_to_use <- your_func
  } else {
    cor_to_use <- cor_func_match(cor_func)
  }

  # Calculating correlation matrix
  cor_mat <- cor_to_use(data_expr %>% as.matrix)

  # If personnal function, checking it returns a matrix with values in [-1;1]
  if (cor_func == "other"){
    if (min(cor_mat) < -1 || max(cor_mat) > 1) stop("your_func must be a function which returns values in [-1;1]")
  }

  # Getting fit
  fit <- get_fit.cor(cor_mat = cor_mat, fit_cut_off = fit_cut_off, network_type = network_type, block_size = block_size, ...)

  # Final list with all infos
  fit$metadata$cor_func = ifelse(cor_func == "other",
                                 paste("other:", deparse(substitute(your_func))),
                                 cor_func)

  return(fit)
}



#' Network build and modules detection
#'
#' Compute the adjacency matrix, then the TOM to build the network. Than detect the modules by hierarchical clustering and thresholding
#'
#' @param data_expr matrix or data.frame, expression data with genes as column and samples as row.
#' @param fit_cut_off integer, cut off by which R^2 (coefficient of determination) will be thresholded.
#' @param cor_func string, name of the correlation function to be used. Must be one of "pearson", "spearman", "bicor", "other". If "other", your_func must be provided
#' @param your_func function returning correlation values. Final values must be in [-1;1]
#' @param power_value integer, power to be applied to the adjacency matrix. If NULL, will be estimated by trying different power law fitting.
#' @param block_size integer, size of blocks by which operations can be proceed. Helping if working with low capacity computers. If null, will be estimated.
#' @param stop_if_no_fit boolean, does not finding a fit above fit_cut_off should stop process, or just print a warning and return the highest fitting power.
#' @param network_type string, type of network to be used. Either "unsigned", "signed", "signed hybrid". See details.
#' @param tom_type string, type of the topological overlap matrix to be computed. Either "none", "unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2"
#' and "signed Nowick 2". See detail at \code{\link[WGCNA]{TOMsimilarityFromExpr}}.
#' @param save_adjacency string, folder's path where adjacency matrix will be saved. If NULL, it is not saved.
#' @param n_threads integer, number of threads that can be used to paralellise the computing
#'
#' @return list containing network matrix, metadata of input parameters and power fitting information.
#' @examples
#' net_building(kuehne_expr[, 1:1000])
#'
#' @importFrom WGCNA adjacency.fromSimilarity TOMsimilarity
#' @importFrom magrittr %>% set_colnames set_rownames
#'
#' @export

net_building <- function(data_expr, fit_cut_off = 0.90, cor_func = c("pearson", "spearman", "bicor", "other"), your_func = NULL,
                         power_value = NULL, block_size = NULL, stop_if_no_fit = FALSE, network_type = c("unsigned", "signed", "signed hybrid"),
                         tom_type = c("unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2", "none"), save_adjacency = NULL,
                         n_threads = 0, ...)  # TODO program the mclapply version
{
  # Checking
  if (!(is.data.frame(data_expr) || is.matrix(data_expr))) stop("data_expr should be a data.frame or a matrix.")
  if (any(is.na(data_expr))) stop("data_expr cannot contain any missing value. To approximate them, see FAQ answer on this subject.")
  if (min(data_expr) < 0) stop("data_expr cannot contain any negative value.")
  if (ncol(data_expr) < nrow(data_expr)) warning("Number of columns inferior to number of rows. Check if columns are the genes name.")
  # No need to check in theory because checked in get_fit.cor
  # if (length(fit_cut_off) != 1 | !is.numeric(fit_cut_off)) stop("power_cut_off should be a single number.")
  # if (fit_cut_off < 0 | fit_cut_off > 1) stop("power_cut_off should be a number between 0 and 1.")
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" && (is.null(your_func) || !is.function(your_func))) stop("If you specify other, your_func must be a function.")
  if (!is.null(power_value) && (power_value < 1 || power_value %% 1 != 0)) stop("If not NULL, power_value must be a whole number >= 1.")
  network_type <- match.arg(network_type)
  tom_type <- match.arg(tom_type)
  if (!is.null(save_adjacency) && grepl(".+\\.\\w+$", save_adjacency)) warning("Provided path in save_adjacency looks like a filename. Remember save_adjacency must be a folder name.")
  if (!is.null(save_adjacency) && !file.exists(save_adjacency)) stop("Provided path in save_adjacency doesn't exists.")
  if (!is.numeric(n_threads) || n_threads < 0 || n_threads %% 1 != 0) stop("n_threads must be a whole number >= 0")

  # Correlation selection and correlation matrix computation
  if (cor_func == "other") {
    cor_to_use <- your_func
  } else {
    cor_to_use <- cor_func_match(cor_func)
  }
  cor_mat <- cor_to_use(data_expr %>% as.matrix)

  # Getting power
  if (is.null(power_value)) {
    fit <- get_fit.cor(cor_mat = cor_mat, fit_cut_off = fit_cut_off, network_type = network_type, ...)
    if (stop_if_no_fit && fit$fit_above_cut_off == FALSE) stop("No fitting power could be found for provided fit_cut_off. You should verify your data (or lower fit_cut_off). See FAQ.")
  } else {
    fit <- list(
      power_value = power_value,
      fit_table = "None. Custom power_value provided")
    }

  # Adjacency
  adj = WGCNA::adjacency.fromSimilarity(similarity = cor_mat, type = network_type, power = fit$power_value)
  if (!is.null(save_adjacency)) {
    write.csv(adj, file = file.path(save_adjacency, "adjacency_matrix.csv"))
  }

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

  return(net)
}



#' Modules detection in a network
#'
#' Detect the modules by hierarchical clustering .
#'
#' @param data_expr matrix of expression data with genes as column and samples as row.
#' @param tom
#' TODO finish
#'
#' @details
#' TODO
#' @return
#' TODO
#'
#' @examples
#' #TODO
#' @importFrom WGCNA mergeCloseModules
#' @importFrom dynamicTreeCut cutreeDynamic
#'
#' @export

modules_detection <- function(data_expr, tom, min_module_size = min(20, ncol(data_expr) / 2), merge_cut_height = 0.25,
                              detailled_result = FALSE, ...) {
  # Checks


  # Detection
  message("Hierarchical clustering")
  gene_tree = stats::hclust(as.dist(tom), method = "average")
  message("Tree cut")
  dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = gene_tree, distM = tom,
                                              deepSplit = 2, pamRespectsDendro = FALSE,
                                              minClusterSize = min_module_size)
  message("Merging closest modules")
  merge = WGCNA::mergeCloseModules(data_expr, dynamicMods, cutHeight = merge_cut_height, relabel = TRUE, ...)

  # Return
  if (detailled_result) {
    detection <- list(
      modules = setNames(merge$colors, colnames(data_expr)),
      modules_premerge = dynamicMods,
      modules_eigengenes = merge$newMEs,
      dendrograms = stats::hclust(as.dist(1 - cor(merge$newMEs)), method = "average")
    )
  } else {
    detection <- list(
      modules = setNames(merge$colors, colnames(data_expr)),
      modules_eigengenes = merge$newMEs
    )
  }

  return(detection)
}



#' Merging modules' plot
#'
#' Plot a bipartite graph to see in which modules all modules have been merged
#'
#' @param
#' TODO finish
#'
#' @details
#' TODO
#' @return
#' TODO
#'
#' @examples
#' #TODO
#' @importFrom igraph graph_from_data_frame V add_layout_ as_bipartite
#' @importFrom magrittr %>%
#'
#' @export

plot_modules_merge <- function(modules) {
  g <- data.frame(before = modules$modules_premerge , after = modules$modules, stringsAsFactors = FALSE) %>%
    distinct %>%
    arrange(before) %>%
    igraph::graph_from_data_frame()
  igraph::V(g)$type <- lapply(igraph::V(g) %>% names %>% as.numeric, function(x) ifelse(x %in% unique(modules$modules), TRUE, FALSE)) %>% unlist
  g %>% igraph::add_layout_(igraph::as_bipartite()) %>%
    plot(vertex.label.color = "gray20",
         vertex.label.family = "Helvetica",
         vertex.label.cex = 0.9,
         vertex.color = "lightskyblue",
         vertex.frame.color = "white")
}


#' Plot all modules expression profiles
#'
#' f
#'
#' @param
#' TODO finish
#'
#' @details
#' TODO
#' @return
#' TODO
#'
#' @examples
#' #TODO
#' @importFrom ggplot2 ggplot aes geom_line facet_grid theme element_text
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#'
#' @export

plot_expression_profiles <- function(data_expr, modules) {
  df <- data.frame(gene = names(modules$modules), module = modules$modules, stringsAsFactors = FALSE)

  eigengenes <- df %>%
    left_join(data_expr %>% t %>% as.data.frame %>%
                tibble::rownames_to_column(var = "gene"), by = "gene") %>%
    split.data.frame(.$module) %>%
    lapply(function(y) y[,c(-1,-2)] %>% t %>% prcomp(center = TRUE, scale.=TRUE) %>% .$x %>% .[, "PC1"] %>% scale) %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column(var="sample") %>%
    mutate(gene = "eigengene") %>%
    tidyr::pivot_longer(c(-gene, -sample), names_to = "module", values_to = "expression", names_ptypes = list(module = numeric()))

  plot_table <- df %>%
    left_join(data_expr %>% scale %>% as.data.frame %>% t %>% as.data.frame %>%
                tibble::rownames_to_column(var = "gene"), by = "gene") %>%
    tidyr::pivot_longer(c(-gene, -module), names_to = "sample", values_to = "expression") %>%
    left_join(eigengenes, c("module", "sample"), suffix = c("_gene", "_eigengene")) %>%
    group_by(gene_gene) %>%
    mutate(cor_sign = ifelse(cor(expression_gene, expression_eigengene) >= 0, "+", "-")) %>%
    mutate(expression_eigengene = ifelse(cor_sign == "-", -(expression_eigengene), expression_eigengene)) %>%
    rename(gene = gene_gene) %>%
    rename(expression = expression_gene) %>%
    select(gene, module, sample, expression, cor_sign, expression_eigengene)

  ggplot2::ggplot(plot_table, ggplot2::aes(x=sample, y=expression, group=gene)) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::facet_grid(cor_sign ~ module) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90,vjust = 0.5)) + ggplot2::theme(legend.position = "none") +
    ggplot2::geom_line(ggplot2::aes(x = sample, y = expression_eigengene), color = "red")
}
