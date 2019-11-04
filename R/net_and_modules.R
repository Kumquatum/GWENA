#' Match a correlation function based on a name
#'
#' Translate a function name into a usable function in code.
#' TODO finish
#'
#' @param cor_func string of the name of the correlation to be use
#'
#' @return A function corresponding to the correlation required
#' @examples
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
#' @param network_type TODO
#' TODO finish params
#'
#' @return A list containing power of the law for best fit, fit table, and metadata about the arguments used.
#' @examples
#' # TODO Add example
#'
#' @importFrom WGCNA pickSoftThreshold.fromSimilarity
#'
#' @export

get_fit.cor <- function(cor_mat, fit_cut_off = 0.90, network_type = c("unsigned", "signed", "signed hybrid", "none"), ...){
  # Checking args
  # TODO Add other checks
  if (length(fit_cut_off) != 1 | !is.numeric(fit_cut_off)) stop("power_cut_off should be a single number") # TODO : A revoir
  if (fit_cut_off < 0 | fit_cut_off > 1) stop("power_cut_off should be a number between 0 and 1")
  network_type <- match.arg(network_type)
  if (any(is.na(cor_mat))) stop("cor_mat must not contain any missing value. To approximate them, see FAQ answer on this subject.")
  if (any(cor_mat > 1) || any(cor_mat < -1)) stop("cor_mat should be filled with value in the [-1,1] range")
  if (nrow(cor_mat) != ncol(cor_mat)) stop("cor_mat should be a squared matrix")

  # Calculating similairty
  similarity <- if (network_type == "unsigned") {
    abs(cor_mat)
  } else if (network_type == "signed") {
    (1 + cor_mat) / 2
  } else if (network_type == "signed hybrid") {
    cor_mat[cor_mat < 0] = 0
  } # else = "none" : do nothing

  # Getting fit
  sft_fit <- WGCNA::pickSoftThreshold.fromSimilarity(similarity = similarity, ...)

  # Final list with all infos
  fit <- list(
    power_value = sft_fit$powerEstimate,
    fit_table = sft_fit$fitIndice,
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
#' @param fit_cut_off integer by which R^2 (coefficient of determination) will be thresholded
#' @param cor_func TODO
#' @param your_func TODO
#' @param network_type TODO
#' @param ...
#' TODO filling info on params
#'
#' @return A list containing power of the law for best fit, fit table, and metadata about the arguments used.
#'
#' @examples
#' #TODO write example
#'
#' @importFrom magrittr %>%
#'
#' @export

get_fit.expr <- function(data_expr, fit_cut_off = 0.90, cor_func = c("pearson", "spearman", "bicor", "other"),
                     your_func = NULL, network_type = c("unsigned", "signed", "signed hybrid", "none"), ...){
  # Checking args
  # TODO ajouter autres checks
  if (any(is.na(data_expr))) stop("data_expr must not contain any missing value. To approximate them, see FAQ answer on this subject.")
  if (ncol(data_expr) < nrow(data_expr)) warning("Number of columns inferior to number of rows. Check if columns are the genes name")
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" && (is.null(your_func) || !is.function(your_func))) stop("If you specify other, your_func must be a function.")


  # Cor selection
  if (cor_func == "other") {
    cor_to_use <- your_func
  } else {
    cor_to_use <- cor_func_match(cor_func)
  }

  # Calculating correlation matrix
  cor_mat <- cor_to_use(data_expr %>% as.matrix)

  # Getting fit
  fit <- get_fit.cor(cor_mat = cor_mat, fit_cut_off = fit_cut_off, network_type = network_type, ...)

  # Final list with all infos
  fit$metadata$func_used = ifelse(cor_func == "other",
                                  paste("other:", deparse(substitute(your_func))),
                                  cor_func)

  return(fit)
}



#' Network build and modules detection
#'
#' Compute the adjacency matrix, then the TOM to build the network. Than detect the modules by hierarchical clustering and thresholding
#'
#' @param data_expr matrix of expression data with genes as column and samples as row.
#' @param power_value soft-thresolding, meaning power which will be applied to the adjacency matrix.
#' @param power_cut_off if no power provided, computes power and thresold it by this value.
#' @param cor_func string specifying correlation function to be used. Must be one of pearson, spearman, bicor.
#' TODO finish
#'
#' @details
#' \enumerate{
#'   \item{data_expr}{Must have been adequatly normalized and filtered if pertinant. ATTENTION : it is not recommended to filter by differential
#'   expression (cf. \href{https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html}{Peter Langfelder and Steve Horvath FAQ on WGCNA}}.
#'   \item{save_adjacency}{Saving adjacency increase the final return size. Working with tom is usually suffisant since it's the only value used for the
#' next step \code{\link{modules_detection}}.}
#' }
#'
#' @return list containing
#' @examples
#'
#' @importFrom WGCNA adjacency.fromSimilarity TOMsimilarity
#' @importFrom magrittr %>% set_colnames set_rownames
#'
#' @export

net_building <- function(data_expr, cor_func = c("pearson", "spearman", "bicor", "other"), your_func = NULL,
                         power_value = NULL, fit_cut_off = 0.90, network_type = c("unsigned", "signed", "signed hybrid", "none"),
                         tom_type = c("unsigned", "signed", "signed Nowick"), save_adjacency = FALSE, n_threads = 0, # TODO program the mclapply version
                         detailled_result = FALSE, ...)
{
  # Checking
  # TODO add check
  cor_func <- match.arg(cor_func)
  network_type <- match.arg(network_type)
  tom_type <- match.arg(tom_type)

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
  } else {fit <- "None. Custom power_value provided"}

  ### Network building
  message("Adjacency")
  adj = WGCNA::adjacency.fromSimilarity(similarity = cor_mat, type = network_type, power = fit$power_value)
  message("TOM")
  tom = 1 - WGCNA::TOMsimilarity(adj, TOMType = tom_type) %>%
    magrittr::set_colnames(colnames(adj)) %>%
    magrittr::set_rownames(rownames(adj))


  # Return
  if (detailled_result) {
    net = list(
      power = fit$power_value,
      fit_power = fit$fit_table[fit$power_value, "SFT.R.sq"],
      cor_func = cor_func,
      network_type = network_type,
      adjacency = if (save_adjacency) adj else NULL,
      tom = tom
    )
  } else {
    net = list(
      cor_func = cor_func,
      network_type = network_type,
      tom_type = tom_type,
      tom = tom
    )
  }

  return(net)
}



#' Modules detection from a network matrix
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
