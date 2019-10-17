#' Filtering on low variating genes
#'
#' Remove low variating genes based on the percentage given and the type of variation specified.
#'
#' @param data_expr matrix of data to be filtered, with genes as column and samples as row
#' @param pct number in 0 and 1 specifying the percentage of gene to keep
#' @param type string which should be either "mean" or "median"
#'
#' @import dplyr
#'
#' @export

filter_low_var <- function(data_expr, pct = 0.8, type = c("mean", "median", "mad")){
  # Checking args
  if (!(is.data.frame(data_expr) || is.matrix(data_expr))) stop("data_expr should be a data.frame or a matrix")
  if (ncol(data_expr) < nrow(data_expr)) warning("Number of columns inferior to number of rows. Check if columns are the genes name")
  if (!is.numeric(pct) || pct <= 0 || pct >= 1) stop("pct should be between 0 and 1")
  type <- match.arg(type)

  # Calculating variation
  var <- lapply(data_expr, function(row) do.call(type, list(row)))

  # Filtering
  top_pct <- data.frame(gene = names(var), var = unlist(var), stringsAsFactors = FALSE) %>% top_frac(pct, var)
  filtered_data_expr <- data_expr %>% select(one_of(top_pct$gene))

  return(filtered_data_expr)
}


#' Filtering of low counts
#'
#' Keeping genes with at least one sample with count above min_count.
#'
#' @param data_expr matrix of RNA-seq data to be filtered, with genes as column and samples as row.
#' @param min_count minimal number of count to be considered in method
#' @param method name of the method for filtering. Must be one of "at least one", "mean", or " all", above min_count.
#'
#' @details Low counts in RNA-seq can bring noise to gene co-expression module building, so filtering them help to improve quality.
#'
#' @import dplyr
#'
#' @export

filter_RNA_seq <- function(data_expr, min_count = 5, method = c("at least one", "mean", "all")){
  # Checking args
  if (!(is.data.frame(data_expr) || is.matrix(data_expr))) stop("data_expr should be a data.frame or a matrix")
  if (ncol(data_expr) < nrow(data_expr)) warning("Number of columns inferior to number of rows. Check if columns are the genes name")
  if (!is.numeric(min_count) || min_count <= 1) stop("min_count should be superior to 1")
  method <- match.arg(method)

  # Filtering
  good_gene <- lapply(data_expr, function(x) {
    {
    if (method == "at least one") {any(x > min_count)}
    else if (method == "mean") {mean(x) > min_count}
    else {all(x > min_count)} # Meaning "all"
    }
    }) %>% unlist %>% .[which(. == TRUE)] %>% names

  filtered_data_expr <- data_expr %>% select(one_of(good_gene))

  return(filtered_data_expr)
}


#' #' Pre-treatment of raw transcriptomic data
#' #'
#' #' Allow specific normalisation and filtration of transcriptomic data, depending of their technology : RNA-seq or microarray
#' #'
#' #' @param data_expr matrix of data only normalized for constructor specificities, with genes as columns and samples as rows
#' #' @param techno string which should be "RNA-seq" or "microarray"
#' #' @param normalisation
#' #' @param filtration
#' #' @param aggregation string designing method to be used for aggregation. Should be either "None", "MaxMean", "MinMean",
#' #'                    "maxRowVariance", "absMaxMean", "absMinMean", "ME", "Average", "Median".
#' #' @param aggreg_group vector of strings of same length as \code{data_expr} id names and in the same order
#' #'
#' #' @return gene expression treated matrix with genes as XXX and samples as XXX
#' #' @author  Gwenaëlle Lemoine <lemoine.gwenaelle@@gmail.com>
#' #' @examples
#' #' exp_mat_path <- system.file("extdata", "h.all.v6.2.entrez.gmt",
#' #'                            package = "gprofiler2.addon", mustWork = TRUE)
#' #'
#' #' @import limma
#' #' @import WGCNA
#' #'
#' #' @export
#'
#' data_preprocessing <- function(data_expr, techno = "RNA-seq", normalisation = "Default", filtration = "Default",
#'                                aggregation = "None", aggreg_group = NULL) {
#'   # Checking args syntax
#'   match.arg(techno, c("RNA-seq", "microarray"))
#'   match.arg(normalisation, c(""))
#'   match.arg(filtration, c(""))
#'   match.arg(aggregation, c("None", "MaxMean", "MinMean", "maxRowVariance", "absMaxMean", "absMinMean", "ME", "Average", "Median"))
#'
#'   # Checking args content
#'   if (!is.vector(aggreg_group) | is.list(aggreg_group) | !is.character(aggreg_group)) stop("aggreg_group should be a vector of strings")
#'   if (aggregation != "None") {
#'     if (is.null(aggreg_group)) stop("aggreg_group should be provided if aggregation method is different of 'None'")
#'     if (length(aggreg_group) != colnames(data_expr)) stop("Length of aggreg_group should be the same as for data_expr id names") # TODO : modifier / verifier si c'est bien les colnames et pas rownames
#'   }
#'
#'   # Process
#'   if (normalisation == FALSE & filtration == FALSE & aggregation == "None") warning("No modification will be done.")
#'   if (techno == "RNA-seq") {
#'     if (filtration == TRUE) {
#'       #
#'     }
#'     if (normalisation == TRUE) {
#'       #
#'     }
#'   } else {# aka microarray
#'     if (filtration == TRUE) {
#'       #
#'     }
#'     if (normalisation == TRUE) {
#'       #
#'     }
#'   }
#'
#'   # Collapsing
#'   if (aggregation) {
#'     WGCNA::collapseRows(data_expr, colnames(data_expr), aggreg_group, method = aggregation) # TODO :
#'   }
#' }
