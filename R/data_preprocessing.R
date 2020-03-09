#' Filtering on low variating genes
#'
#' Remove low variating genes based on the percentage given and the type of variation specified.
#'
#' @param data_expr matrix of data (either microarray or RNA-seq) to be filtered, with genes as column and samples as row
#' @param pct number in 0 and 1 specifying the percentage of gene to keep
#' @param type string which should be either "mean" or "median"
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom dplyr top_frac select
#'
#' @export

filter_low_var <- function(data_expr, pct = 0.8, type = c("mean", "median", "mad")){
  # Checking args
  if (!(is.data.frame(data_expr) || is.matrix(data_expr))) stop("data_expr should be a data.frame or a matrix")
  if (ncol(data_expr) < nrow(data_expr)) warning("Number of columns inferior to number of rows. Check if columns are the genes name")
  if (!is.numeric(pct) || length(pct) != 1) stop("pct should be a single number")
  if (pct <= 0 || pct >= 1) stop("pct should be between 0 and 1")
  type <- match.arg(type)

  # Convert to data.frame if matrix
  if (!is.data.frame(data_expr)) data_expr <- data_expr %>% as.data.frame

  # Calculating variation
  var <- lapply(data_expr, function(row) do.call(type, list(row)))

  # Filtering
  top_pct <- data.frame(gene = names(var), var = unlist(var), stringsAsFactors = FALSE) %>% dplyr::top_frac(pct, var)
  filtered_data_expr <- data_expr %>% dplyr::select(dplyr::one_of(top_pct$gene))

  return(filtered_data_expr)
}


#' Filtering of low counts
#'
#' Keeping genes with at least one sample with count above min_count in RNA-seq data.
#'
#' @param data_expr matrix of RNA-seq data to be filtered, with genes as column and samples as row.
#' @param min_count minimal number of count to be considered in method
#' @param method name of the method for filtering. Must be one of "at least one", "mean", or " all", above min_count.
#'
#' @details Low counts in RNA-seq can bring noise to gene co-expression module building, so filtering them help to improve quality.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#'
#' @export

filter_RNA_seq <- function(data_expr, min_count = 5, method = c("at least one", "mean", "all")){
  # Checking args
  if (!(is.data.frame(data_expr) || is.matrix(data_expr))) stop("data_expr should be a data.frame or a matrix")
  if (ncol(data_expr) < nrow(data_expr)) warning("Number of columns inferior to number of rows. Check if columns are the genes name")
  if (!is.numeric(min_count) || length(min_count) != 1) stop("min_count should be a single number")
  if (min_count <= 1) stop("min_count should be superior to 1")
  method <- match.arg(method)

  # Convert to data.frame if matrix
  if (!is.data.frame(data_expr)) data_expr <- data_expr %>% as.data.frame

  # Filtering
  good_gene <- lapply(data_expr, function(x) {
    {
    if (method == "at least one") {any(x > min_count)}
    else if (method == "mean") {mean(x) > min_count}
    else {all(x > min_count)} # Meaning "all"
    }
    }) %>% unlist %>% .[which(. == TRUE)] %>% names

  filtered_data_expr <- data_expr %>% dplyr::select(one_of(good_gene))

  return(filtered_data_expr)
}

