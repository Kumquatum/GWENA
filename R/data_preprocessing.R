#' Filtering genes with low variability
#'
#' Remove low variating genes based on the percentage given and the type of variation specified.
#'
#' @param data_expr matrix or data.frame or SummarizedExperiment, table of expression values (either microarray or RNA-seq),
#' with genes as column and samples as row
#' @param pct float, percentage of gene to keep, value must be in ]0;1[
#' @param type string, function name used for filtration. Should be either "mean", "median", or "mad"
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom dplyr top_frac select
#' @importFrom SummarizedExperiment assay
#'
#' @return A data.frame of filtered genes
#'
#' @examples
#' df <- matrix(abs(rnorm(15*45)), 15)
#' colnames(df) <- paste0("gene_", 1:ncol(df))
#' rownames(df) <- paste0("sample_", 1:nrow(df))
#' df_filtered <- filter_low_var(df)
#'
#' @export

filter_low_var <- function(data_expr, pct = 0.8, type = c("mean", "median", "mad")){
  # Checking args
  if (is(data_expr, "SummarizedExperiment")) {
    data_expr <- t(SummarizedExperiment::assay(data_expr))
  } else .check_data_expr(data_expr)
  if (!is.numeric(pct) | length(pct) != 1) stop("pct should be a single number")
  if (pct <= 0 | pct >= 1) stop("pct should be between 0 and 1")
  type <- match.arg(type)

  # Convert to data.frame if matrix
  if (!is.data.frame(data_expr)) data_expr <- data_expr %>% as.data.frame

  # Calculating variation
  var <- lapply(data_expr, function(row) do.call(type, list(row)))

  # Filtering
  top_pct <- data.frame(gene = names(var), var = unlist(var), stringsAsFactors = FALSE) %>%
    dplyr::top_frac(pct, var)
  filtered_data_expr <- data_expr %>% dplyr::select(dplyr::one_of(top_pct$gene))

  return(filtered_data_expr)
}


#' Filtering of low counts
#'
#' Keeping genes with at least one sample with count above min_count in RNA-seq data.
#'
#' @param data_expr matrix or data.frame or SummarizedExperiment, table of expression values (either microarray or RNA-seq),
#' with genes as column and samples as row.
#' @param min_count integer, minimal number of count to be considered in method.
#' @param method string, name of the method for filtering. Must be one of "at least one", "mean", or " all"
#'
#' @details Low counts in RNA-seq can bring noise to gene co-expression module building, so filtering them
#' help to improve quality.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr one_of
#' @importFrom matrixStats colMaxs colMeans2 colMins
#' @importFrom SummarizedExperiment assay
#'
#' @return A data.frame of filtered genes
#'
#' @examples
#' df <- matrix(abs(rnorm(15*45)), 15) * 3
#' colnames(df) <- paste0("gene_", 1:ncol(df))
#' rownames(df) <- paste0("sample_", 1:nrow(df))
#' df_filtered <- filter_RNA_seq(df)
#'
#' @export

filter_RNA_seq <- function(data_expr, min_count = 5, method = c("at least one", "mean", "all")){
  # Checking args
  if (is(data_expr, "SummarizedExperiment")) {
    data_expr <- t(SummarizedExperiment::assay(data_expr))
  } else .check_data_expr(data_expr)
  if (!is.numeric(min_count) | length(min_count) != 1) stop("min_count should be a single number")
  if (min_count <= 1) stop("min_count should be superior to 1")
  method <- match.arg(method)

  # Casting to matrix (needed for matrixStats)
  if (!is.matrix(data_expr)) data_expr <- data_expr %>% as.matrix

  # Filtering
  if (method == "at least one") { i <- matrixStats::colMaxs(data_expr) > min_count }
  else if (method == "mean") { i <- matrixStats::colMeans2(data_expr) > min_count }
  else if (method == "all") { i <- matrixStats::colMins(data_expr) > min_count }
  else { stop(paste0("Invalid method value: ", method)) } # Should never be triggered because of check

  # Casting to data.frame
  if (!is.data.frame(data_expr)) data_expr <- data_expr %>% as.data.frame

  return(data_expr[, i])
}
