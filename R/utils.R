# ==== Format ====

#' Muting a function
#'
#' Prevent a function to output multiple message.
#' Source:
#' https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
#'
#' @param func Function who need to be muted.
#'
#' @return Nothing, just mute the called function

quiet <- function(func) {
  sink(tempfile())
  on.exit(sink())
  suppressMessages(invisible(force(func)))
}

#' Mimicking ggplot palette
#' Source : https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' @param n integer, number of colors wanted
#' 
#' @importFrom grDevices hcl
#' 
#' @return character vector, haxadecimal colors of length n
gg_palette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}


# ==== Checks ====

#' Determine if an object is a network
#'
#' Check content of a given object to determine if it's a network, meaning a
#' squared matrix of similarity score between genes.
#'
#' @param network matrix or data.frame, object to test to be a network
#'
#' @return list, a boolean as first element and in second element NULL or the
#' reason why boolean is set to FALSE
#'
#' @examples
#' net <- matrix(runif(40*40), 40)
#' colnames(net) <- paste0("gene_", seq_len(ncol(net)))
#' rownames(net) <- paste0("gene_", seq_len(nrow(net)))
#' is_network(net)
#'
#' @export

is_network <- function(network) {
  if (!(is.data.frame(network) | is.matrix(network)))
    return(list(bool = FALSE,
                reason = "network must be a data.frame or a matrix"))
  if (is.null(colnames(network)) | is.null(rownames(network)))
    return(list(bool = FALSE,
                reason = "network must have colnames and rownames"))
  if (!all(colnames(network) %in% rownames(network)))
    return(list(bool = FALSE,
                reason = "colnames and rownames form network doesn't match"))
  if (ncol(network) != nrow(network))
    return(list(bool = FALSE,
                reason = "network must be a squared matrix"))
  if ((any(network > 1) | any(network < -1)) & !any(is.na(network)))
    return(list(bool = FALSE,
                reason = paste("network values should be between [-1,1]",
                "range")))

  return(list(bool = TRUE, reason = NULL))
}

#' Run checks on an object to test if it's a network
#'
#' Check content of a given object to determine if it's a network, meaning a
#' squared matrix of similarity
#' score between genes.
#'
#' @param network matrix or data.frame, object to test to be a network
#'
#' @return Throw an error if doesn't correspond

.check_network <- function(network) {
  check <- is_network(network)
  if (!check$bool) {
    stop(check$reason)
  }
}


#' Determine if an object is a module or a list of modules
#'
#' Check content of a given object to determine if it's a module or a list of
#' modules, meaning a single vector of characters which are gene names, or a
#' named list of these vectors
#'
#' @param module vector or list, object to test to be a module or list of
#' modules
#' @param is_list boolean, indicate if module must be tested as a single module
#' or a list of modules
#'
#' @return list, a boolean as first element and in second element NULL or the
#' reason why boolean is set to FALSE
#'
#' @examples
#' single_module <- c("BIRC3", "PMAIP1", "CASP8", "JUN", "BCL2L11", "MCL1",
#'                    "IL1B", "SPTAN1", "DIABLO", "BAX", "BIK", "IL1A", "BID",
#'                    "CDKN1A", "GADD45A")
#' is_module(single_module)
#'
#' multi_module <- list(mod1 = single_module,
#'                      mod2 = c("TAF1C", "TARBP2", "POLH", "CETN2", "POLD1",
#'                      "CANT1", "PDE4B", "DGCR8", "RAD51", "SURF1", "PNP",
#'                      "ADA", "NME3", "GTF3C5", "NT5C"))
#' is_module(multi_module$modules, is_list = TRUE)
#'
#' @export

is_module <- function(module, is_list = FALSE) {
  if (is_list) {
    if (!is.list(module))
      return(list(bool = FALSE,
                  reason = "module must be a list"))
    if (is.null(names(module)))
      return(list(bool = FALSE,
                  reason = "module list must have names"))
    if (any(names(module) == "" %>% unlist))
      return(list(bool = FALSE,
                  reason = "module list must have names for all elements"))
    if (any(lapply(module, is.vector, "character") %>% unlist %>% `!`))
      return(list(bool = FALSE,
                  reason = "module list element must be vector of gene names"))
  } else {
    if (!is.vector(module, "character"))
      return(list(bool = FALSE,
                  reason = "module must be vector of gene names"))
  }
  return(list(bool = TRUE, reason = NULL))
}


#' Run checks on an object to test if it's a module or a list of modules
#'
#' Check content of a given object to determine if it's a module or a list of
#' modules, meaning a single vector of characters which are gene names, or a
#' named list of these vectors
#'
#' @param module vector or list, object to test to be a module or list of
#' modules
#' @param is_list boolean, indicate if module must be tested as a single
#' module or a list of modules
#'
#' @return Throw an error if doesn't correspond

.check_module <- function(module, is_list = FALSE) {
  check <- is_module(module, is_list)
  if (!check$bool) {
    stop(check$reason)
  }
}


#' Determine if an object is a gost object
#'
#' Check content of a given object to determine if it's a gost object
#'
#' @param gost_result list, gprofiler2::gost result
#'
#' @return list, a boolean as first element and in second element NULL or the
#' reason why boolean is set to FALSE
#'
#' @examples
#' single_module <- c("BIRC3", "PMAIP1", "CASP8", "JUN", "BCL2L11", "MCL1",
#'                    "IL1B", "SPTAN1", "DIABLO", "BAX", "BIK", "IL1A", "BID",
#'                    "CDKN1A", "GADD45A")
#' single_module_enriched <- bio_enrich(single_module)
#' is_gost(single_module_enriched)
#'
#' @export

is_gost <- function(gost_result) {
  if (!is.list(gost_result))
    return(list(bool = FALSE,
                reason = "gost_result must be a list."))
  if (is.null(gost_result))
    return(list(bool = FALSE,
                reason = "Elements of gost_result cannot be NULL"))
  if (!all(names(gost_result) %in% c("result", "meta")))
    return(list(bool = FALSE,
                reason = paste("gprofiler2::gost first levels should be",
                "'result' and 'meta'")))
  if (!is.data.frame(gost_result$result))
    return(list(bool = FALSE,
                reason = "'result' should be a data.frame"))
  if (any(is.na(match(c("query", "significant", "p_value", "term_size",
                        "query_size", "intersection_size", "precision",
                        "recall", "term_id", "source", "term_name",
                        "effective_domain_size", "source_order", "parents"),
                      colnames(gost_result$result)))))
    return(list(bool = FALSE,
                reason = "'result' is not a gprofiler2::gost result output"))
  if (!is.list(gost_result$meta))
    return(list(bool = FALSE,
                reason = "meta should be a list"))
  if (any(is.na(match(c("query_metadata", "result_metadata", "genes_metadata",
                        "timestamp", "version"),
                      names(gost_result$meta)))))
    return(list(bool = FALSE,
                reason = paste("Bad format: 'meta' is not a gprofiler2::gost",
                " result output")))
  return(list(bool = TRUE,
              reason = NULL))
}


#' Run checks on an object to test if it's a gost result
#'
#' Take a list that should be a gost result and check if format is good.
#'
#' @param gost_result list, gprofiler2::gost result
#'
#' @return Throw an error if doesn't correspond

.check_gost <- function(gost_result) {
  check <- is_gost(gost_result)
  if (!check$bool) {
    stop(check$reason)
  }
}


#' Determine if an object is a data_expr in sens of GWENA
#'
#' Check an object to be a data.frame or a matrix compatible of genes and
#' samples.
#'
#' @param data_expr matrix or data.frame, expression data with genes as column
#' and samples as row.
#'
#' @return list, a boolean as first element and in second element NULL or the
#' reason why boolean is set to FALSE
#'
#' @examples
#' expr <- matrix(runif(15*40), 15)
#' colnames(expr) <- paste0("gene_", seq_len(ncol(expr)))
#' rownames(expr) <- paste0("gene_", seq_len(nrow(expr)))
#' is_data_expr(expr)
#'
#' @export

is_data_expr <- function(data_expr) {
  if (!(is.data.frame(data_expr) | is.matrix(data_expr)))
    return(list(bool = FALSE,
                reason = "data_expr should be a data.frame or a matrix."))
  if (any(is.na(data_expr)))
    return(list(bool = FALSE,
                reason = paste("data_expr cannot contain any missing value.",
                "To approximate them, see FAQ answer on this subject.")))
  if (!(is.numeric(unlist(data_expr))))
    return(list(bool = FALSE),
           reason = "data_expr must contain only numeric values.")
  # HF_remove_limit_neg_value_build_net : commenting this check because
  # data_expr can normalized in negative values (like centering)
  # if (min(data_expr) < 0) return(list(bool = FALSE, reason = "data_expr
  # cannot contain any negative value."))
  if (ncol(data_expr) < nrow(data_expr))
    warning("Number of columns inferior to number of rows. Check if columns",
    " are the genes name.")
  if (is.null(colnames(data_expr)) | is.null(rownames(data_expr)))
    return(list(bool = FALSE,
                reason = "data_expr should have colnames and rownames"))
  return(list(bool = TRUE,
              reason = NULL))
}


#' Run checks on an object to test if it's a data_expr
#'
#' Check an object to be a data.frame or a matrix compatible of genes and
#' samples.
#'
#' @param data_expr matrix or data.frame, expression data with genes as
#' column and samples as row.
#'
#' @return Throw an error if doesn't correspond

.check_data_expr <- function(data_expr){
  check <- is_data_expr(data_expr)
  if (!check$bool) {
    stop(check$reason)
  }
}
