#' Compare modules topology between conditions
#'
#' Take modules built from multiples conditions and search for preservation or non-preservation of them against one or mutliple
#' conditions of reference. Use 7 topological features to perform the differents test, and use permutation to validate results.
#'
#' @param data_expr_list list of data_expr, list of expression data with genes as column and samples as row.
#' @param net_list list of networks, list of square tables representing connectivity between
#' each genes as returned by build_net.
#' @param cor_list list of matrices and/or data.frames, list of square tables representing correlation between
#' each gene. Must be the same used to create networks in \code{\link{net_list}}. If NULL, will be re-calculated
#' according to \code{cor_func}.
#' @param ref string, condition name to be used as reference for permutation tests. Default will be "cross_comparison", meaning
#' each condition will be tested as a reference.
#' @param modules_list list of modules, list by condition of modules detected in each one. If only reference modules are provided,
#' no contingency matrix will be computed.
#' @param cor_func string, name of the correlation function to be used. Must be one of "pearson", "spearman",
#' "bicor", "other". If "other", your_func must be provided
#' @param your_func function returning correlation values. Final values must be in [-1;1]
#'
#' @details
#' \describe{
#'   \item{To avoid recalculation, correlations matrices can be obtain by setting \code{keep_cor_mat} in \code{\link{net_list}} to TRUE.}
#'   \item{Description of the 7 topological features used for preservation testing is available in \code{\link[NetRep]{modulePreservation}}}
#' }
#'
#' @importFrom NetRep modulePreservation
#' @importFrom magrittr %>%
#'
#' @export

compare_modules = function(data_expr_list, net_list, cor_list = NULL, ref = "cross_comparison", modules_list,
                           cor_func = c("pearson", "spearman", "bicor", "other"), your_func = NULL,
                           n_perm = 10000, comparison_modules = c("unpreserved", "preserved"), pvalue_th,
                           n_threads, ...){
  # Checks
  if (!is.list(data_expr_list)) stop("data_expr_list must be a list.")
  if (!is.list(net_list)) stop("net_list must be a list.")
  if (!is.list(cor_list) && !is.null(cor_list)) stop("cor_list must be a list or NULL")
  if (!all(names(data_expr_list) %in% names(net_list))) stop("Names between data_expr_list and net_list don't correspond.")
  if (!all(names(data_expr_list) %in% names(cor_list))) stop("Names between data_expr_list and cor_list don't correspond.")
  lapply(data_expr_list, .check_data_expr)
  lapply(net_list, .check_network)
  if (!is.character(ref)) stop("ref must be a string")
  if (length(ref > 1)) stop("ref must be a single string")
  if (ref != "cross_comparison" && !(ref %in% names(data_expr_list))) stop("ref must be a condition name, one of the names of data_expr_list elements.")
  if (!is.list(modules_list)) stop("modules_list must be a single list of modules or a list by condition of the modules")
  if (all(names(modules_list) %in% names(data_expr_list))) { # Meaning it's a list by condition of modules
    lapply(modules_list, function(cond){ .check_module(cond, is_list = TRUE) })
    is_only_ref_modules = FALSE
  } else { # Meaning it must be a single list of modules
    .check_module(modules_list, is_list = TRUE)
    is_only_ref_modules = TRUE
  }
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" && (is.null(your_func) || !is.function(your_func))) stop("If you specify other, your_func must be a function.")
  if (!is.numeric(n_perm)) stop("n_perm must be a numeric value")
  if (n_perm %% 1 != 0) stop("n_perm must be a whole number")
  comparison_modules <- match.arg(comparison_modules)
  if (!is.numeric(n_threads)) stop("n_threads must be a numeric value")
  if (n_threads %% 1 != 0) stop("n_threads must be a whole number")
  if (!is.numeric(n_threads)) stop("n_threads must be a numeric value")
  if (!is.numeric(pvalue_th)) stop("pvalue_th must be a numeric value")
  if (length(pvalue_th) != 1) stop("pvalue_th must be a single value")
  if (pvalue_th <= 0 || pvalue_th >= 1) stop("pvalue_th must be in ]0;1[")

  # If no cor_list provided, calculating one
  if (is.null(cor_list)) {
    cor_list <- sapply(data_expr_list, function(data_expr){
      if (cor_func == "other") {
        cor_to_use <- your_func
      } else {
        cor_to_use <- .cor_func_match(cor_func)
      }
      cor_mat <- cor_to_use(data_expr %>% as.matrix)
      if (min(cor_mat < -1) || max(cor_mat) > 1) stop("Provided correlation function returned values outside [-1,1].")
    }, simplify = FALSE)
  }

  # Params preparation
  conditions <- names(data_expr_list)

  if (is_only_ref_modules) {
    # ref must be a single condition (because modulePreservation require to have the same conditions passed to moduleAssignements and discovery)
    if (length(ref) > 1) {
      warning("Modules given for multiple condition, keeping only the one matching ref provided")
      modules_list <- modules_list[[ref]]
    }
    modules_reformated <- lapply(names(modules_list), function(x){
      setNames(rep(as.numeric(x), length(modules_list[[x]])), modules_list[[x]]) }) %>% unlist
  } else {
    # modules_reformated must be formated like single condition but for each condition, with the same list names
    modules_reformated <- sapply(modules_list, function(cond){
      lapply(names(cond), function(x){
        setNames(rep(as.numeric(x), length(cond[[x]])), cond[[x]]) }) %>% unlist
    }, simplify = FALSE)
  }

  test_tail_side <- ifelse(comparison_modules != "unpreserved", "greater", "lesser")
  if (ref == "cross_comparison") {
    test_set <- conditions
    ref_set <- conditions
    # Must check that discovery names are the same as moduleAssignments
    # if () stop("")
  } else {
    test_set <- conditions[which(conditions != ref)]
    ref_set <- ref
  }

  # Preservation
  preservation <- quiet(NetRep::modulePreservation(
    network = net_list, data = data_list, correlation = cor_list,
    moduleAssignments = modules_reformated, discovery = ref_set, test = test_set,
    alternative = test_tail_side, nPerm = n_perm, nThreads = n_threads,
    simplify = FALSE, ...))


  # modulePreservation doesn't compute contingency matrix when using only one discovery set, so adding it
  if (ref != "cross_comparison") {
    lapply(test_set, function(cond){
      preservation[[ref]][[cond]][["contingency"]] <- NetRep:::contingencyTable(modules_reformated,
                                                                                modules_reformated[[ref]] %>% table %>% names,
                                                                                modules_reformated[[cond]] %>% names)
    })

  }

  # Formating depending on number of conditions tested
  if (length(test_set > 1)) {

  } else {

  }


  # Modules preserved (those which p value < 0.01)
  max_pval <- apply(preservation$p.value, 1, max)
  preserved_modules = which(max_pval < pvalue_th)

  ## Looking for statistics involved into non concervation of the modules
  unpreserved_modules = which(max_pval > pvalue_th)
  statsImplied = apply(as.data.frame(unpreserved_modules), 1, function(x){
    res = names(x)[which(x>0.01)]
  })




}
