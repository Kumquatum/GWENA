#' Compare modules topology between conditions
#'
#' Take modules built from multiples conditions and search for preservation, non-preservation or one of them, against one or mutliple
#' conditions of reference. Use 7 topological features to perform the differents test, and use permutation to validate results.
#'
#' @param data_expr_list list of data_expr, list of expression data by condition, with genes as column and samples as row.
#' @param net_list list of networks, list of square tables by condition, representing connectivity between
#' each genes as returned by build_net.
#' @param cor_list list of matrices and/or data.frames, list of square tables by condition, representing correlation between
#' each gene. Must be the same used to create networks in \code{\link{build_net}}. If NULL, will be re-calculated
#' according to \code{cor_func}.
#' @param modules_list list of modules or nested list of modules, list of modules in one condition (will be considered as the one
#' from reference) or a condition named list with list of modules built in each one.
#' @param ref string or vector of strings, condition(s) name to be used as reference for permutation tests, or "cross comparison"
#' if you want to compare each condition with the other as reference. Default will be the name of the first element in data_expr_list.
#' @param test string or vector of strings, condition(s) name to be tested for permutation tests. If NULL, all conditions except
#' these in ref will be taken. If ref is set to "cross comparison", any test specified will be ignored.
#' @param cor_func string, name of the correlation function to be used. Must be one of "pearson", "spearman",
#' "bicor", "other". If "other", your_func must be provided
#' @param your_func function returning correlation values. Final values must be in [-1;1]
#' @param n_perm integer, number of permutation, meaning number of random gene name re-assignment inside network to compute all tests
#' and statistics for module comparison between condition.
#' @param comparison_type string, either "unpreserved", "preserved" or "one or the other". Design if the comparison aim to detect
#' preserved modules between condition, unpreserved ones, or modules that are one or the other without specification of which.
#' @param pvalue_th decimal, threshold of pvalue below which comparison_type is considered significant. If "one or the other", then
#' pvalue_th is splitted in two for each side (preserved/not preserved).
#' @param n_threads integer, number of threads that can be used to paralellise the computing
#' @param ... any other parameter compatible with \code{\link[NetRep]{modulePreservation}}
#'
#' @details
#' \describe{
#'   \item{}{Conditions will be based on names of data_expr_list. Please do not use numbers for conditions names as modules are
#'   often named this way}
#'   \item{}{To avoid recalculation, correlations matrices can be obtain by setting \code{keep_cor_mat} in \code{\link[GWENA]{build_net}} to TRUE.}
#'   \item{}{Description of the 7 topological features used for preservation testing is available in \code{\link[NetRep]{modulePreservation}}.}
#' }
#'
#' @importFrom NetRep modulePreservation
#' @importFrom magrittr %>%
#' @importFrom purrr compact
#'
#' @export
#'
#' @return A nested list where first element is each ref provided, second level each condition to test, and then elements containing
#' information on the comparison. See

compare_conditions = function(data_expr_list, net_list, cor_list = NULL, modules_list, ref = names(data_expr_list)[1], test = NULL,
                           cor_func = c("pearson", "spearman", "bicor", "other"), your_func = NULL, n_perm = 10000,
                           comparison_type = c("unpreserved", "preserved", "one or the other"), pvalue_th = 0.05, n_threads = NULL, ...){
  # Basic checks
  if (!is.list(data_expr_list)) stop("data_expr_list must be a list.")
  if (length(data_expr_list) < 2) stop("data_expr_list must have at least 2 elements (2 conditions) to run a comparison.")
  lapply(data_expr_list, .check_data_expr)
  conditions <- names(data_expr_list)
  if (!is.list(net_list)) stop("net_list must be a list.")
  if (!all(conditions %in% names(net_list))) stop("Names in net_list don't match with conditions.")
  lapply(net_list, .check_network)
  if (!is.list(cor_list) && !is.null(cor_list)) stop("cor_list must be a list or NULL")
  if (!is.null(cor_list)) {
    if (!all(conditions %in% names(cor_list))) stop("Names in cor_list don't match with conditions.")
    all_cor <- cor_list %>% unlist # avoid to do a lapply
    if (min(all_cor) < -1 || max(all_cor) > 1) stop("Provided correlation_list contains values outside [-1,1].")
  }
  if (!is.character(ref)) stop("ref must be a string or a vector of strings")
  if (all(ref != "cross comparison") && !(all(ref %in% conditions))) stop("ref must be a condition name, or a vector of it, matching conditions.")
  if (length(ref) > 1 && any("cross comparison" %in% ref)) stop("If multiple ref given, they should be condition names, and none 'cross comparison'")
  if (!is.null(test) && !is.character(test)) stop("test must be null, or a string, or a vector of strings")
  if (is.null(test)) { test <- conditions[which(!(conditions %in% ref))]
  } else if (!(all(test %in% conditions))) { stop("test must be a condition name, or a vector of it, matching conditions.") }
  if (!is.list(modules_list)) stop("modules_list must be a single list of modules or a list by condition of the modules")
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" && (is.null(your_func) || !is.function(your_func))) stop("If you specify other, your_func must be a function.")
  if (!is.numeric(n_perm)) stop("n_perm must be a numeric value")
  if (n_perm %% 1 != 0) stop("n_perm must be a whole number")
  comparison_type <- match.arg(comparison_type)
  if (!is.null(n_threads)) {
    if (!is.numeric(n_threads)) stop("n_threads must be a numeric value")
    if (n_threads %% 1 != 0) stop("n_threads must be a whole number")
    if (!is.numeric(n_threads)) stop("n_threads must be a numeric value")
  }
  if (!is.numeric(pvalue_th)) stop("pvalue_th must be a numeric value")
  if (length(pvalue_th) != 1) stop("pvalue_th must be a single value")
  if (pvalue_th <= 0 || pvalue_th >= 1) stop("pvalue_th must be in ]0;1[")

  # Basic checks on modules_list (single cond or multiple cond list ?) + removing conditions into it which are not defined in data_expr_list
  match_cond_modules_data <- names(modules_list) %in% conditions # Looking at conditions cited or not in modules_list
  if (any(match_cond_modules_data)) { # Meaning it could be a list of modules by condition
    if (!all(match_cond_modules_data)) { # Removing conditions defined in modules_list but not present in data_expr_list
      warning("Conditions defined in modules_list not present in data_expr_list. Removing them.")
      modules_list <- modules_list[which(match_cond_modules_data)]
    }
    lapply(modules_list, function(cond){ .check_module(cond, is_list = TRUE) }) # Checking each cond is a list of modules
    is_multi_cond_modules_list = TRUE
  } else { # Meaning it must be a single list of modules
    .check_module(modules_list, is_list = TRUE)
    is_multi_cond_modules_list = FALSE
  }

  # Checks conditions in ref are defined in modules_list and taking aside conditions defined in modules_list which are not in ref (but
  # still in conditions defined into data_expr_list) for later contingency computing.
  if (all(ref != "cross comparison")) { # Checking adequation of content between conditions described in modules and reference choosen
    if (is_multi_cond_modules_list) {
      if (!all(ref %in% names(modules_list))) stop("All conditions cited in ref must have their modules described in modules_list")
      if (!all(names(modules_list) %in% ref)) {
        additional_modules_list <- modules_list[which(!(names(modules_list) %in% ref))]
        modules_list <- modules_list[ref]
      }
    } else {
      if (length(ref) > 1) stop("modules_list must have at least a list of modules for all condition cited in ref")
    }
  } else {
    if (!all(conditions %in% names(modules_list))) {
      stop("All conditions must have their modules described in modules_list when ref is 'cross comparison'")}
  }


  # Process

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
      return(cor_mat)
    }, simplify = FALSE)
  } else { cor_func <- "unknown"}

  # Adapting modules format
  if (is_multi_cond_modules_list) {
    # list by condition of gene named vector of modules values
    modules_reformated <- sapply(modules_list, function(cond){
      lapply(names(cond), function(x){
        setNames(rep(as.numeric(x), length(cond[[x]])), cond[[x]]) }) %>% unlist
    }, simplify = FALSE)
  } else {
    # Single gene named vector of modules values
    modules_reformated <- lapply(names(modules_list), function(x){
      setNames(rep(as.numeric(x), length(modules_list[[x]])), modules_list[[x]]) }) %>% unlist
  }
  test_tail_side <- switch(comparison_type,
                           "unpreserved" = "less",
                           "preserved" = "greater",
                           "one or the other" = "two.sided"
                           )
  if (all(ref == "cross comparison")) {
    test_set <- conditions
    ref_set <- conditions
  } else {
    test_set <- test
    ref_set <- ref
  }

  # Preservation
  preservation <- quiet(NetRep::modulePreservation(
    network = net_list, data = data_expr_list, correlation = cor_list,
    moduleAssignments = modules_reformated, discovery = ref_set, test = test_set,
    alternative = test_tail_side, nPerm = n_perm, nThreads = n_threads,
    simplify = FALSE, ...))


  # modulePreservation doesn't compute contingency matrix when single ref because it force to also have a single moduleAssignement
  # which prevent to compute contingency matrix (and more generaly, it force to have exactly the same conditions in discovery and
  # moduleAssignment. But I want to allow people to pass moduleAssignement for all cond tested even if they're not in discovery,
  # so adding the functionnality
  if (ref != "cross comparison" && exists("additional_modules_list")) {
    additional_modules_reformated <- sapply(additional_modules_list, function(cond){
      lapply(names(cond), function(x){
        setNames(rep(as.numeric(x), length(cond[[x]])), cond[[x]]) }) %>% unlist
    }, simplify = FALSE)
    for (ref_i in ref_set) {
      for (additional_j in names(additional_modules_list)) {
        tmp_module_labels <- list(ref_i = modules_reformated[[ref_i]],
                             additional_j = additional_modules_reformated[[additional_j]])
        preservation[[ref_i]][[additional_j]][["contingency"]] <- contingencyTable(tmp_module_labels,
                                                                                   tmp_module_labels$ref_i %>% table %>% names,
                                                                                   tmp_module_labels$additional_j %>% names)[["contingency"]]
      }
    }
  }

  # Simplifying final object (yes, modulePreservation has a simplification parameter but I need to get the not simplified
  # version to add contingency matrix correctly)
  prune_list <- function(x) {
    x <- lapply(x, function(y) if (is.list(y)) prune_list(y) else y)
    purrr::compact(x)
  }
  preservation <- prune_list(preservation)

  # Adding summarising table about module's preservation
  for (ref_i in names(preservation)) {
    for (test_j in names(preservation[[ref_i]])) {
      comparison <- preservation[[ref_i]][[test_j]][["p.values"]] %>%
        apply(1, function(mod_stats){
          if (all(mod_stats < pvalue_th/2)) { "unpreserved"
          } else if (all(mod_stats > (1 - pvalue_th/2))) { "preserved"
          } else { "not significant" }
        }) %>% data.frame(module = names(.), comparison = ., stringsAsFactors = FALSE)

      preservation[[ref_i]][[test_j]][["comparison"]] <- comparison
    }
  }

  # Adding metadata about the comparison
  preservation <- list(result = preservation,
                       metadata = list(
                         pvalue_th = pvalue_th,
                         cor_func = cor_func,
                         comparison_type = comparison_type
                       ))

  return(preservation)
}
