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
#' @param ref string, condition name to be used as reference for permutation tests or "cross comparison" if you want to compare
#' each condition with the other as reference. Default will be the name of the first element
#' in data_expr_list.
#' @param modules_list list of modules or nested list of modules, list of modules in one condition (will be considered as the one from reference) or a condition
#' named list with list of modules built in each one.
#' @param cor_func string, name of the correlation function to be used. Must be one of "pearson", "spearman",
#' "bicor", "other". If "other", your_func must be provided
#' @param your_func function returning correlation values. Final values must be in [-1;1]
#'
#' @details
#' \describe{
#'   \item{To avoid recalculation, correlations matrices can be obtain by setting \code{keep_cor_mat} in \code{\link{net_list}} to TRUE.}
#'   \item{Description of the 7 topological features used for preservation testing is available in \code{\link[NetRep]{modulePreservation}}.}
#' }
#'
#' @importFrom NetRep modulePreservation
#' @importFrom magrittr %>%
#'
#' @export

compare_modules = function(data_expr_list, net_list, cor_list = NULL, modules_list, ref = names(data_expr_list)[1],
                           cor_func = c("pearson", "spearman", "bicor", "other"), your_func = NULL,
                           n_perm = 10000, comparison_modules = c("unpreserved", "preserved", "one or the other"), pvalue_th = 0.05,
                           n_threads = NULL, ...){
  # Checks
  if (!is.list(data_expr_list)) stop("data_expr_list must be a list.")
  if (length(data_expr_list) < 2) stop("data_expr_list must have at least 2 elements (2 conditions) to run a comparison")
  if (!is.list(net_list)) stop("net_list must be a list.")
  if (!is.list(cor_list) && !is.null(cor_list)) stop("cor_list must be a list or NULL")
  if (!all(names(data_expr_list) %in% names(net_list))) stop("Names between data_expr_list and net_list don't correspond.")
  if (!all(names(data_expr_list) %in% names(cor_list))) stop("Names between data_expr_list and cor_list don't correspond.")
  lapply(data_expr_list, .check_data_expr)
  lapply(net_list, .check_network)
  if (!is.character(ref)) stop("ref must be a string or a vector of strings")
  if (ref != "cross comparison" && !(all(ref %in% names(data_expr_list)))) stop("ref must be a condition name, or a vector of it, matching names of data_expr_list.")
  if (!is.list(modules_list)) stop("modules_list must be a single list of modules or a list by condition of the modules")
  if (all(names(modules_list) %in% names(data_expr_list))) { # Meaning it's a list of modules by condition
    lapply(modules_list, function(cond){ .check_module(cond, is_list = TRUE) })
    is_modules_by_cond = TRUE
  } else { # Meaning it must be a single list of modules
    .check_module(modules_list, is_list = TRUE)
    is_modules_by_cond = FALSE
  }
  if (ref != "cross comparison") { # Checking adequation of content between conditions described in modules and reference choosen
    if (!all(ref %in% names(modules_list)) && is_modules_by_cond == TRUE) stop("All conditions cited in ref must have their modules described in modules_list")
    if (!all(names(modules_list) %in% ref) && is_modules_by_cond == TRUE) {
      # warning("More conditions were defined in modules_list than those cited in ref. Keeping only the matching ones")
      additionnal_modules_list <- modules_list[which(names(modules_list) != ref)]
      modules_list <- modules_list[ref]
    }
  } else {
      if (!all(names(data_expr_list) %in% names(modules_list))) stop("All conditions in data_expr_list must have their modules described in modules_list")
  }
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" && (is.null(your_func) || !is.function(your_func))) stop("If you specify other, your_func must be a function.")
  if (!is.numeric(n_perm)) stop("n_perm must be a numeric value")
  if (n_perm %% 1 != 0) stop("n_perm must be a whole number")
  comparison_modules <- match.arg(comparison_modules)
  if (!is.null(n_threads)) {
    if (!is.numeric(n_threads)) stop("n_threads must be a numeric value")
    if (n_threads %% 1 != 0) stop("n_threads must be a whole number")
    if (!is.numeric(n_threads)) stop("n_threads must be a numeric value")
  }
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
  if (is_modules_by_cond) {
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
  test_tail_side <- switch(comparison_modules,
                           "unpreserved" = "less",
                           "preserved" = "greater",
                           "one or the other" = "two.sided"
                           )
  if (ref == "cross comparison") {
    test_set <- conditions
    ref_set <- conditions
  } else {
    test_set <- conditions[which(!(conditions %in% ref))]
    ref_set <- ref
  }

  # Preservation
  preservation <- quiet(NetRep::modulePreservation(
    network = net_list, data = data_list, correlation = cor_list,
    moduleAssignments = modules_reformated, discovery = ref_set, test = test_set,
    alternative = test_tail_side, nPerm = n_perm, nThreads = n_threads,
    simplify = FALSE, ...))


  #' modulePreservation doesn't compute contingency matrix when single ref because it force to also have a single moduleAssignement
  #' which prevent to compute contingency matrix (and more generaly, it force to have exactly the same conditions in discovery and
  #' moduleAssignment. But I want to allow people to pass moduleAssignement for all cond tested even if they're not in discovery,
  #' so adding the functionnality
  if (ref != "cross comparison" && exists("additionnal_modules_list")) {
    lapply(test_set, function(cond){
      preservation[[ref]][[cond]][["contingency"]] <- NetRep:::contingencyTable(modules_reformated,
                                                                                modules_reformated[[ref]] %>% table %>% names,
                                                                                modules_reformated[[cond]] %>% names)
    })
  }

  preservation_augmented <- lapply(preservation, function(ref_set){
    lapply(ref_set, function(test_set){
      modules_of_interest <- preservation[[ref_set]][[test_set]][["p.values"]] %>%
        apply(1, function(mod_stats){
          if (all(mod_stats < pvalue_th/2)) { "unpreserved"
          } else if (all(mod_stats > (1 - pvalue_th/2))) { "preserved"
          } else { "" }
        }) %>% data.frame(module = names(.), comparison = .)

      list(modules_of_interest = modules_of_interest,
           detail = test_set)
    })
  })

  # Simplifying final object (yes, modulePreservation has one but I need to get the none simplified version to add contingency matrix correctly)
  prune_list <- function(x) {
    x <- lapply(x, function(y) if (is.list(y)) prune_list(y) else y)
    purrr::compact(x)
  }
  preservation_augmented <- prune_list(preservation_augmented)
}


