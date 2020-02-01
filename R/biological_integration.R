#' Join gprofiler2::gost results
#'
#' Takes list of gprofiler2::gost results and join them. Usefull to join results of gprofiler2::gost with custom gmt to
#' other gprofiler2::gost results.
#'
#' @param gost_result list of gprofiler2::gost result
#'
#' @details First element of the list is taken as reference for checks on gost_result elements compatibility. If warnings
#' returned, value from reference will be used.
#' Also, timestamp is set to timestamp of the join
#'
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' @importFrom rlist list.merge
#'
#' @export

join_gost <- function(gost_result) {
  # Check format
  if (!is.list(gost_result)) stop("gost_result must be a list.")
  if (length(gost_result) < 2) stop("gost_result must contain at least 2 gprofiler2::gost element")
  lapply(gost_result, function(gost){
    if (is.null(gost)) stop("Elements of gost_result cannot be NULL")
    if (!all(names(gost) %in% c("result", "meta"))) stop("Bad format: gprofiler2::gost first levels should be 'result' and 'meta'")
    if (!is.data.frame(gost$result)) stop("Bad format: 'result' should be a data.frame")
    if (any(is.na(match(c("query", "significant", "p_value", "term_size", "query_size", "intersection_size", "precision", "recall",
                          "term_id", "source", "term_name", "effective_domain_size", "source_order", "parents"),
                        colnames(gost$result))))) stop("Bad format: 'result' is not a gprofiler2::gost result output")
    if (!is.list(gost$meta)) stop("Bad format: meta should be a list")
    if (any(is.na(match(c("query_metadata", "result_metadata", "genes_metadata", "timestamp", "version"),
                        names(gost$meta))))) stop("Bad format: 'meta' is not a gprofiler2::gost result output")
  })

  # Checking content is compatible (from the same queries) using 'meta' and the first element of the list as reference
  ref <- gost_result[[1]]
  lapply(2:length(gost_result), function(x){
     addon <- gost_result[[x]]
    lapply(names(ref$meta$query_metadata), function(name) {
      if (!(name %in% c("organism", "sources"))) { # They won't be equivalent anyway
        if (name == "queries") {
          if (length(ref$meta$query_metadata$queries) != length(addon$meta$query_metadata$queries)) stop("Number of queries different.")
          if (!isTRUE(all.equal(lapply(ref$meta$query_metadata$queries, length) %>% unlist %>% sort,
                         lapply(addon$meta$query_metadata$queries, length) %>% unlist %>% sort))) stop("Length of queries different.")
          if (length(setdiff(ref$meta$query_metadata$queries, addon$meta$query_metadata$queries)) > 0) warning(
            "Queries different between reference (first element of gost_result) and gost_result element n°", x, ". It may be due to different type of ID (Ensembl, Entrez, etc.).
            IDs from reference will be kept.")
        } else if (name == "numeric_ns") {
          if (is.vector(ref$meta$query_metadata$numeric_ns, "numeric") || is.vector(addon$meta$query_metadata$numeric_ns, "numeric")) {
            if (ref$meta$query_metadata$numeric_ns != addon$meta$query_metadata$numeric_ns) {
              warning("Different type of IDs. ", ref$meta$query_metadata$numeric_ns, " (reference) and ", addon$meta$query_metadata$numeric_ns,
                      " (list element n°", x, ")")}}
        } else {
          # if (ref$meta$query_metadata[[name]] != addon$meta$query_metadata[[name]]) stop(paste0("Item ", name, " of query_metadata isn't identical"))
          if (!identical(ref$meta$query_metadata[[name]], addon$meta$query_metadata[[name]])) warning(
            name, " from query_metadata isn't identical between reference and gost_result element n°", x, ".")
        }
      }
    })
  })

  # Joining
    tryCatch({
      ref$result <- lapply(gost_result, "[[", "result") %>% bind_rows()
      ref$meta$query_metadata$sources <- lapply(gost_result, "[[", "meta") %>% lapply("[[", "query_metadata") %>% lapply("[[", "sources") %>%
        unlist %>% unique
      ref$meta$result_metadata <- lapply(gost_result, "[[", "meta") %>% lapply("[[", "result_metadata") %>% unlist(recursive = FALSE) %>% rlist::list.merge()
      ref$meta$timestamp <- strftime(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", "GMT") # Note : not exactly the same format as the one from API
    }, error = function(err){
      message("Joining failed with the following error: ", err)
      message("Check the content of the gprofiler2::gost results given in list")
    }, warning = function(warn){
      message("Warning while joining: ", warn)
      message("Check the content of the gprofiler2::gost results given in list")
    })

  return(ref)
}


#' Modules enrichment
#'
#' Enrich genes list from modules.
#'
#' @param module vector or list, vector of gene names representing a module or a named list of this modules.
#' @param custom_gmt string or list, path to a gmt file or a list of these path.
#' @param ... any other parameter you can provide to gprofiler2::gost function.
#'
#' @return a gprofiler2::gost output, meaning a named list containing a 'result' data.frame with enrichement information on the
#' differents databases and custom gmt files, and a 'meta' list containing informations on the input args, the version of gost,
#' timestamp, etc. For more detail, see ?gprofiler2::gost.
#'
#' @examples
#' # TODO
#' @importFrom gprofiler2 gost
#' @importFrom plyr compact
#' @importFrom magrittr %>%
#'
#' @export

bio_enrich <- function(module, custom_gmt = NULL, ...) {
  # Checks
  if (is.list(module)) {
    if (any(!unlist(lapply(module, is.vector, "character")))) {
      stop("If module is a list of modules, all elements of the list must be vectors of gene names") }
    if (is.null(names(module))) warning("No name provided for the list of modules.")
  } else if (!is.vector(module, "character")) {
    stop("module must be either a list of modules or a single module represented by a vector of gene names")
  }
  if (!is.null(custom_gmt)) {
    if (!is.character(custom_gmt) && !is.list(custom_gmt)) stop("custom_gmt must be a path or a list of path to gmt file(s)")
    if (is.list(custom_gmt) && !all(lapply(custom_gmt, is.character) %>% unlist)) stop("all element of the list must be paths")
  }

  # Enrichment with gprofiler internal datasets
  enriched_modules <- gprofiler2::gost(query = module, ...)

  # Enrichment with custom gmt if provided
  if (!is.null(custom_gmt)) {
    if (!is.list(custom_gmt)) custom_gmt <- list(custom_gmt)
    list_res_custom_gmts <- lapply(custom_gmt, function(gmt) {
      gmt_id <- quiet(gprofiler2::upload_GMT_file(gmt))
      enriched_modules_gmt <- quiet(gprofiler2::gost(query = module, organism = gmt_id, ...))
      if (is.null(enriched_modules_gmt)) warning("No enrichment found on gmt ", gmt)
      return(enriched_modules_gmt)
    })

    # If no enrichment with custom_gmt, returning only classic gost enrichment
    if (all(lapply(list_res_custom_gmts, is.null) %>% unlist)) {
      stop("None of the custom_gmt file provided returned an enrichement")
      return(enriched_modules)
    }

    # Removing NULL output from the list to merge
    list_res_custom_gmts <- plyr::compact(list_res_custom_gmts)

    # Joining results
    enriched_modules <- join_gost(c(list(enriched_modules), list_res_custom_gmts))
  }

  return(enriched_modules)
}

#' Plot module from bio_enrich
#'
#' Wrapper of the gprofiler2::gostplot function. Adding support of colorblind palet and selection of subsets if
#' initial multiple query, and/or sources to plot.
#'
#' @param enrich_output list, bio_enrich result which are in fact gprofiler2::gost output.
#' @param modules string or vector of characters designing the modules to plot. "all" by default to plot every module.
#' @param sources string or vector of characters designing the sources to plot. "all" by default to plot every source.
#' @param ... any other parameter you can provide to gprofiler2::gostplot
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gprofiler2 gostplot
#'
#' @export

plot_enrichment <- function(enrich_output, modules = "all", sources = "all", colorblind = TRUE, custom_palette = NULL, ...) {
  # Checks
  if (!is.vector(modules, "character")) stop("modules must be a string or vector of characters")
  if (!is.vector(sources, "character")) stop("sources must be a string or vector of characters")
  if (!("all" %in% modules) && !all(modules %in% enrich_output$meta$query_metadata$queries %>% names)) {
    stop("Specified module must be in enrich_output list of modules names") }
  if (!("all" %in% sources) && !all(sources %in% enrich_output$meta$query_metadata$sources)) {
    stop("Specified sources must be in enrich_output list of sources") }
  if (!is.logical(colorblind)) stop("colorblind should be a boolean")
  if (length(enrich_output$meta$query_metadata$sources) > 12) warning("Cannot use a colorblind palette if more than 12 sources")
  if (!is.null(custom_palette) && colorblind) warning("colorblind palette will not be used if custom_palette is provided")

  # Selection of modules subsets if needed
  if (!("all" %in% modules)) {
    enrich_output$result <- enrich_output$result %>% dplyr::filter(query %in% modules)
  }

  # Selection of sources subsets if needed
  if (!("all" %in% sources)) {
    enrich_output$meta$query_metadata$sources <- sources
    enrich_output$result <- enrich_output$result %>% dplyr::filter(source %in% sources)
  }

  # Colorblind pallet
  if (isTRUE(colorblind) && is.null(custom_palette)) {
    palette <- RColorBrewer::brewer.pal(length(enrich_output$meta$query_metadata$sources), "Paired")
  } else palette <- NULL

  # Plotting
  gprofiler2::gostplot(enrich_output, pal = if (is.null(custom_palette)) palette else custom_palette, ...)
}


#' Modules phenotpic association
#'
#' Compute the correlation between all modules and the phenotypic variables
#'
#' @param eigengenes Matrix of igengenes provided by the output of modules_detection.
#' @param phenotypes Matrix of phenotypes for each sample to associate.
#' TODO finish
#'
#' @details
#' TODO
#' @return
#' TODO
#'
#' @examples
#' #TODO
#' @importFrom WGCNA corPvalueStudent
#' @importFrom dplyr select
#'
#' @export

associate_phenotype <- function(eigengenes, phenotypes) {
  # Checks

  # Design matrix (dummy variable formation for qualitative variables)
  dummies_var <- lapply(colnames(phenotypes), function(dummy_name) {
    df <- phenotypes %>% dplyr::select(dummy_name)
    if (!is.numeric(df[1,])) {
      model_mat <- stats::model.matrix(
        stats::formula(paste("~ ", dummy_name, "+ 0")),
        data = df
      )
      colnames(model_mat) <- gsub(dummy_name, "", colnames(model_mat))
      return(model_mat)
    } else {
      return(df)
    }
  })
  design_mat <- as.data.frame(dummies_var, check.names = FALSE)

  # Correlations
  association <- stats::cor(eigengenes, design_mat)

  # P value associated
  pval_association <- WGCNA::corPvalueStudent(association, nrow(phenotypes))

  # TODO Check if a correction for multiple testing shouldn't be performed here...

  return(list(
    association = association %>% as.data.frame,
    pval = pval_association %>% as.data.frame
  ))
}



#' Heatmap of modules phenotpic association
#'
#' Plot a heatmap of the correlation between all modules and the phenotypic variables and the p value associated
#'
#' @param modules_phenotype List form modules_phenotype function output
#' TODO finish
#'
#' @details
#' TODO
#' @return
#' TODO
#'
#' @examples
#' #TODO
#' @importFrom ggplot2 ggplot geom_tile geom_point scale_color_gradient2 theme_bw xlab ylab
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select bind_cols
#'
#' @export

plot_modules_phenotype <- function(modules_phenotype, signif_th = 0.05){
  # Checks


  # Data preparation
  df_cor <- modules_phenotype$association %>%
    tibble::rownames_to_column(var = "eigengene") %>%
    tidyr::pivot_longer(-eigengene, names_to = "phenotype", values_to = "cor")

  df_pval <- modules_phenotype$pval %>%
    tibble::rownames_to_column(var = "eigengene") %>%
    tidyr::pivot_longer(-eigengene, names_to = "phenotype", values_to = "pval")

  df_total <- dplyr::bind_cols(df_cor, df_pval %>% dplyr::select(pval)) %>%
    dplyr::mutate(signif = ifelse(pval > signif_th, FALSE, TRUE))

  # Plotting
  ggplot2::ggplot(df_total, ggplot2::aes(x = factor(eigengene), y = factor(phenotype))) +
    ggplot2::geom_tile(fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = cor, size = signif)) +
    ggplot2::scale_color_gradient2() +
    ggplot2::theme_bw() +
    ggplot2::xlab("Module") +
    ggplot2::ylab("Phenotype")
}
