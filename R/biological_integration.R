#' Join gprofiler2::gost results
#'
#' Takes list of gprofiler2::gost results and join them. Usefull to join
#' results of gprofiler2::gost with custom gmt to other gprofiler2::gost
#' results.
#'
#' @param gost_result list of gprofiler2::gost result
#'
#' @details First element of the list is taken as reference for checks on
#' gost_result elements compatibility. If warnings returned, value from
#' reference will be used. Also, timestamp is set to timestamp of the join
#'
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' @importFrom rlist list.merge
#' @importFrom purrr compact
#'
#' @return A gprofiler2::gost result
#'
#' @examples
#' query <- c("ENSG00000184349", "ENSG00000158955", "ENSG00000091140",
#'            "ENSG00000163114", "ENSG00000163132", "ENSG00000019186")
#' g1 <- gprofiler2::gost(query, sources = "GO")
#' g2 <- gprofiler2::gost(query, sources = "REAC")
#' gj <- join_gost(list(g1,g2))
#'
#' @export

join_gost <- function(gost_result) {
  # Check format
  if (isTRUE(all.equal(names(gost_result), c("result", "meta")))) {
    stop("You provided a single gprofiler2::gost result, it must be a list",
         "of at least 2 gprofiler2::gost")}
  # All elements of gost_result list must be gprofiler2::gost results
  lapply(gost_result, .check_gost)
  if (length(gost_result) < 2)
    stop("gost_result must contain at least 2 gprofiler2::gost element")

  # Checking content is compatible (from the same queries) using 'meta'
  # and the first element of the list as reference
  ref <- gost_result[[1]]
  lapply(2:length(gost_result), function(x){
     addon <- gost_result[[x]]
    lapply(names(ref$meta$query_metadata), function(name) {
      if (!(name %in% c("organism", "sources"))) {# They'll never be equivalent
        if (name == "queries") {
          if (length(ref$meta$query_metadata$queries) !=
              length(addon$meta$query_metadata$queries))
            stop("Number of queries different.")
          if (!isTRUE(all.equal(
            lapply(ref$meta$query_metadata$queries, length) %>%
            unlist %>% sort,
            lapply(addon$meta$query_metadata$queries, length) %>%
            unlist %>% sort)))
            stop("Length of queries different.")
          if (length(setdiff(ref$meta$query_metadata$queries,
                             addon$meta$query_metadata$queries)) > 0)
            warning("Queries different between reference (first element ",
                    "of gost_result) and gost_result element n\u00B0 ", x,
                    ". It may be due to different type of ID (Ensembl, ",
                    "Entrez, etc.). IDs from reference will be kept. ")
        } else if (name == "numeric_ns") {
          if (is.vector(ref$meta$query_metadata$numeric_ns, "numeric") |
              is.vector(addon$meta$query_metadata$numeric_ns, "numeric")) {
            if (ref$meta$query_metadata$numeric_ns !=
                addon$meta$query_metadata$numeric_ns) {
              warning("Different type of IDs. ",
                      ref$meta$query_metadata$numeric_ns,
                      " (reference) and ",
                      addon$meta$query_metadata$numeric_ns,
                      " (list element \u00B0", x, ")")}}
        } else {
          if (!identical(ref$meta$query_metadata[[name]],
                         addon$meta$query_metadata[[name]]))
            warning(name, " from query_metadata isn't identical between ",
            "reference and gost_result element n\u00B0", x, ".")
        }
      }
    })
  })

  # Joining
    tryCatch({
      ref$result <- lapply(gost_result, "[[", "result") %>% bind_rows()
      ref$meta$query_metadata$sources <- lapply(gost_result, "[[", "meta") %>%
        lapply("[[", "query_metadata") %>%
        lapply("[[", "sources") %>%
        unlist %>% unique
      ref$meta$result_metadata <- lapply(gost_result, "[[", "meta") %>%
        lapply("[[", "result_metadata") %>%
        unlist(recursive = FALSE) %>%
        rlist::list.merge()
      ref$meta$timestamp <- strftime(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", "GMT")
      # Note : date will not be exactly in the same format as the one from API
    }, error = function(err){
      message("Joining failed with the following error: ", err)
      message("Check the content of gprofiler2::gost results given in list")
    }, warning = function(warn){
      message("Warning while joining: ", warn)
      message("Check the content of gprofiler2::gost results given in list")
    })

  return(ref)
}


#' Modules enrichment
#'
#' Enrich genes list from modules.
#'
#' @param module vector or list, vector of gene names representing a module or
#' a named list of this modules.
#' @param custom_gmt string or list, path to a gmt file or a list of these
#' path.
#' @param ... any other parameter you can provide to gprofiler2::gost function.
#'
#' @return A gprofiler2::gost output, meaning a named list containing a
#' 'result' data.frame with enrichement information on the
#' differents databases and custom gmt files, and a 'meta' list containing
#' informations on the input args, the version of gost, timestamp, etc.
#' For more detail, see ?gprofiler2::gost.
#'
#' @importFrom gprofiler2 gost
#' @importFrom purrr compact
#' @importFrom magrittr %>%
#'
#' @examples
#' custom_path <- system.file("extdata", "h.all.v6.2.symbols.gmt",
#'                            package = "GWENA", mustWork = TRUE)
#'
#' single_module <- c("BIRC3", "PMAIP1", "CASP8", "JUN", "BCL2L11", "MCL1",
#'                    "IL1B", "SPTAN1", "DIABLO", "BAX", "BIK", "IL1A", "BID",
#'                    "CDKN1A", "GADD45A")
#' single_module_enriched <- bio_enrich(single_module, custom_path)
#'
#' multi_module <- list(mod1 = single_module,
#'                      mod2 = c("TAF1C", "TARBP2", "POLH", "CETN2", "POLD1",
#'                               "CANT1", "PDE4B", "DGCR8", "RAD51", "SURF1",
#'                                "PNP", "ADA", "NME3", "GTF3C5", "NT5C"))
#' multi_module_enriched <- bio_enrich(multi_module, custom_path)
#'
#' @export

bio_enrich <- function(module, custom_gmt = NULL, ...) {
  # Checks
  if (is.list(module)) {
    if (any(!unlist(lapply(module, is.vector, "character")))) {
      stop("If module is a list of modules, all elements of the list must be ",
           "vectors of gene names") }
    if (is.null(names(module)))
      warning("No name provided for the list of modules.")
  } else if (!is.vector(module, "character")) {
    stop("module must be either a list of modules or a single module ",
         "represented by a vector of gene names")
  } else if (!is.list(module) & length(module) < 2) {
    warning("module represent a single module and only contains one gene name")
  }
  if (!is.null(custom_gmt)) {
    if (is.list(custom_gmt)) {
      if (any(lapply(custom_gmt, function(x) !is.character(x) |
                     length(x) != 1 ) %>% unlist)) {
        stop("all element of the list must be string reprensenting paths")}
      if (!all(lapply(custom_gmt, file.exists)))
        stop("all custom_gmt path provided must exist")
    } else if (is.character(custom_gmt) & length(custom_gmt) == 1) {
      if (!file.exists(custom_gmt))
        stop("custom_gmt path provided does not exists")
    } else stop("custom_gmt must be a path or a list of path to gmt file(s)")
  }

  # Enrichment with gprofiler internal datasets
  enriched_modules <- gprofiler2::gost(query = module, ...)

  # Enrichment with custom gmt if provided
  if (!is.null(custom_gmt)) {
    if (!is.list(custom_gmt)) custom_gmt <- list(custom_gmt)
    list_res_custom_gmts <- lapply(custom_gmt, function(gmt) {
      gmt_id <- quiet(gprofiler2::upload_GMT_file(gmt))
      enriched_modules_gmt <- quiet(gprofiler2::gost(query = module,
                                                     organism = gmt_id, ...))
      if (is.null(enriched_modules_gmt))
        warning("No enrichment found on gmt ", gmt)
      return(enriched_modules_gmt)
    })

    # If no enrichment with custom_gmt, returning only classic gost enrichment
    if (all(lapply(list_res_custom_gmts, is.null) %>% unlist)) {
      warning("None of the custom_gmt file provided returned an enrichement")
      return(enriched_modules)
    }

    # Removing NULL output from the list to merge
    # (exist when at least one of the gmt provided return no enrichment for
    # any module)
    list_res_custom_gmts <- purrr::compact(list_res_custom_gmts)

    # Joining results
    enriched_modules <- join_gost(c(list(enriched_modules),
                                    list_res_custom_gmts))
  }

  return(enriched_modules)
}


# Removing errors about dplyr data-variables
utils::globalVariables(c("query"))

#' Plot module from bio_enrich
#'
#' Wrapper of the gprofiler2::gostplot function. Adding support of colorblind
#' palet and selection of subsets if initial multiple query, and/or sources
#' to plot.
#'
#' @param enrich_output list, bio_enrich result which are in fact
#' gprofiler2::gost output.
#' @param modules string or vector of characters designing the modules to plot.
#' "all" by default to plot every module.
#' @param sources string or vector of characters designing the sources to plot.
#' "all" by default to plot every source.
#' @param colorblind boolean, indicates if a colorblind friendly palette should
#' be used.
#' @param custom_palette vector of character, colors to be used for plotting.
#' @param ... any other parameter you can provide to gprofiler2::gostplot.
#'
#' @details Note: The colorblind friendly palette is limited to maximum 8
#' colors, therefore 8 sources of enrichment.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gprofiler2 gostplot
#'
#' @return A plotly object representing enrichment for specified modules
#'
#' @examples
#' custom_path <- system.file("extdata", "h.all.v6.2.symbols.gmt",
#'                            package = "GWENA", mustWork = TRUE)
#' multi_module <- list(mod1 = c("BIRC3", "PMAIP1", "CASP8", "JUN", "BCL2L11",
#'                               "MCL1", "IL1B", "SPTAN1", "DIABLO", "BAX",
#'                               "BIK", "IL1A", "BID", "CDKN1A", "GADD45A"),
#'                      mod2 = c("TAF1C", "TARBP2", "POLH", "CETN2", "POLD1",
#'                               "CANT1", "PDE4B", "DGCR8", "RAD51", "SURF1",
#'                               "PNP", "ADA", "NME3", "GTF3C5", "NT5C"))
#' multi_module_enriched <- bio_enrich(multi_module, custom_path)
#' plot_enrichment(multi_module_enriched)
#'
#' @export

plot_enrichment <- function(enrich_output, modules = "all", sources = "all",
                            colorblind = TRUE, custom_palette = NULL, ...) {
  # Checks
  .check_gost(enrich_output)
  if (!is.vector(modules, "character"))
    stop("modules must be a string or vector of characters")
  if (!is.vector(sources, "character"))
    stop("sources must be a string or vector of characters")
  if (!("all" %in% modules) & !all(
    modules %in% enrich_output$meta$query_metadata$queries %>% names)) {
    stop("Specified module must be in enrich_output list of modules names") }
  if (!("all" %in% sources) & !all(
    sources %in% enrich_output$meta$query_metadata$sources)) {
    stop("Specified sources must be in enrich_output list of sources") }
  if (!is.logical(colorblind)) stop("colorblind should be a boolean")
  if (length(enrich_output$meta$query_metadata$sources) > 12)
    warning("Cannot use a colorblind palette if more than 12 sources")
  if (!is.null(custom_palette) & colorblind)
    warning("colorblind palette will not be used if custom_palette provided")

  # Selection of modules subsets if needed
  if (!("all" %in% modules)) {
    enrich_output$result <- enrich_output$result %>%
      dplyr::filter(query %in% modules)
  }

  # Selection of sources subsets if needed
  if (!("all" %in% sources)) {
    enrich_output$meta$query_metadata$sources <- sources
    enrich_output$result <- enrich_output$result %>%
      dplyr::filter(source %in% sources)
  }

  # Colorblind pallet
  if (isTRUE(colorblind) & is.null(custom_palette)) {
    palette <- RColorBrewer::brewer.pal(
      length(enrich_output$meta$query_metadata$sources), "Paired")
  } else palette <- NULL

  # Plotting
  gprofiler2::gostplot(
    enrich_output,
    pal = ifelse(is.null(custom_palette), palette, custom_palette), ...)
}


#' Modules phenotpic association
#'
#' Compute the correlation between all modules and the phenotypic variables
#'
#' @param eigengenes matrix or data.frame, eigengenes of the modules.
#' Provided by the output of modules_detection.
#' @param phenotypes matrix or data.frame, phenotypes for each sample to
#' associate.
#' @param cor_func string, name of the correlation function to be used. Must be
#' one of "pearson", "spearman", "kendall", "other". If "other", your_func must
#' be provided
#' @param your_func function returning a correlation matrix. Final values must 
#' be in [-1;1] range
#' @param id_col string or vector of string, optional name of the columns
#' containing the common id between eigengenes and phenotypes.
#' @importFrom WGCNA corPvalueStudent
#' @importFrom dplyr select
#'
#' @return A list of two data.frames : associations modules/phenotype and
#' p.values associated to this associations
#'
#' @examples
#' eigengene_mat <- data.frame(mod1 = rnorm(20, 0.1, 0.2),
#' mod2 = rnorm(20, 0.2, 0.2))
#' phenotype_mat <- data.frame(phenA = sample(c("X", "Y", "Z"), 20,
#'                             replace = TRUE),
#'                             phenB = sample(c("U", "V"), 20, replace = TRUE),
#'                             stringsAsFactors = FALSE)
#' association <- associate_phenotype(eigengene_mat, phenotype_mat)
#'
#' @export

associate_phenotype <- function(
  eigengenes, phenotypes, 
  cor_func = c("pearson", "spearman", "kendall", "other"),
  your_func = NULL, id_col = NULL) {
  # Checks
  if (!(is.data.frame(eigengenes) | is.matrix(eigengenes)))
    stop("eigengenes should be a data.frame or matrix")
  if (is.null(colnames(eigengenes)))
    stop("eigengenes should have modules reference as colnames")
  if (!all(lapply(eigengenes, function(x)
    is.character(x) | is.numeric(x) | is.logical(x)) %>% unlist)) {
    stop("eigengenes content should be only characters, numeric, or ",
         "booleans") }
  if (!(is.data.frame(phenotypes) | is.matrix(phenotypes)))
    stop("phenotypes should be a data.frame or matrix")
  if (!all(lapply(phenotypes, function(x)
    is.character(x) | is.numeric(x) | is.logical(x)) %>% unlist)) {
    stop("eigengenes content should be only characters, numeric, or ",
         "booleans") }
  if (nrow(eigengenes) != nrow(phenotypes))
    stop("Number of row should be the same between eigengene and phenotypes ",
         "(samples)")
  cor_func <- match.arg(cor_func)
  if (cor_func == "other" & (is.null(your_func) | !is.function(your_func)))
    stop("If you specify other, your_func must be a function.")
  if (!is.null(id_col)) {
    if (!is.character(id_col))
      stop("id_col should be a character or a character vector")
    if (!any(id_col %in% c(colnames(eigengenes), colnames(phenotypes))))
      stop("Specified id_col wasn't found in colnames")
    if (length(id_col) > 2)
      stop("More than 2 id_col specified")
    if (length(id_col) == 1 & !all(c(id_col %in% colnames(eigengenes),
                                     id_col %in% colnames(phenotypes))))
      stop("Specified id_col wasn't found in both datasets")
    if (length(id_col) == 2 & id_col[1] == id_col[2]) {
      warning("You specified twice the same id_col")
      id_col <- id_col[1]
    }
  }

  # Ensuring format
  if(is.matrix(eigengenes)) eigengenes <- as.data.frame(eigengenes)
  if(is.matrix(phenotypes)) phenotypes <- as.data.frame(phenotypes)

  # Looking for common id column if none specified (if none found, testing
  # rownames)
  if (is.null(id_col)) {
    matching_id_col <- colnames(eigengenes) %in% colnames(phenotypes)
    nb_matching <- length(matching_id_col[matching_id_col == TRUE])
    if (nb_matching > 1) {
      stop("More than one common column name between eigengenes and phenotypes.",
           "Therefore cannot match rows.")
    } else if (nb_matching == 0) {
      # Looking for matching rownames to use instead
      matching_rownames <- rownames(eigengenes) %in% rownames(phenotypes)
      if (all(matching_rownames)) {
        phenotypes <- phenotypes[match(rownames(eigengenes),
                                       rownames(phenotypes)), ]
      } else {
        warning("No common name found to be used as id, neither matching",
                " rownames. Using both dataframes as is for row matching.")
      }
    } else if (nb_matching == 1) {
      # Ordering rows using the id column and removing it after
      id_col <- colnames(eigengenes)[matching_id_col]
      phenotypes <- phenotypes[match(eigengenes[, id_col], phenotypes$sample), ]
      rownames(phenotypes) <- NULL
    }
  } else {
    if (length(id_col) == 2) {
      # Ordering rows using the id columns and removing them after
      if (id_col[1] %in% colnames(eigengenes)) {
        eigengenes <- eigengenes[order(eigengenes[, id_col[1]]), ]
        phenotypes <- phenotypes[order(phenotypes[, id_col[2]]), ]
        eigengenes[, id_col[1]] <- NULL
        phenotypes[, id_col[2]] <- NULL
      } else {
        eigengenes <- eigengenes[order(eigengenes[id_col[2]]), ]
        phenotypes <- phenotypes[order(phenotypes[id_col[1]]), ]
        eigengenes[, id_col[2]] <- NULL
        phenotypes[, id_col[1]] <- NULL
      }
    } else {
      # Ordering rows using the id column and removing it after
      id_col_sorted <- sort(eigengenes[, id_col])
      eigengenes[order(id_col_sorted),]
      phenotypes[order(id_col_sorted),]
      eigengenes[, id_col] <- NULL
      phenotypes[, id_col] <- NULL
    }
  }


  # Design matrix (dummy variable formation for qualitative variables)
  dummies_var <- lapply(colnames(phenotypes), function(dummy_name) {
    df <- phenotypes %>% dplyr::select(!!dummy_name)
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

  # Correlation + student test to assess correlation significance
  if (cor_func == "other") {
    asso_test_res <- list(cor = your_func(eigengenes, design_mat),
                          p = WGCNA::corPvalueStudent(asso, nrow(phenotypes)))
  } else {
    asso_test_res <- WGCNA::corAndPvalue(eigengenes, design_mat)
  }
  

  # TODO Check if a correction for multiple testing shouldn't be performed
  # here...

  return(list(
    association = asso_test_res$cor %>% as.data.frame,
    pval = asso_test_res$p %>% as.data.frame
  ))
}


# Removing errors about dplyr data-variables
utils::globalVariables(c("eigengene", "pval", "phenotype"))

#' Heatmap of modules phenotpic association
#'
#' Plot a heatmap of the correlation between all modules and the phenotypic
#' variables and the p value associated
#'
#' @param modules_phenotype list, data.frames of correlation and pvalue
#' associated
#' @param pvalue_th float, threshold in ]0;1[ under which module will be
#' considered as significantly associated
#' @param text_angle integer, angle in [0,360] of the x axis labels.
#' @param ... any other parameter you can provide to ggplot2::theme
#'
#' @importFrom ggplot2 ggplot geom_tile geom_point scale_color_gradient2
#' theme_bw xlab ylab theme
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select bind_cols
#' @importFrom stringr str_sort
#'
#' @return A ggplot object representing a heatmap with phenotype association
#' and related pvalues
#'
#' @examples
#' eigengene_mat <- data.frame(mod1 = rnorm(20, 0.1, 0.2),
#' mod2 = rnorm(20, 0.2, 0.2))
#' phenotype_mat <- data.frame(phenA = sample(c("X", "Y", "Z"), 20,
#'                             replace = TRUE),
#'                             phenB = sample(c("U", "V"), 20,
#'                             replace = TRUE),
#'                             stringsAsFactors = FALSE)
#' association <- associate_phenotype(eigengene_mat, phenotype_mat)
#' plot_modules_phenotype(association)
#'
#' @export

plot_modules_phenotype <- function(modules_phenotype, pvalue_th = 0.05, 
                                   text_angle = 90, ...){
  # Checks
  if (!is.list(modules_phenotype)) stop("modules_phenotype must be a list")
  if (!isTRUE(all.equal(names(modules_phenotype),
                        c("association", "pval")))) {
    stop("modules_penotype must have two elements: association and pvalue ",
         "associated")}
  if (!all(lapply(modules_phenotype, is.data.frame) %>% unlist)) {
    stop("All elements of modules_phenotype must be data.frame") }
  if (!isTRUE(all.equal(dim(modules_phenotype$association),
                        dim(modules_phenotype$pval)))) {
    stop("modules_phenotype data.frames must have the same number of rows") }
  if (!modules_phenotype %>% unlist %>% is.numeric) {
    stop("modules_phenotype data.frames must only contains numeric values")}
  if (!(isTRUE(all.equal(rownames(modules_phenotype$association),
                         rownames(modules_phenotype$pval))))) {
    stop("rownames of assocaition and pval must be the same") }
  if (!is.numeric(pvalue_th)) stop("pvalue_th must be a numeric value")
  if (length(pvalue_th) != 1) stop("pvalue_th must be a single value")
  if (pvalue_th <= 0 | pvalue_th >= 1) stop("pvalue_th must be in ]0;1[")
  if (length(text_angle) != 1 | !is.numeric(text_angle))
    stop("text_angle should be a single number")
  if (text_angle < -360 | text_angle > 360) 
    stop("text_angle should be a number between -360 and 360")

  # Data preparation
  df_cor <- modules_phenotype$association %>%
    tibble::rownames_to_column(var = "eigengene") %>%
    tidyr::pivot_longer(-eigengene, names_to = "phenotype", values_to = "cor")

  df_pval <- modules_phenotype$pval %>%
    tibble::rownames_to_column(var = "eigengene") %>%
    tidyr::pivot_longer(-eigengene, names_to = "phenotype", values_to = "pval")

  df_total <- dplyr::bind_cols(df_cor, df_pval %>% dplyr::select(pval)) %>%
    dplyr::mutate(signif = ifelse(pval > pvalue_th, FALSE, TRUE)) %>%
    dplyr::mutate(eigengene = stringr::str_sort(eigengene, numeric = TRUE)) %>%
    dplyr::mutate(eigengene = factor(eigengene, levels = unique(eigengene)))

  # Plotting
  suppressWarnings(quiet(
    # g <- ggplot2::ggplot(df_total, ggplot2::aes(x = factor(eigengene),
    #                                             y = factor(phenotype))) +
    g <- ggplot2::ggplot(df_total, ggplot2::aes(x = eigengene, y = phenotype)) +
      ggplot2::geom_tile(fill = "white") +
      ggplot2::geom_point(ggplot2::aes(colour = cor, size = signif)) +
      ggplot2::scale_color_gradient2() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = text_angle, 
                                                         hjust = 1),
                     ...) +
      ggplot2::xlab("Module") +
      ggplot2::ylab("Phenotype")
  ))
  # Warning message: Using size for a discrete variable is not advised.
  # (because of using TRUE/FALSE as size)
  # TODO : find why suppressWarnings + quiet doesn't make the function shut up
  g
}
