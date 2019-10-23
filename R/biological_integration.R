#' Modules enrichment
#'
#' Enrich genes list from modules.
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
#' TODO
#' @importFrom gprofiler2 gost
#'
#' @export

modules_enrichment <- function(modules, custom_gmt = NULL, ...) {
  # Checks

  # Splitting genes by module
  modules_list <- base::split(names(modules), modules)

  # Enrichment
  enriched_modules <- lapply(modules_list, function(module){
    gprofiler2::gost(query = module, ...)
  })
  # TODO Check why multiple query return less enrichment results than mono query on each module

  return(enriched_modules)
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
#' TODO
#' @importFrom WGCNA corPvalueStudent
#' @importFrom dplyr select
#'
#' @export

modules_phenotype <- function(eigengenes, phenotypes) {
  # Checks

  # Design matrix (dummy variable formation for qualitative variables)
  dummies_var <- lapply(colnames(phenotypes), function(dummy_name) {
    df <- phenotypes %>% dplyr::select(dummy_name)
    if (!is.numeric(df[1,])){
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
  design_mat <- base::as.data.frame(dummies_var, check.names = FALSE)

  # Correlations
  association <- stats::cor(eigengenes, design_mat)

  # P value associated
  pval_association <- WGCNA::corPvalueStudent(association, nrow(phenotypes))

  # TODO Check if a correction for multiple testing shouldn't be performed here...

  return(list(
    association = association %>% base::as.data.frame,
    pval = pval_association %>% base::as.data.frame
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
#' TODO
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
