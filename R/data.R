#' Transcriptomic data from the Kuehne et al. publication
#'
#' A dataset containing the expression levels collapsed to the gene level.
#'
#' @format A data frame with 48 rows and 20672 columns :
#' \describe{
#'   \item{rows}{samples, reference combining the microarray slide, the condition, and the array}
#'   \item{columns}{genes, IDs are gene names (also called gene symbols)}
#' }
#'
#' @source \url{https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3547-3}
"kuehne_expr"


#' Traits data linked to samples in transcriptomic data from the Kuehne et al. publication
#'
#' A dataset containing the phenotype of the donors and technical information about the experiment
#'
#' @format A data frame with 48 rows and 5 columns :
#' \describe{
#'   \item{Slide}{Reference number of the microarray's slide.}
#'   \item{Array}{Array number, 8 by slide usually}
#'   \item{Exp}{Experiment number}
#'   \item{Condition}{Either old (between 55 and 66 years old) or young (between 20 to 25 years old)}
#'   \item{Age}{Real age of the donor}
#' }
#'
#' @source \url{https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3547-3}
"kuehne_traits"

