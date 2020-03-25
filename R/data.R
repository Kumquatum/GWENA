#' Transcriptomic data from the Kuehne et al. publication
#'
#' A dataset containing the expression levels collapsed to the gene level.
#' Obtained from script provided in additional data n°10 runned on GSE85358 and
#' reduced from probe to gene by WGCNA::collapseRows with median as fucntion.
#'
#' @format A data frame with 48 rows (samples) and 15801 columns (genes).
#'
#' @source \url{https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3547-3}
"kuehne_expr"


#' Traits data linked to samples in transcriptomic data from the Kuehne et al. publication
#'
#' A dataset containing the phenotype of the donors and technical information about the
#' experiment
#'
#' @format A data frame with 48 rows (samples) and 5 columns :
#' \describe{
#'   \item{Slide}{Reference number of the microarray's slide.}
#'   \item{Array}{Array number, 8 by slide usually}
#'   \item{Exp}{Experiment number}
#'   \item{Condition}{Either old (between 55 and 66 years old) or young (between 20 to 25
#'   years old)}
#'   \item{Age}{Real age of the donor}
#' }
#'
#' @source \url{https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3547-3}
"kuehne_traits"

#' Transcriptomic muscle data from GTEx consorsium RNA-seq data
#'
#' A subset of GTEx RNA-seq dataset containing read counts collapsed to gene level.
#'
#' @format A data frame with 50 rows (samples) and 15000 columns (genes)
#'
#' @source \url{https://gtexportal.org/home/datasets}
"gtex_expr"

#' Traits data linked to samples in transcriptomic data from GTEx
#'
#' A dataset containing phenotypes of donors. From public data.
#' Note: protected data contain more information but require dbGap accessh
#' (see https://gtexportal.org/home/protectedDataAccess).
#'
#' @format A data frame with 50 rows (samples) and 4 columns :
#' \describe{
#'   \item{SUBJID}{Subject ID, GTEx Public Donor ID}
#'   \item{SEX}{Sex, donor's Identification of sex based upon self-report : 1=Male,
#'   2=Female}
#'   \item{AGE}{Age range, elapsed time since birth in years}
#'   \item{DTHHRDY}{Hardy Scale : 0=Ventilator Case, 1=Violent and fast death, 2=Fast
#'    death of natural causes, 3=Intermediate death, 4=Slow death)}
#' }
#'
#' @source \url{https://gtexportal.org/home/datasets}
"gtex_traits"
