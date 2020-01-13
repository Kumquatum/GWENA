#' Transcriptomic data from the Kuehne et al. publication
#'
#' A dataset containing the expression levels collapsed to the gene level. Obtained from script provided in additional data n°10 runned on GSE85358.
#'
#' @format A data frame with 48 rows (samples) and 20672 columns (probes).
#'
#' @source \url{https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3547-3}
"kuehne_expr"


#' Traits data linked to samples in transcriptomic data from the Kuehne et al. publication
#'
#' A dataset containing the phenotype of the donors and technical information about the experiment
#'
#' @format A data frame with 48 rows (samples) and 5 columns :
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

#' Transcriptomic muscle data from GTEx consorsium RNA-seq data
#'
#' A subset of GTEx RNA-seq dataset containing read counts collapsed to gene level. Obtained from the following script :
#' library(dplyr)
#'
#' # Keeping only 15 000 genes for demo purposes. If you wish to keep everythong, remove "| head -n 15003 "
#' bash_cmds <- "wget -qO- https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz | gunzip | head -n 15003 > gtex_RNAseq_trunc.gct;
#' wget -qO gtex_sample_attributes.txt https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
#' system(bash_cmds)
#'
#' gtex_muscle_sample <- data.table::fread("gtex_sample_attributes.txt") %>%
#'   select(SAMPID, SMTSD) %>% filter(SMTSD == "Muscle - Skeletal")
#'
#' gtex_expr <- data.table::fread("gtex_RNAseq_trunc.gct") %>%
#'   tibble::column_to_rownames("Name") %>% select(-Description) %>% t %>% as.data.frame %>% # formatting
#'   tibble::rownames_to_column("sample") %>% filter(sample %in% gtex_muscle_sample$SAMPID) %>% # keeping only muscle data
#'   tibble::remove_rownames() %>% tibble::column_to_rownames("sample") %>% .[1:50, ] # formatting and keeping only 50 samples
#'
#'
#' @format A data frame with XX TODO rows (samples) and YY TODO columns (genes)
#'
#' @source \url{https://gtexportal.org/home/datasets}
"gtex_expr"

#' Traits data linked to samples in transcriptomic data from GTEx
#'
#' A dataset containing phenotypes of donors. From public data.
#' Note: protected data contain more information but require dbGap accessh (see https://gtexportal.org/home/protectedDataAccess).
#'
#' @format A data frame with 48 rows (samples) and 5 columns :
#' \describe{
#'   \item{SUBJID}{Subject ID, GTEx Public Donor ID}
#'   \item{SEX}{Sex, donor's Identification of sex based upon self-report : 1=Male, 2=Female}
#'   \item{AGE}{Age range, elapsed time since birth in years}
#'   \item{DTHHRDY}{Hardy Scale : 0=Ventilator Case, 1=Violent and fast death, 2=Fast death of natural causes, 3=Intermediate death, 4=Slow death)}
#' }
#'
#' @source \url{https://gtexportal.org/home/datasets}
"gtex_traits"
