library(dplyr)
library(data.table)
library(usethis)

bash_cmds <- "wget -qO- https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt > gtex_phenotype.txt"
system(bash_cmds)

rownames_gtex_expr <- rownames(gtex_expr) %>% sort()

gtex_traits <- data.table::fread("gtex_phenotype.txt") %>%
  dplyr::filter(pmatch(SUBJID, rownames_gtex_expr) > 0) %>%
  dplyr::arrange(SUBJID) %>%
  dplyr::mutate(SUBJID = rownames_gtex_expr)

usethis::use_data(gtex_traits, overwrite = TRUE)
