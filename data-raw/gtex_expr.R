library(dplyr)
library(tibble)
library(data.table)
library(usethis)

# Keeping only 15 000 genes for demo purposes. If you wish to keep everythong, remove "| head -n 15003 "
bash_cmds <- "wget -qO- https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz | gunzip | head -n 15003 > gtex_RNAseq_trunc.gct;
wget -qO gtex_sample_attributes.txt https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
system(bash_cmds)

gtex_muscle_sample <- data.table::fread("gtex_sample_attributes.txt") %>%
  dplyr::select(SAMPID, SMTSD) %>%
  dplyr::filter(SMTSD == "Muscle - Skeletal")

gtex_expr <- data.table::fread("gtex_RNAseq_trunc.gct") %>%
  tibble::column_to_rownames("Name") %>%
  dplyr::select(-Description) %>%
  t %>%
  as.data.frame %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::filter(sample %in% gtex_muscle_sample$SAMPID) %>% # keeping only muscle data
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("sample") %>%
  .[1:50, ] # keeping only 50 samples


usethis::use_data(gtex_expr, overwrite = TRUE)
