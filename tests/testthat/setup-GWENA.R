# ===== setup ====

#### Expression data.frames
df_expr <- list(df_microarray = kuehne_expr[, 3000:5000],
# df_expr <- list(df_microarray = kuehne_expr[, 5000:7000],
                df_rnaseq = gtex_expr[, 5000:7000])

#### SummarizedExperiment version of kuehne data
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(expr = t(df_expr$df_microarray)),
  colData = S4Vectors::DataFrame(kuehne_traits)
)
S4Vectors::metadata(se) <- list(
  GEO_accession_id = "GSE85358",
  URL = "https://www.ncbi.nlm.nih.gov/pubmed/28201987"
)

#### Checking if gprofiler is up. If not, skipping tests.
# Note : couldn't check properly API for gost since no GET method set by
# authors, so checking through a GET query to the POST method and using error
# code as check.
safe_get <- purrr::safely(httr::GET)
i <- 1
while (i <= 4) {
  res_get <- safe_get("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
                      httr::timeout(i*10))
  if (!is.null(res_get$result)) { i <- 5 } else { i <- i + 1 }
}
if (!is.null(res_get$result)) {is_gprofiler_down <- FALSE
} else { is_gprofiler_down <- TRUE }


# ==== pipeline steps needed to avoid redundant operations ====

res_net <- build_net(df_expr$df_microarray, n_threads = 1)
res_detection <- detect_modules(df_expr$df_microarray, res_net$network)
if (!is_gprofiler_down) { res_enrich <- bio_enrich(res_detection$modules) }


