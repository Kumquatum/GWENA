# ===== setup ====

df_expr <- list(df_microarray = kuehne_expr[, 5000:7000],
                df_rnaseq = gtex_expr[, 5000:7000])

# Checking if gprofiler is up. If not, skipping tests.
# Note : couldn't check API for gost itself since no GET endpoint with pin method set by authors
safe_get <- purrr::safely(httr::GET)
i <- 1
while (i <= 4) {
  res_get <- safe_get("https://biit.cs.ut.ee/gprofiler/", httr::timeout(i*10))
  if (!is.null(res_get$result)) { i <- 3 } else { i <- i + 1 }
}
if (!is.null(res_get$result)) {is_gprofiler_down <- FALSE } else { is_gprofiler_down <- TRUE }


# ==== pipeline ====

res_net <- build_net(df_expr$df_microarray)
res_detection <- detect_modules(df_expr$df_microarray, res_net$network)
if (!is_gprofiler_down) { res_enrich <- bio_enrich(res_detection$modules) }


