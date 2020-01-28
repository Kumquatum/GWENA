# library(dplyr)

df_expr <- list(df_microarray = kuehne_expr[, 5000:7000],
                df_rnaseq = gtex_expr[, 5000:7000])

# ==== full_pipeline ====

res_net <- build_net(df_expr$df_microarray)
res_detection <- detect_modules(df_expr$df_microarray, res_net$network)
res_enrich <- bio_enrich(res_detection$modules)

