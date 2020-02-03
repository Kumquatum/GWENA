library(dplyr)
library(magrittr)

query <- res_detection$modules[[5]]
query_entrez <- gprofiler2::gconvert(query, target = "ENTREZGENE_ACC")$target

classic_gost <- gprofiler2::gost(query)

gmt_entrez_path <- system.file("extdata", "h.all.v6.2.entrez.gmt", package = "GWENA", mustWork = TRUE)
gmt_entrez_id <- gprofiler2::upload_GMT_file(gmt_entrez_path)
custom_gost_entrez <- gprofiler2::gost(query_entrez, organism = gmt_entrez_id)

gmt_symbols_path <- system.file("extdata", "h.all.v6.2.symbols.gmt", package = "GWENA", mustWork = TRUE)
gmt_symbols_id <- gprofiler2::upload_GMT_file(gmt_symbols_path)
custom_gost_symbols <- gprofiler2::gost(query, organism = gmt_symbols_id)


# ==== join_gost ====

joined_gost <- join_gost(list(classic_gost, custom_gost_symbols))

test_that("input is a gost result", {
  expect_error(join_gost())
  expect_error(join_gost("this is not a gost result"))
  expect_error(join_gost(42))
  expect_error(join_gost(1:42))
  expect_error(join_gost(list(classic_gost)))
  expect_error(join_gost(list(classic_gost, NULL)))
  expect_error(join_gost(list(a = 1:5, b = "this is not a list of gost results")))
  expect_error(join_gost(list(fake_gost1 = list(result = "this is not a true result",
                                                meta = list(query_metadata = "this",
                                                            result_metadata = "is",
                                                            genes_metadata = "not",
                                                            timestamp = "a",
                                                            version = "true meta")),
                              fake_gost2 = list(result = "this is not a true result",
                                                meta = list(query_metadata = "this",
                                                            result_metadata = "is",
                                                            genes_metadata = "not",
                                                            timestamp = "a",
                                                            version = "true meta")))))
  expect_error(join_gost(list(classic_gost, list(result = "this is not a true result",
                                                 meta = list(query_metadata = "this",
                                                             result_metadata = "is",
                                                             genes_metadata = "not",
                                                             timestamp = "a",
                                                             version = "true meta")))))
})

test_that("gost objects in list are compatible", {
  expect_error(join_gost(list(classic_gost, custom_gost_entrez))) # not same length
  expect_warning(join_gost(list(classic_gost, gprofiler2::gost(query_entrez[28:140], organism = gmt_entrez_id)))) # not same id type (27:140 arbitrairy to get gost result)
  mock_custom_gost <- custom_gost_symbols
  mock_custom_gost$meta$query_metadata$ordered <- TRUE
  expect_warning(join_gost(list(classic_gost, mock_custom_gost))) # element different
})

test_that("return a gost object", {
  expect_false(is.null(joined_gost))
  expect_false(!all(names(joined_gost) %in% c("result", "meta")))
  expect_false(!is.data.frame(joined_gost$result))
  expect_false(any(is.na(match(c("query", "significant", "p_value", "term_size", "query_size", "intersection_size", "precision", "recall",
                        "term_id", "source", "term_name", "effective_domain_size", "source_order", "parents"), colnames(joined_gost$result)))))
  expect_false(!is.list(joined_gost$meta))
  expect_false(any(is.na(match(c("query_metadata", "result_metadata", "genes_metadata", "timestamp", "version"), names(joined_gost$meta)))))
})


# ==== bio_enrich ====

test_that("module input is correctly checked", {
  expect_error(bio_enrich())
  expect_error(bio_enrich(NULL))
  expect_error(bio_enrich(42))
  expect_error(bio_enrich(1:42))
  expect_error(bio_enrich(matrix(1:9, 3)))
  expect_error(bio_enrich(data.frame(a = c(letters[1:5], b = letters[6:10]))))
  expect_error(bio_enrich(list(c(res_detection$modules[4:5], c = 1:5))))
  expect_warning(bio_enrich("this is not modules"))
})

test_that("custom_gmt input is correctly checked", {
  expect_error(bio_enrich(res_detection$modules[[5]], 42))
  expect_error(bio_enrich(res_detection$modules[[5]], "this is not a path"))
  expect_error(bio_enrich(res_detection$modules[[5]], 1:42))
  expect_error(bio_enrich(res_detection$modules[[5]], c("~/.bashrc", "~/.bash_history")))
})

test_that("returns enriched modules", {
  expect_false(is.null(res_enrich))
  expect_false(!all(names(res_enrich) %in% c("result", "meta")))
  expect_false(!is.data.frame(res_enrich$result))
  expect_false(any(is.na(match(c("query", "significant", "p_value", "term_size", "query_size", "intersection_size", "precision", "recall",
                                 "term_id", "source", "term_name", "effective_domain_size", "source_order", "parents"), colnames(res_enrich$result)))))
  expect_false(!is.list(res_enrich$meta))
  expect_false(any(is.na(match(c("query_metadata", "result_metadata", "genes_metadata", "timestamp", "version"), names(res_enrich$meta)))))
})


# ==== plot_enrichment ====

test_that("input enrich_output is correctly checked", {
  expect_error(plot_enrichment())
  expect_error(plot_enrichment("this is not a enrich_output result"))
  expect_error(plot_enrichment(42))
  expect_error(plot_enrichment(1:42))
  expect_error(plot_enrichment(list(a = 1:5, b = "this is not a list of enrich_output result")))
  expect_error(plot_enrichment(list(res_enrich[[5]], NULL)))
  expect_error(plot_enrichment(fake_enrich_output = list(result = "this is not a true result",
                                                         meta = list(query_metadata = "this",
                                                         result_metadata = "is",
                                                         genes_metadata = "not",
                                                         timestamp = "a",
                                                         version = "true meta"))))
})

