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
joined <- join_gost(list(classic_gost, custom_gost_symbols))

test_that("input is a gost result", {
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
  expect_warning(join_gost(list(classic_gost, custom_gost_entrez)))
})

# test_that("correct input return gost_result", {
#   expect_true(all(names(joined) %in% c("result", "meta")))
# })

# test_that("warnings about ")