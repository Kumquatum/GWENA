library(dplyr)

arg_cor_func <- c("pearson", "spearman", "bicor")
arg_network_type = c("unsigned", "signed", "signed hybrid")
fake_mat <-  runif(n = 10*45, max = 100) %>%  matrix(ncol = 45)
fake_cor_mat <- cor(fake_mat)
df_expr <- list(df_microarray = kuehne_expr[,1:100],
                df_rnaseq = gtex_expr[, 1:100])


# ==== cor_func_match ====

test_that("badly formatted input throw an error", {
  # cor_func should be one of the string allowed to design the correlation function
  lapply(arg_cor_func, function(arg){
    expect_error(cor_func_match(cor_func = arg), NA)
    expect_error(cor_func_match(cor_func = 42))
    expect_error(cor_func_match(cor_func = "this is not a function available"))
  })
})
test_that("output format is ok", {
  # Return a data.frame
  lapply(arg_cor_func, function(arg) {
    expect_true(is.function(cor_func_match(cor_func = arg)))
  })
})


# ==== get_fit.cor ====

test_that("cor_mat should be a matrix of correlations", {
  expect_error(get_fit.cor(cor_mat = cor(df_expr$df_microarray)), NA)
  # Type is matrix
  expect_error(get_fit.cor(cor_mat = 42))
  expect_error(get_fit.cor(cor_mat = "this is not a correlation matrix"))
  expect_error(get_fit.cor(cor_mat = fake_cor_mat %>% as.data.frame()))
  # Is a squared matrix
  expect_error(get_fit.cor(cor_mat = fake_cor_mat[, -1]))
  # Values in [0;1]
  expect_error(get_fit.cor(cor_mat = cbind(fake_cor_mat[, -1], 1:nrow(fake_cor_mat))))
  # Contains no NA (need to be checked with a symetric matrix)
  fake_cor_mat_na <- fake_cor_mat
  fake_cor_mat_na[c(2,4), c(2,4)] <- NA
  expect_warning(get_fit.cor(cor_mat = fake_cor_mat_na))
  # Follows a scale free topology
  expect_warning(get_fit.cor(cor_mat = fake_cor_mat))
})
test_that("fit_cut_off should be a number in ]0;1[", {
  expect_error(get_fit.cor(cor_mat = fake_cor_mat, fit_cut_off = -1))
  expect_error(get_fit.cor(cor_mat = fake_cor_mat, fit_cut_off = 42))
  expect_error(get_fit.cor(cor_mat = fake_cor_mat, fit_cut_off = c(0.1, 0.5, 0.9)))
  expect_error(get_fit.cor(cor_mat = fake_cor_mat, fit_cut_off = "this is not a numeric value"))
})
test_that("network_type should be one of the string allowed to design the correlation function", {
  lapply(arg_network_type, function(arg){
    # expect_error(get_fit.cor(cor_mat = cor(df_expr$df_microarray), network_type = arg), NA)
    expect_error(get_fit.cor(cor_mat = cor(df_expr$df_microarray), network_type = 42))
    expect_error(get_fit.cor(cor_mat = cor(df_expr$df_microarray), network_type = "this is not a network available"))
  })
})
test_that("output format is ok", {
  # Return a list with expected elements formated correctly
  res <- get_fit.cor(cor_mat = cor(df_expr$df_microarray))
  expect_true(is.list(res))
  expect_true(all(names(res) == c("power_value", "fit_table", "fit_above_cut_off", "metadata")))
  expect_true(all(names(res$metadata) == c("net_type")))
  expect_true(is.numeric(res$power_value) && length(res$power_value) == 1)
  expect_true(is.data.frame(res$fit_table) && colnames(res$fit_table) == c("Power", "SFT.R.sq", "slope", "truncated.R.sq", "mean.k.", "median.k.", "max.k."))
  expect_true
  expect_true(is.list(res$metadata) && all(names(res$metadata) == c("net_type")) && (res$metadata$net_type %in% arg_network_type))
})


# ==== get_fit.expr ====

test_that("data_expr should be a data.frame or a matrix of values >= 0", {
  expect_error(get_fit.expr(data_expr = cbind(df_expr$df_rnaseq, -10:(nrow(df_expr$df_rnaseq) - 11))))
})
test_that("data_expr should be either a matrix or data frame of intensities or counts", {
  expect_error(get_fit.expr(data_expr = df_expr$df_microarray), NA)
  expect_error(get_fit.expr(data_expr = df_expr$df_microarray %>% as.matrix), NA)
  expect_error(get_fit.expr(data_expr = "this is not a data.frame nor a matrix"))
  expect_error(get_fit.expr(data_expr = list("this is not a data.frame nor a matrix", "this shouldn't work")))
  expect_warning(get_fit.expr(data_expr = df_expr$df_microarray %>% t))
})
test_that("your_func is checked to returns correlation values ([-1;1])", {
  fake_bad_func <- function(df) { matrix(sample(-10:10, nrow(df)*ncol(df), replace = TRUE), ncol = ncol(df)) }
  lapply(df_expr, function(df) {
    expect_error(get_fit.expr(data_expr = df, cor_func = "other", your_func = fake_bad_func))
  })
})
test_that("output format is ok", {
  res <- get_fit.expr(data_expr = df_expr$df_microarray)
  expect_true(is.list(res))
  expect_true(all(names(res) == c("power_value", "fit_table", "fit_above_cut_off", "metadata")))
  expect_true(all(names(res$metadata) == c("net_type", "cor_func")))
  expect_true(is.numeric(res$power_value) && length(res$power_value) == 1)
  expect_true(is.data.frame(res$fit_table) && colnames(res$fit_table) == c("Power", "SFT.R.sq", "slope", "truncated.R.sq", "mean.k.", "median.k.", "max.k."))
  expect_true(is.list(res$metadata) && all(names(res$metadata) == c("net_type", "cor_func")) && (res$metadata$net_type %in% arg_network_type))
})


# ==== net_building ====
test_that("data_expr should be a data.frame or a matrix of values >= 0", {
  expect_error(net_building(data_expr = cbind(df_expr$df_rnaseq, -10:(nrow(df_expr$df_rnaseq) - 11))))
})
test_that("data_expr should be either a matrix or data frame of intensities or counts", {
  expect_error(net_building(data_expr = df_expr$df_microarray), NA)
  expect_error(net_building(data_expr = df_expr$df_microarray %>% as.matrix), NA)
  expect_error(net_building(data_expr = "this is not a data.frame nor a matrix"))
  expect_error(net_building(data_expr = list("this is not a data.frame nor a matrix", "this shouldn't work")))
  expect_warning(net_building(data_expr = df_expr$df_microarray %>% t))
})

