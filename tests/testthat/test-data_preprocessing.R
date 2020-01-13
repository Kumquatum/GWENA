library(dplyr)
library(magrittr)
# fake_microarray_data <- matrix(sample(dlnorm(seq(0,2, 0.05), sdlog = 0.9), 5*10, replace = TRUE), ncol = 5) %>% data.frame %>%
#   set_rownames(paste0("S", 1:10)) %>%
#   set_colnames(paste0("RNA", 1:5))

df_microarray <- kuehne_expr[,1:100]
df_rnaseq <- gtex_expr[, 1:100]
arg_type = c("mean", "median", "mad")
arg_method = c("at least one", "mean", "all")


# ==== filter_low_var ====

test_that("badly formatted input throw an error", {
  # data_expr should be a data.frame or a matrix
  expect_error(filter_low_var(data_expr = df_microarray), NA)
  expect_error(filter_low_var(data_expr = df_microarray %>% as.matrix), NA)
  expect_error(filter_low_var(data_expr = "this is not a data.frame nor a matrix"))
  expect_error(filter_low_var(data_expr = list("this is not a data.frame nor a matrix", "this shouldn't work")))
  # pct should be in ]0;1[
  expect_error(filter_low_var(df_microarray, pct = "this is not a numeric value"))
  expect_error(filter_low_var(df_microarray, pct = -1.5))
  expect_error(filter_low_var(df_microarray, pct = 42))
  expect_error(filter_low_var(df_microarray, pct = 0))
  expect_error(filter_low_var(df_microarray, pct = 1))
  # type should be a char from allowed ones
  expect_error(filter_low_var(df_microarray, type = 42))
  expect_error(filter_low_var(df_microarray, type = "this is not a function available"))
})

test_that("filtration has been done correctly", {
  # Not all columns have been removed
  lapply(arg_type, function(arg) {
    expect_gt(filter_low_var(df_microarray, type = arg) %>% ncol(), 0)
  })
  # Number of RNA transcript/probe have been reduced
  lapply(arg_type, function(arg) {
    expect_lt(filter_low_var(df_microarray, type = arg) %>% ncol(), ncol(df_microarray))
  })
  # No sample have been removed
  lapply(arg_type, function(arg) {
    expect_true(filter_low_var(df_microarray, type = arg) %>% nrow == nrow(df_microarray))
  })
})

test_that("output format is ok", {
  # Return a data.frame
  lapply(arg_type, function(arg) {
    expect_true(is.data.frame(filter_low_var(df_microarray, type = arg)))
  })
})

# ==== filter_RNA_seq ====

test_that("badly formatted input throw an error", {
  # data_expr should be a data.frame or a matrix
  expect_error(filter_RNA_seq(data_expr = df_rnaseq), NA)
  expect_error(filter_RNA_seq(data_expr = df_rnaseq %>% as.matrix), NA)
  expect_error(filter_RNA_seq(data_expr = "this is not a data.frame nor a matrix"))
  expect_error(filter_RNA_seq(data_expr = list("this is not a data.frame nor a matrix", "this shouldn't work")))
  # min_count should be > 1
  expect_error(filter_RNA_seq(df_rnaseq, min_count = "this is not a numeric value"))
  expect_error(filter_RNA_seq(df_rnaseq, min_count = -1.5))
  expect_error(filter_RNA_seq(df_rnaseq, min_count = 1))
  # method should be a char from allowed ones
  expect_error(filter_RNA_seq(df_rnaseq, method = 42))
  expect_error(filter_RNA_seq(df_rnaseq, method = "this is not a function available"))
})


test_that("filtration has been done correctly", {
  # Not all columns have been removed
  lapply(arg_method, function(arg) {
    expect_gt(filter_RNA_seq(df_rnaseq, method = arg) %>% ncol(), 0)
  })
  # Number of RNA transcript/probe have been reduced
  lapply(arg_method, function(arg) {
    expect_lt(filter_RNA_seq(df_rnaseq, method = arg) %>% ncol(), ncol(df_rnaseq))
  })
  # No sample have been removed
  lapply(arg_method, function(arg) {
    expect_true(filter_RNA_seq(df_rnaseq, method = arg) %>% nrow == nrow(df_rnaseq))
  })
})

test_that("output format is ok", {
  # Return a data.frame
  lapply(arg_method, function(arg) {
    expect_true(is.data.frame(filter_RNA_seq(df_rnaseq, method = arg)))
  })
})
