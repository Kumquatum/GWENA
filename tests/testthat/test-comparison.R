library(magrittr)
library(NetRep)

# Datasets
data("NetRep")
data_list <- list(cohort1 = discovery_data,
                  cohort2 = test_data,
                  cohort3 = test_data*0.9) %>% lapply(abs)
data_list_se <- lapply(data_list, function(data_expr) {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = t(data_expr)),
    colData = S4Vectors::DataFrame(kuehne_traits[1:30,])
  )
})
correlation_list <- list(cohort1 = discovery_correlation,
                         cohort2 = test_correlation,
                         cohort3 = test_correlation*0.95)
adja_list <- list(cohort1 = discovery_network,
                  cohort2 = test_network,
                  cohort3 = test_network*0.97)

# Different modules conformations
module_labels_switched1 <- module_labels %>%
  replace(c(which(. == 4), which(. == 3)), .[c(which(. == 3), which(. == 4))])
module_labels_switched2 <- module_labels %>%
  replace(c(which(. == 2), which(. == 3)), .[c(which(. == 3), which(. == 2))])
mod_labels_single_cond <- module_labels %>%
  split(module_labels) %>%
  lapply(names)
mod_labels_multi_cond <- list(cohort1 = module_labels,
                              cohort2 = module_labels_switched1,
                              cohort3 = module_labels_switched2) %>%
  lapply(function(cond){cond %>% split(cond) %>% lapply(names)})



# ==== z_summary ====
mini_comp <- compare_conditions(data_list, adja_list, correlation_list,
                               mod_labels_single_cond, ref = "cohort1",
                               n_perm = 10)
test_that("Valid input doesn't thow errors", {
  expect_error(z_summary(mini_comp$result$cohort1$cohort2$observed,
                         mini_comp$result$cohort1$cohort2$nulls),
               NA)
})
test_that("Check on input format", {
  expect_error(z_summary(NULL, mini_comp$result$cohort1$cohort2$nulls))
  expect_error(z_summary("a", mini_comp$result$cohort1$cohort2$nulls))
  expect_error(z_summary(1, mini_comp$result$cohort1$cohort2$nulls))
  expect_error(z_summary(letters[1:5], mini_comp$result$cohort1$cohort2$nulls))
  expect_error(z_summary(1:5, mini_comp$result$cohort1$cohort2$nulls))
  expect_error(z_summary(matrix(letters[1:9], 3),
                         mini_comp$result$cohort1$cohort2$nulls))
  expect_error(z_summary(mini_comp$result$cohort1$cohort2$observed, NULL))
  expect_error(z_summary(mini_comp$result$cohort1$cohort2$observed, "a"))
  expect_error(z_summary(mini_comp$result$cohort1$cohort2$observed, 1))
  expect_error(z_summary(mini_comp$result$cohort1$cohort2$observed, letters[1:5]))
  expect_error(z_summary(mini_comp$result$cohort1$cohort2$observed, 1:5))
  expect_error(z_summary(mini_comp$result$cohort1$cohort2$observed,
                         matrix(1:9, 3)))
  expect_error(z_summary(mini_comp$result$cohort1$cohort2$observed,
                         array(letters[1:3], c(1:3))))
})

# ==== compare_conditions ====
test_that("Valid input doesn't thow errors", {
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond,
                                  ref = "cohort1"), NA) # single ref // single mod in vector
  expect_error(compare_conditions(data_list, adja_list, correlation_list, list(cohort1 = mod_labels_single_cond),
                                  ref = "cohort1"), NA) # single ref // single mod in list
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_multi_cond,
                                  ref = "cohort1"), NA) # single ref // multi mod
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_multi_cond,
                                  ref = c("cohort1", "cohort2")), NA) # multi ref // multi mod
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_multi_cond,
                                  ref = "cross comparison"), NA) # cross comparison
  expect_error(compare_conditions(data_list, adja_list, NULL, mod_labels_single_cond,
                                  cor_func = "other", your_func = cor), NA)
  expect_error(compare_conditions(data_list_se, adja_list, correlation_list, mod_labels_single_cond,
                                  ref = "cohort1"), NA) # single ref // single mod in vector // with SummarizedExperiment
})

test_that("Checks on compatibility conditions/ref/modules_list are correctly done", {
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond,
                                  ref = c("cohort1", "cohort2")))  # multi ref // single mod) in vector
  expect_error(compare_conditions(data_list, adja_list, correlation_list, list(cohort1 = mod_labels_single_cond),
                                  ref = c("cohort1", "cohort2")))  # multi ref // single mod) in list
  expect_warning(compare_conditions(data_list[c("cohort1", "cohort2")], adja_list[c("cohort1", "cohort2")], correlation_list[c("cohort1", "cohort2")],
                                    mod_labels_multi_cond, ref = "cohort1"))  # Warning about cohort2 et cohort3 defined in modules but not elsewhere
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_multi_cond[c("cohort1", "cohort2")],
                                  ref = "cross comparison"))
})


test_that("Checks on missing values", {
  expect_error(compare_conditions())
  expect_error(compare_conditions(data_list))
  expect_error(compare_conditions(data_list, adja_list))
  expect_error(compare_conditions(data_list, adja_list, correlation_list))
})

test_that("Checks on input format", {
  # HF_remove_limit_neg_value_build_net : removing negative value prohibition
  # data_list_negative_values <- list(cohort1 = discovery_data, cohort2 = test_data, cohort3 = test_data*0.9)
  # expect_error(compare_conditions(data_list_negative_values, adja_list, correlation_list, mod_labels_single_cond))
  net_list_not_sq <- adja_list %>% lapply(function(x) { x[,-1] })
  expect_error(compare_conditions(data_list, net_list_not_sq, correlation_list, mod_labels_single_cond))
  cor_list_wrong_range <- correlation_list %>% lapply(function(x) { x + sample(c(-2,2), 1) })
  expect_error(compare_conditions(data_list, adja_list, cor_list_wrong_range, mod_labels_single_cond))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, ref = NULL))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, ref = 42))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, ref = "inexisting condittion"))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, ref = c(42, "inexisting condition")))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, ref = c("cohort1", "inexisting condition")))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, test = 42))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, test = "inexisting condition"))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, test = c(42, "inexisting condition")))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, mod_labels_single_cond, test = c("cohort1", "inexisting condition")))
  expect_error(compare_conditions(data_list, adja_list, NULL, mod_labels_single_cond, cor_func = 42))
  expect_error(compare_conditions(data_list, adja_list, NULL, mod_labels_single_cond, cor_func = c(42, "not a func")))
  expect_error(compare_conditions(data_list, adja_list, NULL, mod_labels_single_cond, cor_func = "not a func"))
  expect_error(compare_conditions(data_list, adja_list, NULL, mod_labels_single_cond, cor_func = "other"))
  expect_error(compare_conditions(data_list, adja_list, NULL, mod_labels_single_cond, cor_func = "other", your_func = function(x){ cor(x) + 3 }))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, "this is not a module"))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, c(cohort1 = "this is not a module", cohort2 = "this is not either")))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, list(cohort1 = "this is not a module", cohort2 = "this is not either")))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, list(cohort1 = mod_labels_single_cond, cohort2 = "this is not either")))
  expect_error(compare_conditions(data_list, adja_list, correlation_list, list(cohort1 = mod_labels_single_cond,
                                                                                  cohort2 = c(`0` = "nodex", `1` = "nodey", `2` = "nodez"))))
})

