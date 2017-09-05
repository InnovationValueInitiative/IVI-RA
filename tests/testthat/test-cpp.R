context("C++ unit tests")

test_that("sim_acr_test", { 
  expect_equal(iviRA:::sim_acr_test(), 0)
})
test_that("sim_lm_test", { 
  expect_equal(iviRA:::sim_lm_test(), 0)
})