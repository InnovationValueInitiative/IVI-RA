context("Model structures")

# select_model_structure ------------------------------------------------------
test_that("select_model_structures", {
  
  # default options are constant when selected options are greater than length 1
  expect_error(
    select_model_structures(tx_ihaq = c("acr-eular-haq", "acr-haq"),
                            tx_iswitch = c("acr-eular-switch", "acr-das28-switch")), NA
  )
  
  expect_warning(
    select_model_structures(tx_ihaq = c("acr-eular-haq", "acr-haq"),
                            tx_iswitch = c("acr-eular-switch", "acr-das28-switch")), NA
  )
  
  # return error when options are of different length
  expect_error(
    select_model_structures(tx_ihaq = c("acr-eular-haq", "acr-haq", "haq"),
                            tx_iswitch = c("acr-eular-switch", "acr-das28-switch"))
  )
  
  # combinations of model structures must be correct
  expect_error(
    select_model_structures(tx_ihaq = c("acr-eular-haq", "acr-haq"),
                            tx_iswitch = c("acr-eular-switch", "acr-eular-switch"))
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                            tx_iswitch = "acr-switch")
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                            tx_iswitch = "acr-das28-switch")
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                            tx_iswitch = "acr-sdai-switch")
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                            tx_iswitch = "acr-sdai-switch")
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                            tx_iswitch = "acr-eular-switch")
  )
  
})