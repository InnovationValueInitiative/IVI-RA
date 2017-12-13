context("pop.R")
library("data.table")

# sim_iviRA -------------------------------------------------------------------
test_that("sample_pop", {
  expect_error(sample_pop(age = 1000))
  
  # heterogeneous population
  expect_error(sample_pop(type = "heterog"), NA)
  pop <- sample_pop(n = 5)
  expect_true(inherits(pop, "matrix"))
})

test_that("get_input_data", {
  pop <- sample_pop(n = 5)
  expect_true(inherits(get_input_data(pop), "input_data"))
  expect_error(get_input_data(pop = list(x = 1)))
  
  # checking error messages given variables in pop
  pop.names <- c("age", "male", "weight", "prev_dmards", "das28", "sdai",
                 "cdai", "haq0")
  for (n in pop.names){
    expect_error(get_input_data(pop[, -which(colnames(pop) == n)]))
  }
  
  # design matrcices
  test_dm <- function(var_arg, var_ret){
    n <- 3
    pop <- sample_pop(n = n)
    args <- as.list(formals(get_input_data))
    args$pop <- pop
    
    # correct matric passed
    args[[var_arg]] <- matrix(rep(3, n)) 
    input.data <- do.call("get_input_data", args)
    expect_equal(input.data[[var_ret]], args[[var_arg]])
    
    # incorrect matrix passed
    args[[var_arg]] <- matrix(rep(1, n - 1)) 
    expect_error(do.call("get_input_data", args))
  }
  test_dm("x_acr", "x.acr")
  test_dm("x_haq", "x.haq")
  test_dm("x_ttd_all", "x.ttd.all")
  test_dm("x_ttd_da", "x.ttd.da")
  test_dm("x_ttd_eular", "x.ttd.eular")
  test_dm("x_mort", "x.mort")

})

