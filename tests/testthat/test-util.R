context("util.R unit tests")
library("data.table")

# sim_survtime -----------------------------------------------------------------
test_that("sim_survtime", {
  expect_error(sim_survtime(n = 1), NA)
  expect_error(sim_survtime(n = 2), NA)
  sim <- sim_survtime(n = 1000, male = 0, haq0 = 0, cycle_length = 6,
                      logor = 0, loghr = rep(0, 5), age = 55)
  expect_equal(mean(sim) - 55, iviRA::lifetable.female[age == 55, ex],
               tolerance = 1, scale = 1)
})
