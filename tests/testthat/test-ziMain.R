context("Test ziMain function")


test_that("ZiMain function works", {
  data(mtx)
  Zi <- ziMain(mtx)
  expect_equal(sum(is.na(Zi@output)), 657)
  Zi_poisson <- ziMain(mtx, dist = "poisson")
  expect_equal(sum(is.na(Zi_poisson@output)), 676)
})



