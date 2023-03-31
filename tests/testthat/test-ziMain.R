context("Test ziMain function")


test_that("ZiMain function works", {
  n <- 1000
  male <- sample(c(0,1), size = n, replace = TRUE)
  z <- rbinom(n = n, size = 1, prob = 0.3)
  mean(z == 0)
  y_sim <- ifelse(z == 0, 0,
                  rnbinom(n = n,
                          mu = exp(1.3 + 1.5 * (male == 1)),
                          size = 2))
  mtx <- matrix(y_sim, 100, 10)
  Zi <- ziMain(mtx)
  expect_equal(sum(is.na(Zi@output)), 634)
  Zi_poisson <- ziMain(mtx, dist = "poisson")
  expect_equal(sum(is.na(Zi_poisson@output)), 660)
})

