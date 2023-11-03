# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(zitools)


data(mtx)
Zi <- ziMain(mtx)


## Unit tests for comparing the case where all weights are 1 to the case where zi@inputcounts is used
# Zi object with all weights equal to 1
ZiW1 <- Zi
ZiW1@weights <- Zi@weights * 0 + 1

test_that(
  "Compare results for sd, mean, rowMeans2 and colMeans2 on Zi object with weights equal to 1 to the case where Zi@inputcounts is used",
  {
    expect_equal(sd(ZiW1), sd(ZiW1@inputcounts))
    expect_equal(mean(ZiW1), mean(ZiW1@inputcounts))
    expect_equal(unname(rowMeans2(ZiW1)), unname(rowMeans2(ZiW1@inputcounts)))
    expect_equal(unname(colMeans2(ZiW1)), unname(colMeans2(ZiW1@inputcounts)))
  }
)




## Testing whether weights = 1/2 lead to correct results

# Zi object with nrows*2 but weights*(1/2)
mtx2 <- rbind(mtx, mtx)
Zi2 <- ziMain(mtx2)
Zi2@weights <- rbind(Zi@weights, Zi@weights)/2


test_that(
  "Testing whether results for colMeans2, mean and sd on small dataset is the same as for datasets with nrows*2 but weights*(1/2)  ",
  {
    expect_equal(colMeans2(Zi2), colMeans2(Zi))
    expect_equal(mean(Zi2), mean(Zi)) #Equal hat eine gewisse Toleranz
    expect_equal(sd(Zi2), sd(Zi))
  }
)



## Testing if in case of all counts >0 all weights are set to 1

# Input data with all values>0
mtxAll <- mtx
mtxAll[mtxAll==0] <- 1
ZiAll <- ziMain(mtxAll)

test_that("Testing if weights of Zi object are 1",{
  expect_equal(sum(ZiAll@weights<1), 0)  #Wenn man expect_identical nimmt, dann failed der Test, da actual ein integer und expected ein double
})




